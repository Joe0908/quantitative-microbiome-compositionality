"""Small, streaming XLSX reader for very large supplementary workbooks.

Only cached cell values are required here. Avoiding openpyxl's full style model
cuts peak memory substantially for the heavily formatted MetaCardis workbook.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CELL_REFERENCE = re.compile(r"([A-Z]+)")


def _column_index(reference: str) -> int:
    match = CELL_REFERENCE.match(reference)
    if not match:
        raise ValueError(f"Invalid Excel cell reference: {reference}")
    result = 0
    for character in match.group(1):
        result = result * 26 + ord(character) - 64
    return result - 1


class LiteXlsx:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.archive = zipfile.ZipFile(self.path)
        self.sheet_paths = self._sheet_paths()
        self.shared_strings = self._shared_strings()

    def close(self) -> None:
        self.archive.close()

    def __enter__(self) -> "LiteXlsx":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _sheet_paths(self) -> dict[str, str]:
        workbook = ET.fromstring(self.archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(
            self.archive.read("xl/_rels/workbook.xml.rels")
        )
        targets = {
            relationship.attrib["Id"]: relationship.attrib["Target"]
            for relationship in relationships.findall(f"{{{PKG_REL_NS}}}Relationship")
        }
        result: dict[str, str] = {}
        for sheet in workbook.findall(f".//{{{MAIN_NS}}}sheet"):
            relationship_id = sheet.attrib[f"{{{REL_NS}}}id"]
            target = targets[relationship_id].lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            result[sheet.attrib["name"]] = target
        return result

    def _shared_strings(self) -> list[str]:
        if "xl/sharedStrings.xml" not in self.archive.namelist():
            return []
        strings: list[str] = []
        with self.archive.open("xl/sharedStrings.xml") as handle:
            for event, element in ET.iterparse(handle, events=("end",)):
                if element.tag == f"{{{MAIN_NS}}}si":
                    text = "".join(
                        node.text or ""
                        for node in element.iter()
                        if node.tag == f"{{{MAIN_NS}}}t"
                    )
                    strings.append(text)
                    element.clear()
        return strings

    def _cell_value(self, cell: ET.Element) -> object:
        cell_type = cell.attrib.get("t")
        value_node = cell.find(f"{{{MAIN_NS}}}v")
        if cell_type == "inlineStr":
            return "".join(
                node.text or ""
                for node in cell.iter()
                if node.tag == f"{{{MAIN_NS}}}t"
            )
        if value_node is None or value_node.text is None:
            return None
        raw = value_node.text
        if cell_type == "s":
            return self.shared_strings[int(raw)]
        if cell_type in {"str", "e"}:
            return raw
        if cell_type == "b":
            return raw == "1"
        try:
            numeric = float(raw)
            return int(numeric) if numeric.is_integer() else numeric
        except ValueError:
            return raw

    def read_sheet(
        self,
        sheet_name: str,
        header: int = 0,
        usecols_end: int | None = None,
    ) -> pd.DataFrame:
        if sheet_name not in self.sheet_paths:
            raise KeyError(f"Unknown sheet {sheet_name!r}")
        rows: list[list[object]] = []
        header_values: list[object] | None = None
        row_position = -1
        with self.archive.open(self.sheet_paths[sheet_name]) as handle:
            for _, element in ET.iterparse(handle, events=("end",)):
                if element.tag != f"{{{MAIN_NS}}}row":
                    continue
                row_position += 1
                values: dict[int, object] = {}
                for cell in element.findall(f"{{{MAIN_NS}}}c"):
                    index = _column_index(cell.attrib["r"])
                    if usecols_end is None or index < usecols_end:
                        values[index] = self._cell_value(cell)
                width = usecols_end or (max(values, default=-1) + 1)
                dense = [values.get(index) for index in range(width)]
                if row_position == header:
                    header_values = dense
                elif row_position > header:
                    if header_values is None:
                        raise ValueError("Header row was not found")
                    width = len(header_values)
                    rows.append((dense + [None] * width)[:width])
                element.clear()
        if header_values is None:
            raise ValueError(f"Sheet {sheet_name!r} has no header row {header}")
        return pd.DataFrame(rows, columns=header_values)

