import numpy as np
import pandas as pd

from compositionality.common import bh_fdr, close_rows, clr_transform


def test_close_rows_sums_to_one_and_preserves_zeros():
    source = pd.DataFrame([[2.0, 0.0, 6.0], [1.0, 1.0, 2.0]])
    closed = close_rows(source)
    assert np.allclose(closed.sum(axis=1), 1.0)
    assert closed.iloc[0, 1] == 0.0
    assert np.allclose(closed.iloc[0].to_numpy(), [0.25, 0.0, 0.75])


def test_clr_uses_training_replacements_for_held_out_zeros():
    train = pd.DataFrame(
        {"a": [0.10, 0.20, 0.0], "b": [0.90, 0.80, 1.0]},
        index=["t1", "t2", "t3"],
    )
    test = pd.DataFrame(
        {"a": [0.0], "b": [1.0]},
        index=["held_out"],
    )
    _, transformed, replacements = clr_transform(train, test)
    assert replacements["a"] == 0.10
    expected = np.log([0.10, 1.0])
    expected = expected - expected.mean()
    assert np.allclose(transformed.loc["held_out"].to_numpy(), expected)


def test_bh_fdr_known_example_and_nan_handling():
    observed = bh_fdr([0.01, 0.04, 0.03, np.nan])
    assert np.allclose(observed[:3], [0.03, 0.04, 0.04])
    assert np.isnan(observed[3])
