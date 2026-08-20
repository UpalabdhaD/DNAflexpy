import numpy as np

import DNAflexpy
from DNAflexpy.lookup import default_table


def test_one_liner_needs_no_object():
    out = DNAflexpy.profile("ATGCGTACGTAGCTAGCGTAGCTAGT", feature="DNaseI",
                            window_size=10)
    assert isinstance(out, np.ndarray)
    assert len(out) == 17


def test_defaults_match_the_documented_cli_defaults():
    assert len(DNAflexpy.profile("ATGCGTACGTAGCTAGCGTAGCTAGT")) == 17


def test_matches_the_class_api():
    seq = "ATGCGTACGTAGCTAGCGTAGCTAGT"
    expected = DNAflexpy.FlexProfiler("trx", window_size=5).profile(seq)
    assert np.array_equal(DNAflexpy.profile(seq, feature="trx", window_size=5), expected)


def test_repeated_calls_parse_the_yaml_once(monkeypatch):
    """The 8ms defect must not sneak back in through the convenience path."""
    default_table.cache_clear()
    calls = {"n": 0}
    original = DNAflexpy.lookup.FeatureTable.from_yaml

    def counting(*args, **kwargs):
        calls["n"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(DNAflexpy.lookup.FeatureTable, "from_yaml", counting)
    for _ in range(1000):
        DNAflexpy.profile("ATGCGTACGT")
    assert calls["n"] == 1
