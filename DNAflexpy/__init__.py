"""DNAflexpy — DNA flexibility profiling from sequence."""

from DNAflexpy import lookup
from DNAflexpy.core import FlexProfiler, ProfileSet
from DNAflexpy.encode import FeatureMatrix
from DNAflexpy.lookup import FeatureTable, default_table
from DNAflexpy.results import FlexProfile

# The `encode` FUNCTION is deliberately not re-exported here. `DNAflexpy.encode`
# stays the submodule, so `from DNAflexpy.encode import encode` always works.
# Re-exporting the function would shadow the module it lives in -- the same
# trap `results.py` is named to avoid for `profile`. The normal way to reach
# it is the method: `prof.encode([...])`.

__version__ = "0.3.0.dev0"
__all__ = [
    "FlexProfiler",
    "FlexProfile",
    "ProfileSet",
    "FeatureTable",
    "FeatureMatrix",
    "default_table",
    "profile",
    "lookup",
    "__version__",
]


def profile(seq: str, feature: str = "DNaseI", window_size: int = 10):
    """Profile one sequence without constructing anything.

    The packaged table is memoised by `default_table`, so the first call
    pays the YAML parse and later calls do not.
    """
    return FlexProfiler(feature, window_size=window_size).profile(seq)
