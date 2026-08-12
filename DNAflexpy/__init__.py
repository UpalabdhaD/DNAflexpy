"""DNAflexpy — DNA flexibility profiling from sequence."""

from DNAflexpy import lookup
from DNAflexpy.core import FlexProfiler, ProfileSet
from DNAflexpy.lookup import FeatureTable, default_table
from DNAflexpy.profile import FlexProfile

__version__ = "0.3.0.dev0"
__all__ = [
    "FlexProfiler",
    "FlexProfile",
    "ProfileSet",
    "FeatureTable",
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
