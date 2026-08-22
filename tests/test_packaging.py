"""The archive and the new package must coexist in one install."""
import pathlib
import tomllib

ROOT = pathlib.Path(__file__).resolve().parent.parent


def test_setup_py_is_gone():
    """setup.py beside a [project] table makes setuptools error."""
    assert not (ROOT / "setup.py").exists()


def test_pyproject_declares_packages_explicitly():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    packages = cfg["tool"]["setuptools"]["packages"]
    # Explicit list, never find_packages(): find_packages() would sweep the
    # archive into the distribution while shipping none of its YAML.
    assert "rxv" in packages
    assert "rxv.DNAflexpy" in packages
    assert "DNAflexpy" in packages


def test_pyproject_ships_both_lookup_tables():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    data = cfg["tool"]["setuptools"]["package-data"]
    assert data["rxv.DNAflexpy"] == ["data/*.yaml"]
    # The new package also ships example data; the archive deliberately does
    # not, so the two entries are not expected to match.
    assert "data/*.yaml" in data["DNAflexpy"]
    assert "data/examples/*" in data["DNAflexpy"]


def test_manifest_points_at_the_archive():
    assert "rxv/DNAflexpy/data/lookupNEW.yaml" in (ROOT / "MANIFEST.in").read_text()


def test_plot_profiles_invokes_the_archived_cli():
    """It shells out to a module path; unfixed it would hit the new CLI."""
    src = (ROOT / "scripts/plot_profiles.py").read_text()
    assert "rxv.DNAflexpy.cli" in src
    assert '"DNAflexpy.cli"' not in src
