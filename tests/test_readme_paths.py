"""Check that the README points only at files that exist, and that every script is documented.

Run from the repository root:  python -m pytest tests/
"""
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
README = (ROOT / "README.md").read_text(encoding="utf-8")
REPRODUCE = (ROOT / "reproduce.sh").read_text(encoding="utf-8")
DOCS = README + "".join((ROOT / d / "README.md").read_text(encoding="utf-8") for d in ("scripts", "results"))

SCRIPT_SUFFIXES = (".R", ".py", ".cpp")


def test_full_paths_exist():
    paths = set(re.findall(r"(?:scripts|models|data)/[\w./-]+\.(?:R|py|cpp)", README))
    assert paths, "no paths found"
    missing = sorted(p for p in paths if not (ROOT / p).is_file())
    assert not missing, f"README paths that do not exist: {missing}"


def test_bare_script_names_exist():
    names = set(re.findall(r"`(?:memo_0919/)?([\w]+\.(?:R|py))`?", README))
    shipped = {p.name for p in (ROOT / "scripts").rglob("*") if p.suffix in SCRIPT_SUFFIXES}
    shipped |= {p.name for p in (ROOT / "models").rglob("*") if p.suffix in SCRIPT_SUFFIXES}
    shipped.add("fetch_data.py")
    missing = sorted(n for n in names if n not in shipped)
    assert not missing, f"script names in README that are not shipped: {missing}"


def test_cited_run_folders_exist():
    # paragraphs that record removed folders cite names that are meant to be absent
    cited = "\n\n".join(p for p in DOCS.split("\n\n") if "were removed" not in p)
    runs = set(re.findall(r"([a-zA-Z][\w]*?_20\d{6}_\d{6})", cited))
    assert runs, "no run folders found"
    missing = sorted(r for r in runs if not (ROOT / "results" / "tDDM" / r).is_dir())
    assert not missing, f"run folders cited in README that are not shipped: {missing}"


def test_every_script_is_documented():
    scripts = [p for p in (ROOT / "scripts").rglob("*") if p.suffix in SCRIPT_SUFFIXES]
    undocumented = sorted(str(p.relative_to(ROOT)) for p in scripts if p.name not in DOCS)
    assert not undocumented, f"scripts that no README mentions: {undocumented}"


def test_every_shipped_run_folder_is_documented():
    folders = [p.name for p in (ROOT / "results" / "tDDM").iterdir()
               if p.is_dir() and re.search(r"_20\d{6}_\d{6}$", p.name)]
    # folders cited by a wildcard family in the README (e.g. diag_quadrature_refit_20260927_*)
    families = set(re.findall(r"([a-zA-Z][\w]*?_20\d{6})_\*", DOCS)) | set(re.findall(r"`([a-z_]+)_\*`", DOCS))
    def documented(name):
        return name in DOCS or any(name.startswith(f) for f in families)
    tracked = [f for f in folders if documented(f)]
    undocumented = sorted(set(folders) - set(tracked))
    assert not undocumented, f"run folders that no README mentions: {undocumented}"


def test_readme_images_exist():
    images = set(re.findall(r"!\[[^\]]*\]\(([^)]+)\)", README)) | set(re.findall(r'<img src="([^"]+)"', README))
    assert images, "no images found"
    missing = sorted(i for i in images if not (ROOT / i).is_file())
    assert not missing, f"README images that do not exist: {missing}"


def test_reproduce_script_paths_exist():
    paths = set(re.findall(r"(?:scripts|data)/[\w./-]+\.(?:R|py)", REPRODUCE))
    assert len(paths) > 30, f"too few script paths found in reproduce.sh: {len(paths)}"
    missing = sorted(p for p in paths if not (ROOT / p).is_file())
    assert not missing, f"reproduce.sh paths that do not exist: {missing}"


@pytest.mark.parametrize("stage, n_steps", [("regenerate", 20), ("fits", 6), ("checks", 4), ("validation", 14)])
def test_reproduce_dry_run_lists_every_step(stage, n_steps):
    if stage == "validation" and shutil.which("Rscript") is None:
        pytest.skip("validation reads the fitted onsets with R")
    out = subprocess.run(["bash", str(ROOT / "reproduce.sh"), stage, "--dry-run", "--cores", "3"],
                         capture_output=True, text=True, check=True).stdout
    commands = [line for line in out.splitlines() if line.startswith(("Rscript ", "python "))]
    assert len(commands) == n_steps, out
    assert all((ROOT / c.split()[1]).is_file() for c in commands), out
    if stage in ("fits", "validation"):
        assert all("--cores 3" in c for c in commands if "--cores" in c)
    if stage == "validation":
        assert "--tau-common 0.891" in out and "--tau-base 0.804" in out
