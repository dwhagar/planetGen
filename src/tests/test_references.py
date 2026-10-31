"""
docs/REFERENCES.md and docs/REFERENCES.html are written from one cache by scripts/build_references.py and must match it;
every listed paper carries a URL.
"""

import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_references as builder  # noqa: E402

CACHE = json.loads((ROOT / "docs" / "reference" / "references-data.json").read_text(encoding="utf-8"))


def test_pages_match_the_cache():
    result = subprocess.run([sys.executable, str(ROOT / "scripts" / "build_references.py"), "render", "--check"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_every_entry_has_a_url():
    missing = [key for key, entry in CACHE["entries"].items() if builder.entry_urls(entry)[0] is None]
    assert missing == []


def test_the_opening_paragraph_thanks_the_researchers():
    text = (ROOT / "docs" / "REFERENCES.md").read_text(encoding="utf-8")
    assert "would not be possible without the dedicated work of the researchers" in text
