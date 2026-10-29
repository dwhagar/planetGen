# tests/test_spin_scope.py
"""MAP.328: bodies turn in the views that draw a body as a sphere and in no galactic view."""

from pathlib import Path

STATIC = Path(__file__).resolve().parents[1] / "html" / "static"

GALACTIC = ("galaxymap3d.js", "galaxystageview.js", "galaxyglow.js", "galaxyglowlayer.js", "galaxynebulae.js", "galaxyblocks.js",
            "galaxymarkerlayer.js", "galaxymarkers.js", "galaxysector.js")


def test_no_galactic_view_imports_the_turning_function():
    for name in GALACTIC:
        text = (STATIC / name).read_text()
        assert "spin.js" not in text and ".spin(" not in text, f"{name} turns something: the galactic views turn nothing"


def test_the_system_views_use_the_one_turning_function_not_a_rate_of_their_own():
    # galaxysystem.js is the system opened in the sector of the Galaxy Map: it draws its stars and bodies as spheres.
    for name in ("systemmap.js", "systemview3d.js", "galaxysystem.js"):
        text = (STATIC / name).read_text()
        assert "spin.js" in text, name
    assert "marker.rotation" not in (STATIC / "systemmap.js").read_text(), "the old rate for every body is gone"
