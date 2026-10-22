# tests/test_info_box.py

"""
MAP.139, MAP.254: every map's selection panel is the one shared component (static/mappick.js), and each names the
info box over its map that shows the same spec; no map builds a panel or a box of its own.
"""

import re
from pathlib import Path

MAPS = Path(__file__).resolve().parents[1] / "planetgen" / "web" / "maps"
STATIC = Path(__file__).resolve().parents[1] / "html" / "static"
PANEL = re.compile(r'<aside class="starmap-info"([^>]*)>')


def _panels():
    for path in sorted(MAPS.glob("*.py")):
        text = path.read_text(encoding="utf-8")
        for match in PANEL.finditer(text):
            yield path, text, match.group(1)


def test_every_map_panel_names_its_info_boxes():
    found = list(_panels())
    assert {path.name for path, _, _ in found} >= {"galaxymap3d.py", "systemmap.py"}
    for path, text, attributes in found:
        named = re.search(r'data-infobox="([^"]+)"', attributes)
        assert named, f"{path.name}: a starmap-info panel with no data-infobox"
        for box_id in named.group(1).split():
            assert f'class="map-infobox" id="{box_id}"' in text, f"{path.name}: no info box {box_id}"


def test_only_the_shared_component_builds_the_info_box():
    for path in sorted(STATIC.glob("*.js")):
        if path.name == "mappick.js":
            continue
        text = path.read_text(encoding="utf-8")
        assert "map-infobox" not in text, f"{path.name} builds an info box of its own: use mappick.js"
        assert "starmap-info" not in text, f"{path.name} builds an info panel of its own: use mappick.js"
