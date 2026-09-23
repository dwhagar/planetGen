#!/usr/bin/env python3
# html/galaxy3d.py

"""
Galaxy Map (3D): the interactive, free-flying counterpart to the flat
overview map (`galaxy.py`) -- a real perspective camera
(`lib/galaxymap3d.py`/`static/galaxymap3d.js`) instead of a fixed, face-on
SVG projection, so a visitor can rotate, dolly, and click through the
galaxy's actual 3D shape rather than only ever seeing it from directly
"above" the disk.

This page itself only makes the two `apiclient` calls its first paint
needs (`get_galaxy_shape`, then one `get_galaxy_view` call for the
zoomed-all-the-way-out starting view) and hands both to
`lib.galaxymap3d.render_galaxy_map3d_panel` -- every subsequent view (as
the visitor's camera moves) is fetched directly by the page's own
client-side JS from `galaxy_view.py`, never through this handler again.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import get_galaxy_shape, get_galaxy_view
from fmt import esc, post_link
from galaxymap3d import render_galaxy_map3d_panel, view_radius_bounds
from page import nav_params, run

try:
    from stellarObjects.program_constants import DEFAULT_SECTOR_EDGE_LY
    from stellarObjects.utils import ly_to_pc
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching every other page's identical pattern
    # (see e.g. galaxy.py's own top-of-file try/except).
    DEFAULT_SECTOR_EDGE_LY = 11.5
    def ly_to_pc(ly):
        return ly / 3.2616


def handler():
    params = nav_params()
    db_name = params.get("db", "")

    galaxy_shape = get_galaxy_shape(db_name)
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    _min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)

    initial_view = get_galaxy_view(db_name, 0.0, 0.0, 0.0, max_radius)

    map_html = render_galaxy_map3d_panel(db_name, galaxy_shape, edge_pc, initial_view)

    title = "Galaxy Map (3D)"
    breadcrumb_html = (
        f'<p class="breadcrumb">{post_link("browse.py", {"db": db_name}, esc(db_name))} &rarr; '
        f'{post_link("galaxy.py", {"db": db_name}, "Galaxy Map")} &rarr; {esc(title)}</p>'
    )
    body = f"""
<div class="page-subhead">{breadcrumb_html}</div>
{map_html}
<script type="module" src="static/galaxymap3d.js"></script>
"""
    return title, body


run(handler)
