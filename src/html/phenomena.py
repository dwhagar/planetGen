#!/usr/bin/env python3
# html/phenomena.py

"""
Phenomena list: every exotic phenomenon (nebula/asteroid field/black hole/
neutron star/supernova remnant/rogue planet/interstellar comet --
`queryDb._PHENOMENON_TABLES` plus `_UNPLACED_PHENOMENON_TABLES`) in the
current database, one flat table across every sector and regardless of
galaxy placement -- `GET /api/phenomena`, the paginated counterpart to
`galaxy.py`'s own plotted-only view (`GET /api/galaxy/phenomena`, which
only shows the galaxy-placed subset as dots on the map -- a supernova
remnant/rogue planet/interstellar comet never appears there, or in its "On
Galaxy Map" column below, since none of those three tables have any
galaxy-frame placement columns at all; see `_SUPERNOVA_REMNANT_TABLE`'s
docstring).

Each row links to `phenomenon.py`, this project's first detail/info page
for a standalone phenomenon -- until now these had no page of their own at
all (`lib/starmap.py`/`lib/galaxymap.py`'s own tooltips said so
explicitly); both of those maps' phenomenon markers now link here too (see
their own modules).

Paged through the API (`GET /api/phenomena`'s `limit`/`offset`),
`pagination.PAGE_SIZE` rows at a time with the site's shared pager
(`lib/pagination.py`), `page` carrying the current page number.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from apiclient import get_phenomena
from fmt import esc, post_link
from page import nav_params, run
from pagination import fetch_page, parse_page, render_pagination

_TYPE_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field",
    "black_hole": "Black Hole", "neutron_star": "Neutron Star",
    "supernova_remnant": "Supernova Remnant",
    "rogue_planet": "Rogue Planet", "interstellar_comet": "Interstellar Comet",
}


def handler():
    params = nav_params()
    db_name = params.get("db", "")

    envelope, page = fetch_page(
        lambda limit, offset: get_phenomena(db_name, limit=limit, offset=offset),
        parse_page(params.get("page")),
    )
    items, total = envelope["items"], envelope["total"]
    pager_html = render_pagination("phenomena.py", {"db": db_name}, "page", page, total, label="Phenomena pages")

    def _row_html(row):
        radius_text = f"{row['radius_ly']:,.2f} ly" if row["radius_ly"] else "&ndash;"
        if row["sector_id"] is not None:
            sector_link = post_link("sector.py", {"db": db_name, "id": row["sector_id"]}, esc(row["sector_name"]))
            sector_cell = f"<td>{sector_link}</td>"
        else:
            sector_cell = "<td><em>None</em></td>"
        phenomenon_link = post_link(
            "phenomenon.py", {"db": db_name, "type": row["type"], "id": row["id"]}, esc(row["name"])
        )
        return (
            "<tr>"
            f'<td>{phenomenon_link}</td>'
            f'<td>{esc(_TYPE_LABELS.get(row["type"], row["type"]))}</td>'
            f'<td>{esc((row["descriptor"] or "").replace("_", " ").capitalize())}</td>'
            f'<td>{radius_text}</td>'
            f'{sector_cell}'
            f'<td>{"Yes" if row["placed"] else "No"}</td>'
            "</tr>"
        )

    rows = "".join(_row_html(row) for row in items) or '<tr><td colspan="6"><em>None</em></td></tr>'

    count_badge = f"{total} phenomen{'a' if total != 1 else 'on'} found"
    body = f"""
<p class="breadcrumb">{post_link("browse.py", {"db": db_name}, esc(db_name))} &rarr; Phenomena</p>
<p class="badges"><span class="badge">{count_badge}</span></p>
<section class="panel">
<h2>Phenomena</h2>
<div class="table-scroll"><table>
  <thead>
    <tr><th>Name</th><th>Type</th><th>Descriptor</th><th>Radius</th><th>Sector</th><th>On Galaxy Map</th></tr>
  </thead>
  <tbody>
    {rows}
  </tbody>
</table></div>
{pager_html}
</section>
"""
    return "Phenomena", body


run(handler)
