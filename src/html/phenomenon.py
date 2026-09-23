#!/usr/bin/env python3
# html/phenomenon.py

"""
Phenomenon detail: one exotic phenomenon's full info -- this project's
first detail/info page for a standalone nebula/asteroid field/black
hole/neutron star (until now these had no page of their own at all, only
a hover tooltip on the Sector Map/Galaxy Map -- see `lib/starmap.py`/
`lib/galaxymap.py`'s own docstrings). Reached from `phenomena.py`'s
listing, or directly from a Sector Map/Galaxy Map marker.

`GET /api/phenomena/<type>/<id>` returns the phenomenon's own row as-is
(each type's genuinely different column set, see `queryDb.
phenomenon_detail`'s docstring) rather than `phenomena.py`'s normalized
list shape -- `_FIELD_SPECS` below picks and formats the columns worth
showing per type.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import get_phenomenon
from fmt import esc, post_link
from page import nav_params, run
from phenomenonmap import render_phenomenon_map_panel

try:
    from stellarObjects.utils import pc_to_ly
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching galaxy.py's own identical pattern.
    pc_to_ly = None


def _ly(value):
    return f"{value:,.2f} ly"


def _pc_to_ly_text(radius_pc):
    if radius_pc is None:
        return None
    radius_ly = pc_to_ly(radius_pc) if pc_to_ly else radius_pc * 3.2616
    return _ly(radius_ly)


def _bool_text(value):
    return "Yes" if value else "No"


def _title_case(value):
    return (value or "").replace("_", " ").replace("-", " ").capitalize()


def _progenitor_text(value):
    # "Type Ia" is already correctly cased as stored -- _title_case would
    # mangle it (no "_"/"-" to split on, so its own .capitalize() just
    # lowercases the "Ia" to "ia"). "core-collapse" has no such special
    # casing to preserve, so it goes through _title_case normally.
    return value if value == "Type Ia" else _title_case(value)


_TYPE_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field",
    "black_hole": "Black Hole", "neutron_star": "Neutron Star",
    "supernova_remnant": "Supernova Remnant",
}

# (column, label, formatter) per type -- formatter takes the raw column
# value and returns display text, or None to omit the row entirely (a
# NULL galactic-motion column on an anchored remnant, which never applies
# here since phenomenon_detail's own callers only ever link a standalone
# row -- see that function's docstring -- but kept defensive rather than
# assuming).
_FIELD_SPECS = {
    "nebula": [
        ("nebula_type", "Nebula Type", _title_case),
        ("radius_ly", "Radius", _ly),
        ("composition", "Composition", esc),
        ("formation_cause", "Formation", esc),
        ("galactic_orbital_speed_kms", "Galactic Orbital Speed", lambda v: f"{v:,.1f} km/s"),
        ("galactic_orbital_period_gy", "Galactic Orbital Period", lambda v: f"{v:,.2f} Gy"),
    ],
    "asteroid_field": [
        ("density", "Density", _title_case),
        ("radius_ly", "Radius", _ly),
        ("composition_summary", "Composition", esc),
        ("galactic_orbital_speed_kms", "Galactic Orbital Speed", lambda v: f"{v:,.1f} km/s"),
        ("galactic_orbital_period_gy", "Galactic Orbital Period", lambda v: f"{v:,.2f} Gy"),
    ],
    "black_hole": [
        ("mass_solar", "Mass", lambda v: f"{v:,.2f} solar masses"),
        ("event_horizon_radius_km", "Event Horizon Radius", lambda v: f"{v:,.1f} km"),
        ("spin", "Spin (dimensionless)", lambda v: f"{v:.3f}"),
        ("has_accretion_disk", "Accretion Disk", _bool_text),
        ("temperature_k", "Hawking Temperature", lambda v: f"{v:.2e} K"),
        ("luminosity_w", "Luminosity", lambda v: f"{v:.2e} W"),
        ("age_gy", "Age", lambda v: f"{v:,.2f} Gy"),
        ("galactic_orbital_speed_kms", "Galactic Orbital Speed", lambda v: f"{v:,.1f} km/s" if v is not None else None),
        ("galactic_orbital_period_gy", "Galactic Orbital Period", lambda v: f"{v:,.2f} Gy" if v is not None else None),
    ],
    "neutron_star": [
        ("mass_solar", "Mass", lambda v: f"{v:,.2f} solar masses"),
        ("radius_km", "Radius", lambda v: f"{v:,.2f} km"),
        ("spin_period_ms", "Spin Period", lambda v: f"{v:,.2f} ms"),
        ("magnetic_field_gauss", "Magnetic Field", lambda v: f"{v:.2e} G"),
        ("pulsar_type", "Pulsar Type", _title_case),
        ("surface_temperature_k", "Surface Temperature", lambda v: f"{v:,.0f} K"),
        ("luminosity_w", "Luminosity", lambda v: f"{v:.2e} W"),
        ("age_gy", "Age", lambda v: f"{v:,.2f} Gy"),
        ("galactic_orbital_speed_kms", "Galactic Orbital Speed", lambda v: f"{v:,.1f} km/s" if v is not None else None),
        ("galactic_orbital_period_gy", "Galactic Orbital Period", lambda v: f"{v:,.2f} Gy" if v is not None else None),
    ],
    "supernova_remnant": [
        ("morphology", "Morphology", _title_case),
        ("progenitor_type", "Progenitor Type", _progenitor_text),
        ("age_years", "Age", lambda v: f"{v:,.0f} years"),
        ("radius_ly", "Radius", _ly),
        ("compact_remnant_kind", "Compact Remnant Left Behind", _title_case),
        ("galactic_orbital_speed_kms", "Galactic Orbital Speed", lambda v: f"{v:,.1f} km/s"),
        ("galactic_orbital_period_gy", "Galactic Orbital Period", lambda v: f"{v:,.2f} Gy"),
    ],
}


def _fields_html(phenomenon_type, detail):
    rows = []
    for column, label, formatter in _FIELD_SPECS.get(phenomenon_type, []):
        raw = detail.get(column)
        if raw is None:
            continue
        text = formatter(raw)
        if text is None:
            continue
        rows.append(f"<tr><td>{esc(label)}</td><td>{text}</td></tr>")
    return "".join(rows)


def handler():
    params = nav_params()
    db_name = params.get("db", "")
    phenomenon_type = params.get("type", "")
    phenomenon_id = params.get("id", "")

    detail = get_phenomenon(db_name, phenomenon_type, phenomenon_id)

    type_label = _TYPE_LABELS.get(phenomenon_type, phenomenon_type)
    phenomena_link = post_link("phenomena.py", {"db": db_name}, "Phenomena")
    breadcrumb = f'<p class="breadcrumb">{phenomena_link} &rarr; {esc(detail["name"])}</p>'

    badge_bits = [type_label]
    distance_text = _pc_to_ly_text(detail.get("galactic_radius_pc"))
    if distance_text:
        badge_bits.append(f"{distance_text} from Galactic Center")
    if detail.get("sector_id") is not None:
        sector_link = post_link("sector.py", {"db": db_name, "id": detail["sector_id"]}, esc(detail["sector_name"]))
        badge_bits.append(f"Sector: {sector_link}")
    badges_html = "<p class=\"badges\">" + "".join(f'<span class="badge">{bit}</span>' for bit in badge_bits) + "</p>"

    # Offered unconditionally for every OTHER type, unlike system.py's own
    # sector-gated nav buttons -- nav.py itself renders a clear "not
    # available" message for a phenomenon never placed in the galaxy, so
    # this page doesn't have to duplicate that same placement check just
    # to decide whether to show the button at all. A supernova remnant is
    # the one exception: its own table has no galaxy-frame placement
    # columns at ALL (see queryDb._SUPERNOVA_REMNANT_TABLE's docstring),
    # not merely "not placed yet" -- nav_between rejects it outright
    # (ValueError, queryDb._load_nav_phenomenon_endpoint), so offering the
    # button here would just walk a visitor into that error instead of
    # nav.py's normal unplaced-phenomenon message.
    if phenomenon_type == "supernova_remnant":
        nav_html = ""
        nav_hint_html = (
            '<p class="hint">A supernova remnant has no known position in the galaxy '
            "yet, so it can't be used for navigation.</p>"
        )
    else:
        nav_html = (
            post_link(
                "nav.py", {"db": db_name, "from": detail["id"], "from_kind": "phenomenon", "from_type": phenomenon_type},
                "Navigate from here", css_class="btn",
            )
            + post_link(
                "nav.py", {"db": db_name, "to": detail["id"], "to_kind": "phenomenon", "to_type": phenomenon_type},
                "Navigate to here", css_class="btn",
            )
        )
        nav_hint_html = ""

    map_html = render_phenomenon_map_panel(phenomenon_type, detail["name"], detail.get("radius_ly") or 0)
    fields_html = _fields_html(phenomenon_type, detail)
    body = f"""
<div class="page-subhead">{breadcrumb}{badges_html}{nav_html}</div>
{nav_hint_html}
{map_html}
<section class="panel">
<h2>{esc(type_label)} Data</h2>
<div class="table-scroll"><table>
  <tbody>
    {fields_html}
  </tbody>
</table></div>
</section>
"""
    return detail["name"], body


run(handler)
