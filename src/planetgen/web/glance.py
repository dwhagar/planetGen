# planetgen/web/glance.py

"""
The "At a glance" rows of the pages that are not an item page (UX.188 to UX.194): a system's, and later a sector's and the
galaxy's. An item page reads its rows from the property sheet (`describe/sheet.py` `GLANCE`); a system or a sector has no
sheet, so its rows are worked out here from what its page already holds, in the order a captain asks: what is there,
where is the best place to stop, what is the hazard.

A row that would say none or 0 is left out unless its absence is the answer (a system with no planets says so).
"""

from planetgen.describe import model
from planetgen.describe.html import inline_html
from planetgen.web import item_pages
from planetgen.web.helpers import trusted_html, page_url
from planetgen.web.lib import habitability_symbols
from planetgen.web.lib.fmt import esc

GOOD_TIER = 1
"""int: The worst overall PHI-6 tier counted as a habitable world: Blue (0) or Green (1)."""


def _count(word, plural, number):
    return f"{number:,} {word if number == 1 else plural}"


def _worlds(detail):
    """`[(kind, body)]` for every planet and moon of a system's detail."""
    worlds = []
    for planet in detail.get("planets") or []:
        worlds.append(("planet", planet))
        worlds.extend(("moon", moon) for moon in planet.get("moons") or [])
    return worlds


def system_rows(detail, facilities, symbols):
    """
    The rows of a system page's panel: `[{"label", "html"}]`.

    Args:
        detail (dict): `GET /api/systems/<id>`.
        facilities (dict or list): `GET /api/systems/<id>/facilities` (its `items`, or the list itself).
        symbols (callable): Draws the symbols of a scored body (`web.habitability_symbol_row`).
    """
    rows = []
    stars = detail.get("stars") or []
    if stars:
        primary = next((s for s in stars if s.get("role") in ("single", "primary")), stars[0])
        rows.append({"label": "Stars", "html": esc(f"{_count('star', 'stars', len(stars))}, the main one {primary.get('star_type') or 'unclassified'}")})
    planets = detail.get("planets") or []
    rows.append({"label": "Planets", "html": esc(_count("planet", "planets", len(planets)) if planets else "none")})
    moons = sum(len(planet.get("moons") or []) for planet in planets)
    for label, number in (("Moons", moons), ("Asteroid belts", len(detail.get("belts") or [])),
                          ("Comets", len(detail.get("comets") or []))):
        if number:
            rows.append({"label": label, "html": esc(f"{number:,}")})
    scored = [(habitability_symbols.overall_tier(body), kind, body) for kind, body in _worlds(detail)
              if habitability_symbols.overall_tier(body) is not None]
    good = [item for item in scored if item[0] <= GOOD_TIER]
    rows.append({"label": "Habitable worlds", "html": esc(f"{len(good):,}" if good else "none")})
    if scored:
        tier, kind, body = min(scored, key=lambda item: item[0])
        link = f'<a href="{esc(item_pages.item_url(f"{kind}:{body["id"]}"))}">{esc(body["name"])}</a>'
        rows.append({"label": "Best world", "html": f"{link} {symbols(body)}"})
    held = facilities.get("items") if isinstance(facilities, dict) else facilities
    if held:
        rows.append({"label": "Facilities", "html": esc(f"{len(held):,}")})
    perimeters = [s["system_perimeter_km"] for s in stars if s.get("system_perimeter_km")]
    if perimeters:
        rows.append({"label": "System perimeter", "html": inline_html([model.qty("length", max(perimeters), "km")])})
    return [{"label": row["label"], "html": trusted_html(row["html"])} for row in rows]
