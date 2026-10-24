# planetgen/web/lib/systempage.py

"""
The HTML pieces the system page (`/system/<id>`, `web/system_pages.py`) and the object pages (`web/item_pages.py`)
share about a system's facilities (schema v42: starbases, colonies, outposts): `facility_row` formats one for the
system page's Facilities panel and `facilities_html` lists the ones a body hosts on that body's own page.

The system's own stars, planets, moons, belts and comets are the summary tree (`web/lib/objecttree.py`, UX.173); the
rows that used to open onto each body's description here are gone, the description is on the body's own page.

Every value taken from the database is escaped with `fmt.esc`; the only raw markup comes from `fmt`'s unit formatters
(built from floats and fixed unit literals, e.g. `<sup>`). So the page wraps what these return in `trusted_html`.
"""

from planetgen.web.lib.fmt import distance_km_html, esc, period_years_html, speed_kms_html


_FACILITY_PLACEMENT_LABELS = {
    "terrestrial": "On the surface",
    "orbital": "In orbit",
    "asteroid": "Among the asteroids",
    "standalone": "Parked in space",
}
"""dict: How each `facilities.placement` reads on the page."""


def facility_kind_label(kind):
    """A `facilities.kind` as a label: "mining-colony" -> "Mining colony"."""
    return (kind or "").replace("-", " ").capitalize()


def facility_row(facility):
    """
    One facility's display values, every one already escaped or built from
    numbers (the page wraps them in `trusted_html`): `name`, `kind`,
    `placement` and, for an orbital facility, `distance`, `period` and
    `speed` (else `None`). The orbit is the stored one the API returns,
    worked out from the host's mass when it was saved.
    """
    orbital = facility.get("orbit_distance_km") is not None
    return {
        "name": esc(facility["name"]),
        "kind": esc(facility_kind_label(facility["kind"])),
        "placement": _FACILITY_PLACEMENT_LABELS.get(facility["placement"], esc(facility["placement"])),
        "distance": distance_km_html(facility["orbit_distance_km"]) if orbital else None,
        "period": period_years_html(facility["orbit_period_years"])
        if orbital and facility.get("orbit_period_years") is not None else None,
        "speed": speed_kms_html(facility["orbital_speed_kms"])
        if orbital and facility.get("orbital_speed_kms") is not None else None,
    }


def facilities_html(facilities, link=None):
    """The facilities a body hosts as a list for its own page, or `""` for none. `link` turns a facility into the URL of
    its own page (UX.183), and its name into a link to it."""
    if not facilities:
        return ""
    items = []
    for facility in facilities:
        row = facility_row(facility)
        where = row["placement"]
        if row["distance"]:
            where += f' at {row["distance"]}'
            where += "".join(f", {part}" for part in (row["period"], row["speed"]) if part)
        name = f'<a href="{esc(link(facility))}">{row["name"]}</a>' if link else row["name"]
        items.append(f'<li><span class="facility-name">{name}</span> '
                     f'<span class="stat">{row["kind"]}</span> {where}</li>')
    return f'<ul class="facility-list" aria-label="Facilities">{"".join(items)}</ul>'


def facilities_by_host(facilities):
    """`{(host_type, host_id): [facility, ...]}`."""
    by_host = {}
    for facility in facilities or ():
        by_host.setdefault((facility["host_type"], facility["host_id"]), []).append(facility)
    return by_host
