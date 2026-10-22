# planetgen/web/object_lists.py

"""
The lists of every kind of object (UX.141): `/lists` names them, `/lists/<kind>` is one as a data table (UX.41, keyset
paging as in PERF.75) that reads the object index through `GET /api/object-lists/<kind>`, so the thousandth page costs
what the first does. Systems have their own richer list at `/systems`, phenomena theirs at `/phenomena`; they are
named here beside the rest.
"""

from flask import abort, request

from planetgen.web.lib import apiclient
from planetgen.web.lib.datatable import Column, Result, Table, read_from_start
from planetgen.web.lib.fmt import distance_km_html, measured_html

from . import bp, tables
from .helpers import crumb, db_name, page_url, render_page

KIND_LABELS = {
    "star": ("star", "stars", "Stars"),
    "planet": ("planet", "planets", "Planets"),
    "moon": ("moon", "moons", "Moons"),
    "belt": ("asteroid belt", "asteroid belts", "Asteroid belts"),
    "comet": ("comet", "comets", "Comets"),
    "nebula": ("nebula", "nebulae", "Nebulae"),
    "asteroid_field": ("asteroid field", "asteroid fields", "Asteroid fields"),
    "black_hole": ("black hole", "black holes", "Black holes"),
    "neutron_star": ("neutron star", "neutron stars", "Neutron stars"),
    "supernova_remnant": ("supernova remnant", "supernova remnants", "Supernova remnants"),
    "rogue_planet": ("rogue planet", "rogue planets", "Rogue planets"),
    "interstellar_comet": ("interstellar comet", "interstellar comets", "Interstellar comets"),
    "quasar": ("quasar", "quasars", "Quasars"),
}
"""dict: Each kind with a list page, to its singular, plural and heading (facilities have no page of their own)."""

COLUMNS = [Column("name", "Name"), Column("class", "Class"), Column("mass", "Mass"), Column("radius", "Radius"),
           Column("sector", "Sector", sortable=False)]

_SORT_KEYS = {"name": "name", "class": "class", "mass": "mass", "radius": "radius"}


def _title(value):
    return (value or "").replace("_", " ").replace("-", " ").capitalize()


def _row(item):
    name = {"text": item["name"], "href": page_url("object_page", ref=f"{item['kind']}:{item['id']}")}
    sector = ({"text": f"Sector {item['sector_id']}", "href": page_url("sector", sector_id=item["sector_id"])}
              if item["sector_id"] is not None else {"text": "None", "muted": True})
    return [name,
            {"text": _title(item["object_class"])} if item["object_class"] else {"text": "None", "muted": True},
            {"text": measured_html("kg", item["mass_kg"])},
            {"text": distance_km_html(item["radius_km"])},
            sector]


def _loader(kind):
    def load(state, limit, offset, want_facets, after=None):
        envelope, note = read_from_start(
            lambda start: apiclient.get_object_list(
                db_name(), kind, limit=limit, sort=_SORT_KEYS.get(state.sort), descending=state.descending,
                after=after, start=start or None),
            state.start)
        items = envelope["items"]
        return Result([_row(item) for item in items], envelope["total"], None, envelope.get("total_capped", False),
                      envelope["next"], note)
    return load


LIST_TABLES = {
    kind: tables.register(Table(
        f"list-{kind}", labels[2], COLUMNS, _loader(kind), prefix=f"{kind}_", noun=labels[:2], keyset=True,
        default_sort="name", jump={"name": "A letter or the start of a name"}))
    for kind, labels in KIND_LABELS.items()
}
"""dict: Each kind's table, registered under `list-<kind>`."""


@bp.route("/lists")
def lists():
    """Every list on the site, one line each with how many objects it holds."""
    entries = [{"label": "Systems", "url": page_url("systems"), "total": None},
               {"label": "Sectors", "url": page_url("sectors"), "total": None},
               {"label": "Phenomena", "url": page_url("phenomena"), "total": None}]
    for kind, labels in KIND_LABELS.items():
        envelope = apiclient.get_object_list(db_name(), kind, limit=1)
        entries.append({"label": labels[2], "url": page_url("object_list", kind=kind), "total": envelope["total"],
                        "capped": envelope.get("total_capped", False)})
    return render_page("lists.html", title="Lists", section="lists", breadcrumbs=[crumb("Lists")],
                       description="Every list of objects in this charted galaxy.", entries=entries)


@bp.route("/lists/<kind>")
def object_list(kind):
    """One kind of object as a data table: sort by name, class, mass or radius and scroll through the pages."""
    table = LIST_TABLES.get(kind)
    if table is None:
        abort(404)
    heading = KIND_LABELS[kind][2]
    rendered = tables.render(table, request.path, anchor="object-list")
    return render_page("object_list.html", title=heading, section="lists",
                       breadcrumbs=[crumb("Lists", "lists"), crumb(heading)],
                       description=f"Every {KIND_LABELS[kind][0]} in this charted galaxy.", table=rendered,
                       total=rendered["total"], noun=KIND_LABELS[kind])
