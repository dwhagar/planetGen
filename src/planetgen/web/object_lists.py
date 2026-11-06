# planetgen/web/object_lists.py

"""
The lists of every kind of object (UX.141): `/lists` names them, `/lists/<kind>` is one as a data table (UX.41, keyset
paging as in PERF.75) that reads the object index through `GET /api/object-lists/<kind>`, so the thousandth page costs
what the first does. Systems have their own richer list at `/systems`, phenomena theirs at `/phenomena`; they are
named here beside the rest.
"""

import copy

from flask import abort, redirect, request

from planetgen.db import objectlist, search_values
from planetgen.util.format import format_number
from planetgen.web.lib import apiclient
from planetgen.web.lib.datatable import Column, Result, Table, read_from_start
from planetgen.web.lib.fmt import distance_km_html, measured_html
from planetgen.web.lib.habitability_cells import COLUMN as HABITABILITY_COLUMN, habitability_cell

from . import bp, tables
from .helpers import crumb, db_name, page_url, render_page
from .searchpage import SearchState, tidy_value_args, value_chips, value_form

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
    "hypervelocity_star": ("hypervelocity star", "hypervelocity stars", "Hypervelocity stars"),
}
"""dict: Each kind with a list page, to its singular, plural, and heading (facilities have no page of their own)."""

COLUMNS = [Column("name", "Name"), Column("class", "Class"), Column("mass", "Mass"), Column("radius", "Radius"),
           Column("sector", "Sector", sortable=False)]

HABITABILITY_KINDS = ("planet", "moon")
"""tuple: The kinds whose list has a Habitability column of the PHI-6 symbols (UX.100), before the Sector column."""

_SPIN_COLUMNS = (Column("rotation_period", "Rotation period", sortable=False), Column("axial_tilt", "Axial tilt", sortable=False),
                 Column("internal_heat_flux", "Internal heat flow", sortable=False))

PHYSICS_COLUMNS = {
    "star": (Column("temperature", "Temperature"), Column("luminosity", "Luminosity"), Column("age", "Age"),
             Column("rotation_period", "Rotation period", sortable=False)),
    "planet": (Column("orbit_distance", "Orbit"), Column("eccentricity", "Eccentricity", sortable=False),
               Column("surface_temperature", "Surface temperature"), *_SPIN_COLUMNS),
    "moon": (Column("orbit_distance", "Orbit"), Column("eccentricity", "Eccentricity", sortable=False),
             Column("surface_temperature", "Surface temperature"), *_SPIN_COLUMNS),
}
"""dict: The stored physics values a kind's list shows after the radius (UX.187), in the order of `db.objectlist.PHYSICS_COLUMNS`.
A column named for a key of `db.search_values` sorts by that value (UX.228); the eccentricity has no key, so it does not sort."""

VALUES_COLUMN = Column("values", "Values", sortable=False)
"""Column: The cell that carries the physical values a range or an order asked for and no column of the list shows (UX.228); it is on
the list only while one is in use."""

_SORT_KEYS = {"name": "name", "class": "class", "mass": "mass", "radius": "radius"}


def _title(value):
    return (value or "").replace("_", " ").replace("-", " ").capitalize()


def _columns(kind):
    physics = list(PHYSICS_COLUMNS.get(kind, ()))
    if kind not in HABITABILITY_KINDS:
        return [*COLUMNS[:-1], *physics, COLUMNS[-1]]
    return [*COLUMNS[:-1], *physics, Column("habitability", HABITABILITY_COLUMN, sortable=False), COLUMNS[-1]]


def _plain(value, formatter):
    """A physics cell: the value in the reader's units, or a muted dash for one the row does not hold."""
    return {"text": formatter(value)} if value is not None else {"text": "\u2013", "muted": True}


def _eccentricity(value):
    return f"{value:.3f}"


def _physics_cells(item, kind):
    """The cells of `PHYSICS_COLUMNS[kind]` for one item (UX.187)."""
    if kind == "star":
        return [_plain(item.get("temperature_k"), lambda v: measured_html("K", v, companions=False)),
                _plain(item.get("luminosity_w"), lambda v: measured_html("W", v)),
                _plain(item.get("age_gy"), lambda v: measured_html("Gyr", v)),
                _plain(item.get("rotation_period_hours"), lambda v: measured_html("h", v))]
    if kind in ("planet", "moon"):
        return [_plain(item.get("distance_km"), distance_km_html),
                {"text": _eccentricity(item["eccentricity"])} if item.get("eccentricity") is not None else {"text": "\u2013", "muted": True},
                _plain(item.get("surface_temperature_k"), lambda v: measured_html("K", v, companions=False)),
                _plain(item.get("rotation_period_hours"), lambda v: measured_html("h", v)),
                _plain(item.get("axial_tilt_deg"), lambda v: f"{format_number(v, '.4g')}\u00b0"),
                _plain(item.get("internal_heat_flux_w_m2"), lambda v: measured_html("W_m2", v))]
    return []


def _row(item, kind):
    if item.get("uncharted"):
        # MAP.400: no sector has stored it yet. A scattered row has the page of its scatter row; a field cloud has none.
        name = {"text": item["name"]}
        if item.get("scatter_id") is not None:
            name["href"] = page_url("uncharted_phenomenon", scatter_id=item["scatter_id"])
    else:
        name = {"href": page_url("object_page", ref=f"{item['kind']}:{item['id']}"), "text": item["name"]}
    sector = ({"text": f"Sector {item['sector_id']}", "href": page_url("sector", sector_id=item["sector_id"])}
              if item["sector_id"] is not None else {"text": "None", "muted": True})
    cells = [name,
             {"text": _title(item["object_class"])} if item["object_class"] else {"text": "None", "muted": True},
             {"text": measured_html("kg", item["mass_kg"])},
             {"text": distance_km_html(item["radius_km"])}]
    cells += _physics_cells(item, kind)
    if kind in HABITABILITY_KINDS:
        cells.append(habitability_cell(item, name["href"]))
    return cells + [sector]


def _value_state():
    """The request's physical-value ranges and order, read as the search page reads them (UX.227): a bound that is not a number, an inverted range
    and a unit that does not measure the value are ignored."""
    return SearchState.from_args(request.args)


def _ranges(kind):
    """`{key: (least, most, unit or None)}` of the physical-value ranges on this kind's list, or `None` for a kind that takes none."""
    if kind not in objectlist.VALUE_KINDS:
        return None
    return _value_state().value_ranges().get(kind) or None


def _order(kind, state):
    """`(sort, descending)` for the API: a column's own order, or the value order of `<kind>_value_sort` while the Name column is at its default (UX.228)."""
    ordered = _value_state().value_sorts().get(kind)
    if ordered and state.sort == "name" and not state.descending:
        return ordered.lstrip("-"), ordered.startswith("-")
    own = set(objectlist.value_sort_keys(kind))
    return _SORT_KEYS.get(state.sort) or (state.sort if state.sort in own else None), state.descending


def _extra_keys(kind, state):
    """The `(key, unit)` pairs the Values column carries: the values a range or an order is on that no column of the list shows already."""
    ranges = _ranges(kind)
    sort = _order(kind, state)[0]
    if kind not in objectlist.VALUE_KINDS or not (ranges or sort in objectlist.value_sort_keys(kind)):
        return []
    resolved = {key: (least, most, unit or search_values.KEYS[key].default_unit(kind)) for key, (least, most, unit) in (ranges or {}).items()}
    shown = {"mass", *(column.key for column in PHYSICS_COLUMNS.get(kind, ()))}
    return [(key, unit) for key, unit in search_values.used(kind, resolved, sort) if key not in shown]


def _values_cell(item, keys):
    """The Values cell: each value the list was asked about, with its label and unit, or a muted dash for a body without any."""
    parts = []
    for key, _unit in keys:
        found = (item.get("values") or {}).get(key)
        if found and found["value"] is not None:
            tail = f" {found['symbol']}" if found["symbol"] else ""
            parts.append(f"{search_values.KEYS[key].label.split(' (')[0]} {found['value']:.4g}{tail}")
    return {"text": "; ".join(parts)} if parts else {"text": "\u2013", "muted": True}


def _loader(kind):
    def load(state, limit, offset, want_facets, after=None):
        sort, descending = _order(kind, state)
        ranges = _ranges(kind)
        extra = _extra_keys(kind, state)
        envelope, note = read_from_start(
            lambda start: apiclient.get_object_list(
                db_name(), kind, limit=limit, sort=sort, descending=descending,
                after=after, start=start or None, ranges=ranges),
            state.start)
        items = envelope["items"]
        rows = [_row(item, kind) for item in items]
        if extra:
            for row, item in zip(rows, items):
                row.insert(len(row) - 1, _values_cell(item, extra))
        return Result(rows, envelope["total"], None, envelope.get("total_capped", False), envelope["next"], note)
    return load


LIST_TABLES = {
    kind: tables.register(Table(
        f"list-{kind}", labels[2], _columns(kind), _loader(kind), prefix=f"{kind}_", noun=labels[:2], keyset=True,
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


VALUE_PARAMS = {kind: tuple(name for key in search_values.KEYS if kind in search_values.KEYS[key].entities
                            for name in (f"{kind}_min_{key}", f"{kind}_max_{key}", f"{kind}_{key}_unit"))
                for kind in objectlist.VALUE_KINDS}
"""dict: Kind -> every range parameter of its list (UX.228), kept in the address while the table is paged and sorted."""


def _kept_params(kind):
    """The physical-value parameters present in the address, to carry through the table's links and its JSON address (UX.228)."""
    if kind not in objectlist.VALUE_KINDS:
        return ()
    wanted = {*VALUE_PARAMS[kind], f"{kind}_value_sort", f"{kind}_other_gas", f"{kind}_other_gas_min", f"{kind}_other_gas_max",
              f"{kind}_other_gas_unit"}
    return tuple(name for name in request.args if name in wanted)


def _chips(kind, ranges):
    """The "remove this filter" chips of the ranges and the order in use: `[{label, url}]` (UX.228)."""
    pairs = list(request.args.items(multi=True))
    chips = value_chips({kind: ranges} if ranges else {}, pairs, page_url("object_list", kind=kind), with_entity=False)
    ordered = _value_state().value_sorts().get(kind)
    if ordered:
        label = search_values.KEYS[ordered.lstrip("-")].label.split(" (")[0].lower()
        remaining = [(key, value) for key, value in pairs if key != f"{kind}_value_sort"]
        chips.append({"label": f"ordered by {label}, {'high to low' if ordered.startswith('-') else 'low to high'}",
                      "url": page_url("object_list", kind=kind, **dict(remaining))})
    return chips


def _canonical_url(kind):
    """The short form of this request's address (a submitted form sends every box, filled or not, and every row's unit), or `None` when it already is
    that (`searchpage.tidy_value_args`, UX.228)."""
    if kind not in objectlist.VALUE_KINDS:
        return None
    pairs, changed = tidy_value_args(request.args, (kind,))
    return page_url("object_list", kind=kind, **dict(pairs)) if changed else None


@bp.route("/lists/<kind>")
def object_list(kind):
    """One kind of object as a data table: sort by name, class, mass or radius (a star, planet or moon list also by its physical values and filtered by ranges on them, UX.228) and scroll through the pages."""
    table = LIST_TABLES.get(kind)
    if table is None:
        abort(404)
    canonical = _canonical_url(kind)
    if canonical:
        return redirect(canonical)
    heading = KIND_LABELS[kind][2]
    state = tables.parse_state(table, request.args)
    if _extra_keys(kind, state):
        table = copy.copy(table)
        table.columns = [*table.columns[:-1], VALUES_COLUMN, table.columns[-1]]
    rendered = tables.render(table, request.path, anchor="object-list", keep=_kept_params(kind), **{
        name: request.args[name] for name in _kept_params(kind)})
    groups = [group for group in value_form(_value_state()) if group["entity"] == kind]
    return render_page("object_list.html", title=heading, section="lists",
                       breadcrumbs=[crumb("Lists", "lists"), crumb(heading)],
                       description=f"Every {KIND_LABELS[kind][0]} in this charted galaxy.", table=rendered,
                       total=rendered["total"], noun=KIND_LABELS[kind],
                       value_group=groups[0] if groups else None, value_chips=_chips(kind, _ranges(kind)),
                       clear_url=page_url("object_list", kind=kind),
                       value_hidden=[(name, request.args[name]) for name in (f"{table.prefix}sort", f"{table.prefix}order")
                                     if request.args.get(name)],
                       filter_url=page_url("search", type=kind, **{f"{kind}s_sort": "phi4", f"{kind}s_order": "desc"}) if kind in HABITABILITY_KINDS else None)
