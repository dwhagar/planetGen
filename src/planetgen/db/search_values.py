# planetgen/db/search_values.py

"""
The numeric keys of the search: range filters and sorts on the stored physical values (UX.226)
==============================================================================================

`SEARCH_TAG_FACETS` are tags and the old sorts are the habitability scores; this registry adds the physical quantities a person
asks for in numbers: "gravity 0.8 to 1.2 g and a dose under 2 mSv a year". A key names one stored quantity of a star, a planet or a
moon (or a count or a gas read from a child table), the unit its bounds are given in, and the SQL that reads it. The engine is one:

* `ranges_from_args` reads the request: `<entity>_min_<key>` and `<entity>_max_<key>` (`planet_min_gravity=0.8`), in the key's unit
  or in another unit of the same kind named by `<entity>_<key>_unit` (`planet_surface_pressure_unit=kPa`). A bound may be left out.
* `clauses` turns the bounds into SQL conditions on the stored column, converted from the unit given to the unit the column is
  stored in (`db.column_units`, `physics.unit_registry`); a body with no value for the key (NULL) is not in the range.
* `order_by` sorts a panel by a key (`planet_sort=gravity` smallest first, `planet_sort=-gravity` largest first; a body with no value
  comes last either way).
* `select_columns` and `row_values` carry the value of each key in use on the rows, in the unit the request used.

The values are read from the stored columns as they are (no value is derived or changed here), and every bound and every number
returned has its unit: the key's own, or the one the request named.

Keys (`KEYS`): gravity, surface temperature, surface pressure, orbit distance, day length (the solar day), rotation period, year
length, surface dose, magnetic moment and ocean depth of a planet or a moon; the number of moons of a planet; the partial
pressure of each gas of the gas table (`gas_<formula>`, from the child tables `planet_gases` and `moon_gases`; no row is 0); the
mass, luminosity, age and temperature of a star (`mass` is the key of all three: a star in solar masses, a body in Earth masses). A radius keeps its own parameters (`<entity>_min_radius_km`).
"""

import dataclasses
import math

from planetgen.db import column_units
from planetgen.physics import gases, unit_registry

ENTITIES = ("star", "planet", "moon")
"""tuple: The result panels that take numeric keys."""

TABLE_OF = {"star": "stars", "planet": "planets", "moon": "moons"}
"""dict: The table of each entity."""

GAS_PREFIX = "gas_"
"""str: The front of a gas key (`gas_CO2`)."""


@dataclasses.dataclass(frozen=True)
class SearchValue:
    """
    One key of the search.

    Attributes:
        key (str): The name in the parameters (`gravity`).
        label (str): What the person reads.
        entities (tuple): The panels that have it.
        column (str or None): The column of the entity's table, or `None` for a quantity read by `expression`.
        units (dict): Entity -> the unit key the bounds are given in when the request names none.
        stored (dict): Entity -> the unit key the quantity is stored in, for a quantity with no column of its own.
        expression (callable or None): `(entity, alias) -> SQL` for a quantity read from somewhere else than a column.
    """

    key: str
    label: str
    entities: tuple
    column: str | None
    units: dict
    stored: dict = dataclasses.field(default_factory=dict)
    expression: object = None

    def sql(self, entity, alias):
        """The SQL that reads the quantity of `entity` (its table alias `alias`): the stored value, in the stored unit."""
        if self.expression is not None:
            return self.expression(entity, alias)
        return f"{alias}.{self.column}"

    def stored_unit(self, entity):
        """The unit key the quantity is stored in."""
        if entity in self.stored:
            return self.stored[entity]
        return column_units.unit_of(TABLE_OF[entity], self.column)

    def default_unit(self, entity):
        """The unit key the bounds are given in when the request names none."""
        return self.units[entity]


def _same(unit, entities):
    return {entity: unit for entity in entities}


def _moon_count(entity, alias):
    return f"(SELECT COUNT(*) FROM moons mc WHERE mc.planet_id = {alias}.id)"


def _gas(formula):
    quoted = formula.replace("'", "''")

    def expression(entity, alias):
        table, owner = ("planet_gases", "planet_id") if entity == "planet" else ("moon_gases", "moon_id")
        return f"COALESCE((SELECT g.partial_pa FROM {table} g WHERE g.{owner} = {alias}.id AND g.gas = '{quoted}'), 0)"

    return SearchValue(f"{GAS_PREFIX}{formula}", f"Partial pressure of {formula}", ("planet", "moon"), None,
                       _same("bar", ("planet", "moon")), _same("Pa", ("planet", "moon")), expression)


_BODIES = ("planet", "moon")


def _body(key, label, column, unit, entities=_BODIES):
    return SearchValue(key, label, entities, column, _same(unit, entities))


KEYS = {value.key: value for value in (
    _body("gravity", "Surface gravity", "gravity_g", "g0"),
    _body("surface_temperature", "Surface temperature", "surface_temperature_k", "K"),
    _body("surface_pressure", "Surface pressure", "atmospheric_pressure_pa", "bar"),
    SearchValue("orbit_distance", "Orbit distance (from the star for a planet, from its planet for a moon)", _BODIES,
                "distance_km", {"planet": "au", "moon": "km"}),
    _body("day_length", "Length of the solar day", "solar_day_hours", "h"),
    _body("rotation_period", "Rotation period", "rotation_period_hours", "h"),
    _body("year_length", "Length of the year (the orbital period)", "period_years", "yr"),
    SearchValue("mass", "Mass", ENTITIES, "mass_kg", {"star": "Msun", "planet": "Mearth", "moon": "Mearth"}),
    _body("surface_dose", "Surface radiation dose", "surface_dose_msv_yr", "mSv/yr"),
    _body("magnetic_moment", "Magnetic moment", "magnetic_moment_a_m2", "A*m2"),
    _body("ocean_depth", "Ocean depth", "ocean_depth_km", "km"),
    SearchValue("moon_count", "Number of moons", ("planet",), None, {"planet": "count"}, {"planet": "count"}, _moon_count),
    SearchValue("luminosity", "Luminosity", ("star",), "luminosity_w", {"star": "Lsun"}),
    SearchValue("age", "Age", ("star",), "age_gy", {"star": "Gyr"}),
    SearchValue("temperature", "Surface temperature", ("star",), "temperature_k", {"star": "K"}),
    *(_gas(gas.formula) for gas in gases.GASES),
)}
"""dict: Key -> `SearchValue`, every numeric key the search takes."""

def _info(unit):
    return column_units.unit_info(unit)


def symbol(unit):
    """What is written after a number in `unit`."""
    if unit == "count":
        return ""
    if unit in unit_registry.UNITS:
        return unit_registry.UNITS[unit].symbol
    return unit


def convert(value, from_unit, to_unit):
    """`value` in `from_unit` as a number in `to_unit` (both unit keys of the unit registry or of `column_units.LOCAL_UNITS`).

    Raises:
        ValueError: The two measure different things, or one has no factor.
    """
    if from_unit == to_unit:
        return float(value)
    if from_unit in unit_registry.UNITS and to_unit in unit_registry.UNITS:
        return unit_registry.convert(value, from_unit, to_unit)
    source, target = _info(from_unit), _info(to_unit)
    if source.kind != target.kind or source.factor is None or target.factor is None:
        raise ValueError(f"{from_unit!r} cannot be turned into {to_unit!r}")
    return float(value) * source.factor / target.factor


def _unit_choice(value, entity, name, requested):
    """The unit key the request named, checked to measure what the key measures; the key's own when none was named."""
    default = value.default_unit(entity)
    if not requested:
        return default
    if requested != default and requested not in unit_registry.UNITS and requested not in column_units.LOCAL_UNITS:
        raise ValueError(f"{name} names the unit {requested!r}, which is not a unit")
    try:
        convert(1.0, requested, default)
    except (ValueError, KeyError, TypeError):
        raise ValueError(f"{name} must be a unit of {_info(default).kind} (the default is {default!r}), got {requested!r}") from None
    return requested


COMMON_UNITS = {
    "acceleration": ("g0", "m_s2"),
    "temperature": ("K", "degC", "degF"),
    "pressure": ("bar", "mbar", "kPa", "Pa", "atm", "psi"),
    "length": ("km", "m", "au", "Rearth", "Rjup", "Rsun", "mi"),
    "duration": ("h", "d", "yr", "Gyr", "min", "s"),
    "mass": ("Mearth", "Mjup", "Mmoon", "Msun", "kg"),
}
"""dict: The units a form offers for a kind of quantity beside the key's own (the unit registry holds many more, and any of them is taken by the API)."""


def unit_choices(entity, key):
    """The unit keys a form offers for `key` of `entity`: the default first, then the common units of its kind that the registry holds."""
    default = KEYS[key].default_unit(entity)
    if default not in unit_registry.UNITS:
        return (default,)
    kind = unit_registry.UNITS[default].kind
    others = [unit for unit in COMMON_UNITS.get(kind, ()) if unit != default and unit in unit_registry.UNITS
              and unit_registry.UNITS[unit].kind == kind]
    return (default, *others)


def _bound(args, name):
    raw = args.get(name)
    if raw is None or raw.strip() == "":
        return None
    try:
        number = float(raw)
    except ValueError:
        raise ValueError(f"{name} must be a number, got {raw!r}") from None
    if not math.isfinite(number):
        raise ValueError(f"{name} must be a finite number, got {raw!r}")
    return number


def ranges_from_args(args):
    """
    The numeric ranges of a request.

    Args:
        args (Mapping): The query parameters (`request.args`).

    Returns:
        dict: `{entity: {key: (least, most, unit)}}`, each bound `None` when left out and `unit` the unit key both are in; an
            entity with no range has no entry.

    Raises:
        ValueError: A bound is not a finite number, the least is over the most, or a unit is not one of the quantity.
    """
    found = {}
    for entity in ENTITIES:
        for key, value in KEYS.items():
            if entity not in value.entities:
                continue
            least = _bound(args, f"{entity}_min_{key}")
            most = _bound(args, f"{entity}_max_{key}")
            if least is None and most is None:
                continue
            if least is not None and most is not None and least > most:
                raise ValueError(f"{entity}_min_{key} must not exceed {entity}_max_{key}")
            unit = _unit_choice(value, entity, f"{entity}_{key}_unit", args.get(f"{entity}_{key}_unit"))
            found.setdefault(entity, {})[key] = (least, most, unit)
    return found


def sorts_in(sort, entity):
    """`(SearchValue, descending)` when `sort` (`gravity` or `-gravity`) names a key of `entity`, else `None`."""
    key = (sort or "").lstrip("-")
    value = KEYS.get(key)
    if value is None or entity not in value.entities:
        return None
    return value, (sort or "").startswith("-")


def clauses(entity, alias, ranges, out_clauses, out_params):
    """Adds one condition per bound of `ranges` (`{key: (least, most, unit)}` of the entity) to `out_clauses` and `out_params`, the
    bounds turned into the unit the quantity is stored in."""
    for key, (least, most, unit) in (ranges or {}).items():
        value = KEYS[key]
        stored = value.stored_unit(entity)
        sql = value.sql(entity, alias)
        if least is not None:
            out_clauses.append(f"{sql} >= ?")
            out_params.append(convert(least, unit, stored))
        if most is not None:
            out_clauses.append(f"{sql} <= ?")
            out_params.append(convert(most, unit, stored))


def order_by(entity, alias, sort, default):
    """The ORDER BY of a panel: by the key `sort` names (smallest first, or with `-` largest first; no value last), then `default`.
    Any other `sort` keeps `default`."""
    found = sorts_in(sort, entity)
    if found is None:
        return default
    value, descending = found
    sql = value.sql(entity, alias)
    return f"{sql} IS NULL, {sql} {'DESC' if descending else 'ASC'}, {default}"


def used(entity, ranges, sort):
    """The keys whose value the rows of the panel carry: those with a range and the one sorted by, as `[(key, unit)]`."""
    keys = {key: unit for key, (_l, _m, unit) in (ranges or {}).items()}
    found = sorts_in(sort, entity)
    if found is not None:
        keys.setdefault(found[0].key, found[0].default_unit(entity))
    return list(keys.items())


def select_columns(entity, alias, keys):
    """The SELECT columns (`, expr AS value_<n>`) that read the value of each `(key, unit)` of `keys`, to put after the panel's own."""
    return "".join(f", {KEYS[key].sql(entity, alias)} AS value_{number}" for number, (key, _unit) in enumerate(keys))


def row_values(entity, row, keys):
    """The `values` of a result row: `{key: {"value", "unit", "symbol"}}` for each `(key, unit)` of `keys`, the stored value
    turned into that unit (`None` when the body has none); the `value_<n>` columns are taken out of `row`."""
    values = {}
    for number, (key, unit) in enumerate(keys):
        raw = row.pop(f"value_{number}", None)
        stored = KEYS[key].stored_unit(entity)
        values[key] = {"value": None if raw is None else convert(raw, stored, unit), "unit": unit, "symbol": symbol(unit)}
    return values


def describe():
    """The registry as the page reads it: a list of `{key, label, entities, units}`, `units` the default unit key of each entity
    with its symbol and the kind the other units must be of."""
    described = []
    for key, value in KEYS.items():
        units = {entity: {"unit": value.default_unit(entity), "symbol": symbol(value.default_unit(entity)),
                          "kind": _info(value.default_unit(entity)).kind} for entity in value.entities}
        described.append({"key": key, "label": value.label, "entities": list(value.entities), "units": units})
    return described


def active(ranges):
    """The ranges of a request as a flat list for a filter chip: `{entity, key, label, least, most, unit, symbol}`."""
    flat = []
    for entity, by_key in (ranges or {}).items():
        for key, (least, most, unit) in by_key.items():
            flat.append({"entity": entity, "key": key, "label": KEYS[key].label, "least": least, "most": most,
                         "unit": unit, "symbol": symbol(unit)})
    return flat
