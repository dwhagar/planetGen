# planetgen/physics/oddities.py

"""
The oddity table: rare, physically plausible outcomes that enter the generator as a cause at one input (GEN.1036, docs/design/systems/weird-worlds.md)
===================================================================================================================================================

An oddity is never painted on a finished planet. Each row of `physics/data/oddities.json` is a **cause** that changes one input of the chain (a
lever of design section 8.1: the spin axis and rate, the orbit, the bulk composition, the interior heat, the history, the moon and ring set, or an
air state that has two branches); the solvers already built then make the surface. A row holds:

- `id`, `name`, `kind` (planet, moon, giant, rogue, small body or star system), and the `lever` and `input` it changes (the input must belong to
  the lever, so a row can only change what a lever allows);
- `requires`: the gate, conditions on stored quantities in the shape of `planet_classes.json` (`quantity` with `min`, `lt`, `equals` or
  `not_equals`), never empty: a rate is a rate on a gated population, not a global percentage;
- `rate`: the chance given the gate, `value` with its `source`, `label` (READ in the source, 2ND-HAND in a paper that quotes it, RECALLED), GEN.725
  `tier` and `tweak` level and the `registry` id of its stand-in row;
- `apply`: the named `effect` with its `params`, which the module that implements the row defines;
- `limits`: the conditions that must hold on the body after the effect (the plausibility harness checks them), never empty;
- `source` and `text`: where it is read and the one sentence the page shows ("This world turns backwards because ...").

The loader refuses a row without any of these, a quantity the file does not declare, an input outside its lever, a rate outside 0 to 1, a repeated
id, or a tier or tweak that is not in the registry's scale. An empty table loads.

**Draws.** Every roll uses a stream of its own, keyed on the oddity id and the body (`stream`), never the run's stream, so adding, removing or
turning off a row changes no other draw (design 4.4; the isolation test holds this). `probability` multiplies the row's rate by the world's
weirdness dial (`tuning.WEIRDNESS_DIAL`, 0 turns every oddity off) and the row's own scale (`tuning.WEIRDNESS_RATE_SCALE`).

**The wiring (GEN.1047).** `apply(body, levers)` is what the generator calls at the point of the chain where a lever's input is settled: it asks
for the rows of the body's kind whose gate holds (`quantities_of`), rolls each on its own stream keyed by `body_key`, runs the effect the row
names (`EFFECTS`, filled by the module that implements the effect through `register`), and records the ids on `body.oddities`. A row an
admin or a request forces (`SystemConfig.FORCE_ODDITIES`) skips the roll but not the gate: a world cannot be given what its conditions
forbid. `body_key` is made from the star and the body's own orbit and mass, so an empty table, or a row the dial turns off, draws nothing
from the run's stream and the golden seed does not move. Use:

    applied = apply(body, levers=("L3",))             # at the bulk step: the ids that were given to the body
    broken = limit_failures(row, quantities_of(body))  # [] when the body is still physical
"""

import json
import os
import re
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import constants

EFFECTS = {}
"""dict: The effect functions by name, `function(body, row) -> None`, which change the one input of the row's lever on the body; registered by the module
that implements the effect (`register`)."""


def register(effect):
    """Decorator: `@register("strip_mantle")` makes the function the effect of every row that names it."""
    def wrap(function):
        if effect in EFFECTS and EFFECTS[effect] is not function:
            raise OddityError(f"the effect {effect!r} is registered twice")
        EFFECTS[effect] = function
        return function
    return wrap


# (`EFFECTS` and `register` come before the imports below on purpose: a module that registers an effect at import time, such as `rings`, can be reached
# through `planet_classes` while this module is still loading, and needs them to exist by then.)

from planetgen.physics import elements, planet_classes
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "oddities.json")
"""str: The data file."""

LABELS = ("READ", "2ND-HAND", "RECALLED")
"""tuple: How far a rate was checked: read in the source, read in a paper that quotes it, or recalled and not opened (a stand-in until confirmed)."""

TIERS = ("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9")
TWEAKS = ("L0", "L1", "L2", "L3", "L4")
CONDITION_KEYS = {"quantity", "min", "lt", "equals", "not_equals"}
_ID = re.compile(r"^[A-Z]{1,2}[0-9]{1,3}$")
_NAME = re.compile(r"^[a-z][a-z0-9_]*$")


class OddityError(DataError):
    """An oddity file that is refused, or a body the table cannot be asked about."""


@dataclass(frozen=True)
class Rate:
    """The chance of an oddity given its gate, and where the number comes from."""
    value: float
    source: str
    label: str
    tier: str
    tweak: str
    registry: str


@dataclass(frozen=True)
class Oddity:
    """One row of the table."""
    id: str
    name: str
    kind: str
    lever: str
    input: str
    effect: str
    params: dict
    requires: tuple
    limits: tuple
    rate: Rate
    source: str
    text: str


def _text(row, field, owner):
    value = row.get(field)
    if not isinstance(value, str) or not value.strip():
        raise OddityError(f"{owner}: '{field}' is missing or empty")
    return value


def _conditions(row, field, owner, quantities):
    found = row.get(field)
    if not isinstance(found, list) or not found:
        raise OddityError(f"{owner}: '{field}' needs at least one condition")
    for condition in found:
        if not isinstance(condition, dict) or "quantity" not in condition or not (set(condition) - {"quantity"}):
            raise OddityError(f"{owner}: {field} condition {condition!r} is not a quantity with min, lt, equals or not_equals")
        if set(condition) - CONDITION_KEYS:
            raise OddityError(f"{owner}: {field} condition {condition!r} uses a key outside {sorted(CONDITION_KEYS)}")
        if condition["quantity"] not in quantities:
            raise OddityError(f"{owner}: {field} reads {condition['quantity']!r}, which the file does not declare")
    return tuple(found)


def _rate(row, owner):
    rate = row.get("rate")
    if not isinstance(rate, dict):
        raise OddityError(f"{owner}: 'rate' is missing")
    value = rate.get("value")
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0.0 <= value <= 1.0:
        raise OddityError(f"{owner}: rate.value {value!r} is not a number from 0 to 1")
    if rate.get("label") not in LABELS:
        raise OddityError(f"{owner}: rate.label {rate.get('label')!r} is not one of {LABELS}")
    if rate.get("tier") not in TIERS:
        raise OddityError(f"{owner}: rate.tier {rate.get('tier')!r} is not one of T1 to T9")
    if rate.get("tweak") not in TWEAKS:
        raise OddityError(f"{owner}: rate.tweak {rate.get('tweak')!r} is not one of L0 to L4")
    return Rate(float(value), _text(rate, "source", f"{owner} rate"), rate["label"], rate["tier"], rate["tweak"],
                _text(rate, "registry", f"{owner} rate"))


def parse(document):
    """
    The table of a loaded document: `(rows, levers, quantities)`.

    Raises:
        OddityError: For any row or vocabulary entry the schema refuses (see the module docstring).
    """
    kinds = tuple(document.get("kinds", ()))
    levers = document.get("levers")
    quantities = document.get("quantities")
    if not kinds or not isinstance(levers, dict) or not levers or not isinstance(quantities, dict):
        raise OddityError("the file needs 'kinds', 'levers' and 'quantities'")
    for key, lever in levers.items():
        if not lever.get("inputs") or not lever.get("label"):
            raise OddityError(f"lever {key}: needs a label and inputs")
    rows = []
    seen = set()
    for row in document.get("rows", ()):
        owner = f"oddity {row.get('id')!r}"
        identifier = row.get("id")
        if not isinstance(identifier, str) or not _ID.match(identifier):
            raise OddityError(f"{owner}: the id is not a letter or two and a number (W10, G1, R4)")
        if identifier in seen:
            raise OddityError(f"{owner}: the id is used by two rows")
        seen.add(identifier)
        if row.get("kind") not in kinds:
            raise OddityError(f"{owner}: kind {row.get('kind')!r} is not one of {kinds}")
        apply = row.get("apply")
        if not isinstance(apply, dict):
            raise OddityError(f"{owner}: 'apply' is missing")
        lever = apply.get("lever")
        if lever not in levers:
            raise OddityError(f"{owner}: apply.lever {lever!r} is not one of {sorted(levers)}")
        if apply.get("input") not in levers[lever]["inputs"]:
            raise OddityError(f"{owner}: apply.input {apply.get('input')!r} is not an input of lever {lever} ({levers[lever]['inputs']})")
        effect = apply.get("effect")
        if not isinstance(effect, str) or not _NAME.match(effect):
            raise OddityError(f"{owner}: apply.effect {effect!r} is not a lower-case name")
        params = apply.get("params", {})
        if not isinstance(params, dict):
            raise OddityError(f"{owner}: apply.params is not an object")
        rows.append(Oddity(
            identifier, _text(row, "name", owner), row["kind"], lever, apply["input"], effect, dict(params),
            _conditions(row, "requires", owner, quantities), _conditions(row, "limits", owner, quantities),
            _rate(row, owner), _text(row, "source", owner), _text(row, "text", owner)))
    return tuple(rows), levers, quantities


def load(path=DATA_FILE):
    """The table in the file at `path`: `(rows, levers, quantities)` (`parse` of its JSON)."""
    with open(path, encoding="utf-8") as handle:
        return parse(json.load(handle))


ROWS, LEVERS, QUANTITIES = load()
"""tuple, dict, dict: The shipped table's `Oddity` rows, the levers (each with its inputs) and the quantities a gate or a limit may read."""

BY_ID = {row.id: row for row in ROWS}


def candidates(kind, found, rows=None):
    """The rows of `kind` whose gate holds for the named quantities in `found` (a quantity the body lacks fails the gate)."""
    return [row for row in (ROWS if rows is None else rows)
            if row.kind == kind and planet_classes.holds(row.requires, found, owner=f"oddity {row.id}", missing=None)]


def limit_failures(row, found):
    """The limits of `row` that do not hold on the body after the effect, as the conditions themselves ([] when it is still physical). A limit
    on a quantity the body lacks counts as failed: an oddity must leave every quantity its limits read."""
    return [limit for limit in row.limits if not planet_classes.holds([limit], found, owner=f"oddity {row.id}", missing=None)]


def probability(row, dial=None, scale=None):
    """The chance that a body that passes the gate gets `row`: its rate times the weirdness dial and the row's own scale, held to 0..1.
    `dial` and `scale` (a dict by id) default to `tuning.WEIRDNESS_DIAL` and `tuning.WEIRDNESS_RATE_SCALE`."""
    dial = tuning.WEIRDNESS_DIAL if dial is None else dial
    scale = tuning.WEIRDNESS_RATE_SCALE if scale is None else scale
    return min(1.0, max(0.0, row.rate.value * dial * scale.get(row.id, 1.0)))


def site_dial():
    """The dial the settings file or page sets (`generation.weirdness_dial`), or `None` when it names none (the tuning's applies). Read live: the settings are re-read only
    when the file changes. A settings file that does not fit the model leaves the tuning's dial, since generation must not stop for it."""
    from planetgen.util import settings

    try:
        section = settings.get_settings().generation
    except settings.SettingsError:
        return None
    return section.weirdness_dial if "weirdness_dial" in section.model_fields_set else None


def stream(row_id, body_key):
    """The stream of one oddity on one body: seeded by the id and the body's key alone (a string seed, hashed by SHA-512, so it does not depend on
    `PYTHONHASHSEED`), so no other row and no other draw of the run moves when a row is added, removed or turned off."""
    return draw.Stream(f"oddity:{row_id}:{body_key}")


def roll(row, body_key, dial=None, scale=None):
    """Whether the body keyed `body_key` gets `row`. The draw comes from the row's own stream and is made even at dial 0, so the stream a body
    has for a row never depends on the dial."""
    draw_value = stream(row.id, body_key).random()
    return draw_value < probability(row, dial, scale)


def words(row):
    """The sentence the page shows for `row`: 'Odd: <name>. <text>'."""
    return f"Odd: {row.name}. {row.text}"


BOND_ALBEDO_FOR_GATES = 0.3
"""float: The Bond albedo the equilibrium temperature of a gate is taken at (the Earth's, as `atmosphere_retention` and `envelope_escape` take it; stand-in HC1 note)."""

KIND_OF_BODY = {"t": "planet", "g": "giant"}
"""dict: The row `kind` a body is asked about, by its `body_type`; a moon is asked as a moon, a rogue planet as a rogue."""


EFFECT_MODULES = ("planetgen.physics.mantle_stripping", "planetgen.physics.venus_spin", "planetgen.physics.seasonal_swing",
                  "planetgen.physics.resonant_chains", "planetgen.physics.orbit_oddities",
                  "planetgen.physics.disintegrating_planet", "planetgen.physics.rings",
                  "planetgen.physics.air_branches", "planetgen.physics.giant_tilt", "planetgen.physics.hycean_candidate", "planetgen.physics.helium_fractionation", "planetgen.physics.obliquity_variation")
"""tuple: The modules that register an effect (`register`), imported the first time `apply` has a row to run; a module that implements a row's effect is listed here."""


def _import_effects():
    import importlib

    for name in EFFECT_MODULES:
        importlib.import_module(name)


def kind_of(body):
    """The row kind that applies to `body`: moon, rogue, giant or planet."""
    if getattr(body, "is_moon", False):
        return "moon"
    if getattr(body, "is_rogue", False) or getattr(body, "rogue", False):
        return "rogue"
    return KIND_OF_BODY.get(getattr(body, "body_type", None), "planet")


def body_key(body):
    """
    A key for `body` that no other body shares and that needs no draw: the star's key, whether it is a moon, the orbit's distance and the mass. A body made
    again with the same star, orbit and mass gets the same key, so its oddities do not move.
    """
    from planetgen.physics import radiogenic

    star = getattr(body, "star", None)
    star_key = radiogenic._star_key(star) if star is not None else None
    parent = getattr(body, "primary_mass_kg", None) if getattr(body, "is_moon", False) else None
    return f"{star_key!r}:{int(bool(getattr(body, 'is_moon', False)))}:{getattr(body, 'distance', None)!r}:{getattr(body, 'mass', None)!r}:{parent!r}"


def quantities_of(body):
    """
    The quantities of `body` that the table's gates and limits read (`oddities.json` `quantities`), for the ones the body has. A quantity the body
    lacks is left out, and a gate that reads it fails: a row never applies on a guess.
    """
    found = {}

    def put(name, value):
        if value is not None:
            found[name] = value

    body_type = getattr(body, "body_type", None)
    put("body_type", body_type)
    mass, radius = getattr(body, "mass", None), getattr(body, "radius", None)
    put("mass_earth", None if mass is None else mass / constants.EARTH_MASS_TO_KG)
    put("radius_earth", None if radius is None else radius / constants.EARTH_RADIUS_KM)
    put("semi_major_axis_au", getattr(body, "distance", None))
    if mass and radius:       # GEN.1523: the bulk density (`hycean_candidate.density_g_cm3`, the one the box reads) and the equilibrium temperature a gate reads
        from planetgen.physics import hycean_candidate

        put("density_g_cm3", hycean_candidate.density_g_cm3(mass / constants.EARTH_MASS_TO_KG, radius / constants.EARTH_RADIUS_KM))
    put("equilibrium_temperature_k", _equilibrium_temperature_k(body))
    put("envelope_helium_y", getattr(body, "_envelope_helium_y", None))      # GEN.1527: set by `envelope_capture.helium_share` when the loss changed it
    host_radius_km = getattr(getattr(body, "star", None), "radius", None)
    if getattr(body, "distance", None) and host_radius_km and not getattr(body, "is_moon", False):
        put("distance_over_host_radius", body.distance * constants.AU_TO_KM / host_radius_km)
    put("eccentricity", getattr(body, "eccentricity", None))
    period = getattr(body, "period", None)
    put("orbital_period_days", None if period is None else period * 365.25)
    star = getattr(body, "star", None)
    put("age_gyr", getattr(star, "age", None))
    put("host_teff_k", getattr(star, "temperature", None))
    star_mass = getattr(star, "mass", None)
    put("host_mass_sun", None if star_mass is None else star_mass / constants.SOLAR_MASS_TO_KG)
    star_teff, star_radius_km, distance = getattr(star, "temperature", None), getattr(star, "radius", None), getattr(body, "distance", None)
    if star_teff and star_radius_km and distance and not getattr(body, "is_moon", False):
        put("dayside_temperature_k", star_teff * (star_radius_km * 1000.0 / (distance * constants.AU_M)) ** 0.5)     # `disintegrating_planet.dayside_temperature_k`
    if hasattr(body, "ring_origin"):
        put("ring_from_fallen_moon", 1 if body.ring_origin == "disrupted_moon" else 0)
    put("surface_pressure_pa", getattr(body, "atmospheric_pressure", None))
    put("obliquity_deg", getattr(body, "axial_tilt_deg", None))
    hours = getattr(body, "rotation_period_hours", None)
    put("spin_period_days", None if hours is None else abs(hours) / 24.0)
    put("core_mass_fraction", getattr(body, "metal_mass_fraction", None))
    put("water_mass_fraction", getattr(body, "water_mass_fraction", None))
    put("envelope_mass_fraction", getattr(body, "envelope_mass_fraction", None))
    put("collision_count", getattr(body, "collision_count", None))
    moons = getattr(body, "moons", None)
    if moons is not None and mass:
        put("moon_mass_ratio", max((float(moon.mass) for moon in moons if getattr(moon, "mass", None)), default=0.0) / mass)
    for name, value in (getattr(body, "oddity_quantities", None) or {}).items():
        put(name, value)       # quantities the chain keeps with the body for a gate (the architecture's, the moon set's)
    return found


def _equilibrium_temperature_k(body):
    """The equilibrium temperature of a planet at the Bond albedo `BOND_ALBEDO_FOR_GATES` under the light of every star at its distance (`radiation_field.starlight_w_m2`, GEN.1441);
    `None` for a moon, or a body with no star, orbit or luminosity."""
    from planetgen.physics import climate_draws, orbit_elements, radiation_field

    star, distance = getattr(body, "star", None), getattr(body, "distance", None)
    if star is None or not distance or not getattr(star, "luminosity", None) or getattr(body, "is_moon", False):
        return None
    flux = radiation_field.starlight_w_m2(star, orbit_elements.flux_distance_au(distance, getattr(body, "eccentricity", 0.0) or 0.0))
    return climate_draws.equilibrium_temperature_k(BOND_ALBEDO_FOR_GATES, flux)


def ids_of(body):
    """The ids of the oddities `body` has, a tuple (empty for none)."""
    return tuple(getattr(body, "oddities", ()) or ())


def forced_ids(body):
    """The oddity ids the body's system asks for (`SystemConfig.FORCE_ODDITIES`), a tuple."""
    config = getattr(body, "system_config", None)
    return tuple(getattr(config, "FORCE_ODDITIES", None) or ())


def apply(body, levers=None, rows=None, dial=None, scale=None):
    """
    Gives `body` the oddities of the table that belong at this step of the chain.

    Args:
        body: A planet, moon or rogue planet (anything `quantities_of` can read).
        levers (iterable of str, optional): Only rows whose `apply` lever is one of these (`'L3'`...): the step names the levers whose inputs it has settled.
            All rows when omitted.
        rows (iterable of Oddity, optional): The table to use instead of the loaded one (tests).
        dial, scale: The weirdness dial and the per-row scale to use instead of the tuning's.

    Returns:
        list[str]: The ids given at this call (also added to `body.oddities`, in table order).
    """
    table = ROWS if rows is None else tuple(rows)
    if not table:
        return []
    _import_effects()
    if dial is None:
        dial = getattr(getattr(body, "system_config", None), "WEIRDNESS_DIAL", None)       # a run's own dial; the settings file's, then the tuning's, when it has none
        if dial is None:
            dial = site_dial()
    found = quantities_of(body)
    forced = set(forced_ids(body))
    key = body_key(body)
    given = []
    for row in table:
        if row.kind != kind_of(body) or not candidates(row.kind, found, (row,)):      # the gate is read on the body as it is after the rows before this one
            continue
        if levers is not None and row.lever not in levers:
            continue
        if row.id in ids_of(body):
            continue
        wanted = row.id in forced
        if not (wanted or roll(row, key, dial, scale)):
            continue
        if row.effect not in EFFECTS:
            raise OddityError(f"oddity {row.id}: no module registered the effect {row.effect!r}")
        if EFFECTS[row.effect](body, row) is False:      # the effect found the body cannot have it after all (it would break a limit): not given, not recorded
            continue
        body.oddities = ids_of(body) + (row.id,)
        given.append(row.id)
        found = quantities_of(body)         # the effect moved an input: the next gate reads the body as it now is
    return given


def drop(body, levers):
    """
    Takes from `body` the oddities whose lever is in `levers`, so that a step that runs again (the spin is made again when the age moves or a moon is added) rolls
    its rows again on the body as it now is. The roll is the row's own stream on the body's key, so a row whose gate still holds is given again.
    """
    keep = tuple(identifier for identifier in ids_of(body) if identifier not in BY_ID or BY_ID[identifier].lever not in levers)
    if keep != ids_of(body):
        body.oddities = keep


def column_text(ids):
    """The text a body's ids are stored as: comma separated, `None` for none."""
    return ",".join(ids) if ids else None


def from_column(text):
    """The ids of a stored column value, a tuple (empty for `None` and for an empty string)."""
    return tuple(part for part in (text or "").split(",") if part)
