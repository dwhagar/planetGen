# tests/test_stored_outputs.py

"""
DB.101: every public number a generated planet or moon carries has a home in the database, or says why it has none.

The wiring audit of 2026-10-24 (reports/wiring-db-ux-audit-2026-10-24.md) asked whether each physics output reaches a column. `tests/test_db_column_round_trip.py` holds the other direction (every column is
written and read back); this holds this one: a new public attribute on a generated body that no `object_rows` field of the body tables loads, and that `HOMES` does not place elsewhere, fails here. The
fix is a column, its `db/column_units.py` entry and a migration shipped with the item that makes the value (the standing rule), or a line in `HOMES` that says where the value is kept or why it is not.
"""

import collections

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

HOMES = {
    # kept in another table
    "deposit_types": "`planet_deposit_types` and `moon_deposit_types`, the child rows `store._body_child_map` loads",
    "interior_layers": "`planet_interior_layers` and `moon_interior_layers` (DB.89)",
    "life_class_scores": "`planet_life_class_scores` and `moon_life_class_scores`",
    "life_living": "`planet_life` and `moon_life`, the living species a body carries",
    "reflection_spectrum_visible": "`planet_reflection_spectrum` and `moon_reflection_spectrum`",
    "reflection_spectrum_non_visible": "`planet_reflection_spectrum` and `moon_reflection_spectrum`",
    "isotopes": "`planet_isotopes` and `moon_isotopes`",
    "gas_partials_pa": "`planet_gases` and `moon_gases` (DB.71): the air as partial pressures",
    "moons": "the `moons` table, one row a moon with its planet's row id as parent",
    "formation_distance_au": "`body_history.formation_distance_au` (GEN.1753): where the body formed, which keys its end-stage draw streams",
    # not a column of the body: the row says it, or the parent does
    "is_moon": "which table the row is in (`planets` or `moons`)",
    "host_eccentricity": "not stored (DB.123): a moon is told its parent's eccentricity when the system is loaded",
    "primary_mass_kg": "not stored: what the body orbits, read from its parent when the system is loaded",
    "moon_disk_system": "not stored: the circumplanetary disk of a disk moon's giant (`generation/planet.py`), drawn again with the moons",
    "fallen_moons": "not stored: the moons the history took away, worked out again by `moon_history.apply` from the born moons (`lost_moons`, GEN.1753); the ring they fed is stored (`ring_*` columns)",
    "phi_bio": "not stored: the likelihood the microbial life score is made from (`life_scores`); the score is stored (`planet_life_class_scores`, migration 0150 dropped the likelihood columns)",
    "phi_cpx": "not stored: the likelihood the complex life score is made from (`life_scores`); the score is stored (`planet_life_class_scores`)",
    "oddity_quantities": "not stored: the numbers an oddity row gave the body (`physics/oddities.py`) and the eccentricity a resonant chain forces on it (`resonant_chains`, GEN.1459); the oddity ids are stored (`oddities`)",
    "belt_member": "not stored: the flag of a body drawn as a belt member; a promoted body keeps its belt in `belt_order` and `belt_rank`",
}
"""dict: Public attribute of a generated planet or moon -> where it is kept when it is not loaded by a field of the body tables (`store.BODY_FIELDS`, `PLANET_ONLY_FIELDS`), or why it is not kept."""


def _loaded():
    """The attribute names the field lists of the planet and moon rows load back."""
    class Row(collections.defaultdict):
        def __missing__(self, key):
            return 1.0

    names = set()
    for field in (*store.PLANET_FIELDS, *store.MOON_FIELDS):
        if field.load is None:
            continue
        try:
            names |= set(field.load(Row()))
        except Exception:      # noqa: BLE001 -- a loader that wants a real row (a tri-state flag) is read by its columns below
            names |= set(field.columns)
    return names


def _systems(count=5):
    for seed in range(60):
        config = SystemConfig()
        config.PLANETS = True
        config.MAX_PLANETS = True
        config.MOONS = True
        with draw.bound(seed):
            system = StarSystem(system_config=config)
        yield system
        count -= 1
        if count == 0:
            return


def _public(body):
    return {name for name, value in vars(body).items()
            if not name.startswith("_") and (value is None or isinstance(value, (int, float, str, bool, dict, list, tuple)))}


def test_every_public_number_of_a_generated_body_is_loaded_from_a_column_or_has_a_home():
    loaded = _loaded()
    found = collections.defaultdict(set)
    for system in _systems():
        for planet in system.planets:
            if planet.body_type == "a":
                continue
            for body in (planet, *planet.moons):
                for name in _public(body):
                    found[name].add(body.name)
    homeless = sorted(name for name in found if name not in loaded and name not in HOMES)
    assert not homeless, (f"public attributes of a generated body with no column that loads them and no line in HOMES (DB.101): {homeless}; add the column, its column_units entry and a migration "
                          "in the item's PR, or say in HOMES where it is kept or why it is not")


def test_a_home_that_the_body_tables_load_anyway_or_that_no_body_has_is_taken_out_of_the_list():
    loaded = _loaded()
    found = set()
    for system in _systems():
        for planet in system.planets:
            if planet.body_type == "a":
                continue
            for body in (planet, *planet.moons):
                found |= _public(body)
    stale = sorted(name for name in HOMES if name in loaded or name not in found)
    assert not stale, f"HOMES lists attributes the body tables load or no generated body has: {stale}"


STAR_HOMES = {
    "galactic_center_dist_ly": "not stored: worked out from the stored galactic position when the star is loaded",
    "galactic_orbital_speed_kms": "not stored: worked out from the stored galactic velocity; the orbit's period and phase are stored (`galactic_orbital_period_gy`, `galactic_orbital_phase_deg`)",
}
"""dict: Public attribute of a generated star -> where it is kept when no field of the `stars` row loads it, or why it is not kept."""


def _star_loaded():
    class Row(collections.defaultdict):
        def __missing__(self, key):
            return 1.0

    names = set()
    for field in store.STAR_FIELDS:
        if field.load is None:
            continue
        try:
            names |= set(field.load(Row()))
        except Exception:      # noqa: BLE001 -- as `_loaded`
            names |= set(field.columns)
    return names


def test_every_public_number_of_a_generated_star_is_loaded_from_a_column_or_has_a_home():
    loaded = _star_loaded()
    found = collections.defaultdict(set)
    for system in _systems():
        for star in system.stars:
            for name in _public(star):
                found[name].add(star.name)
    homeless = sorted(name for name in found if name not in loaded and name not in STAR_HOMES)
    assert not homeless, f"public attributes of a generated star with no column that loads them and no line in STAR_HOMES (DB.101): {homeless}"
