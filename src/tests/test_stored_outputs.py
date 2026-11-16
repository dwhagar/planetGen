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
    "oddity_quantities": (
        "stored in part: the eccentricity a chain forces on a planet is `body_history.forced_eccentricity` (DB.132); `tilt_chaotic` (GEN.1535) is not stored: it is the gate that oddities W6 and W61 are drawn by when the "
        "tilt is made, read at generation only; the oddity ids and the tilt they led to are stored (`oddities`, `axial_tilt_deg`, `obliquity_min_deg`, `obliquity_max_deg`)"),
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
    "initial_mass_sol": "not stored on `stars` (`bright_stars` keeps its own): the mass the star was born with, from the population model; a restored star has None and `remnant_orbits` checks for that",
    "phase_end_age_gy": "not stored on `stars` (`bright_stars` keeps its own): the age the drawn evolutionary phase ends at; a restored star has None",
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


# --- DB.134: the same check for belts, comets, the objects that stand alone, and the private attributes -------------------------------------------------------------------------

NAME_GIVEN = "not stored: only says whether the name was given by hand, which decides whether it is reserved when the object is first saved (`store._take_name`); the name is stored"
DERIVED_SPEED = "not stored: worked out from the stored galactic velocity when the object is loaded; the orbit's period and phase are stored (`galactic_orbital_period_gy`, `galactic_orbital_phase_deg`)"

KIND_HOMES = {
    "asteroid belts": {
        "body_type": "which table the row is in (`asteroid_belts`): a belt is always 'a'",
        "composition": "`asteroid_belt_composition`, the child rows of the belt",
        "cleared_spans": "not stored: the stretches of the belt the giants cleared (GEN.1092), used only to place the members; a belt read back keeps the members it was saved with",
        "scatter_giant_au": "not stored: the distance of the giant that swept an outer disk (GEN.1092), used only to place the members; the members it placed are stored",
        "oddity_quantities": "not stored: the numbers an oddity row gave the belt, read when the oddity is drawn at generation; the oddity ids are stored (`oddities`)",
    },
    "comets": {
        "composition": "`comet_composition`, the child rows of the comet",
        "name_given": NAME_GIVEN,
    },
    "rogue planets": {"galactic_orbital_speed_kms": DERIVED_SPEED, "name_given": NAME_GIVEN},
    "interstellar comets": {
        "composition": "`interstellar_comet_composition`, the child rows of the comet",
        "galactic_orbital_speed_kms": DERIVED_SPEED,
        "name_given": NAME_GIVEN,
    },
    "asteroid fields": {
        "composition": "`asteroid_field_composition`, the child rows of the field",
        "galactic_orbital_speed_kms": DERIVED_SPEED,
        "name_given": NAME_GIVEN,
    },
    "black holes": {
        "galactic_center_dist_ly": "not stored: worked out from the stored galactic position when the object is loaded",
        "galactic_orbital_speed_kms": DERIVED_SPEED,
        "name_given": NAME_GIVEN,
    },
    "neutron stars": {
        "galactic_center_dist_ly": "not stored: worked out from the stored galactic position when the object is loaded",
        "galactic_orbital_speed_kms": DERIVED_SPEED,
        "name_given": NAME_GIVEN,
    },
}
"""dict: Kind -> {public attribute of a generated object -> where it is kept when no field of the kind's rows loads it, or why it is not kept}."""

KIND_FIELDS = {
    "asteroid belts": store.BELT_FIELDS,
    "comets": store.COMET_FIELDS,
    "rogue planets": (*store.ROGUE_PLANET_FIELDS, *store.ROGUE_BODY_ROW_FIELDS),
    "interstellar comets": store.INTERSTELLAR_COMET_FIELDS,
    "asteroid fields": store.ASTEROID_FIELD_FIELDS,
    "black holes": (*store.STAR_FIELDS, *store.BLACK_HOLE_FIELDS),
    "neutron stars": (*store.STAR_FIELDS, *store.NEUTRON_STAR_FIELDS),
}
"""dict: Kind -> the `Field`s of the rows that hold it (a remnant's `stars` row and its own)."""

PRIVATE_HOMES = {
    # names and the placement of an object before it has its motion
    "_name": "not stored: the name drawn on first read (`LazySaladName`, PERF.43); the name is stored",
    "_name_seed": "not stored: the seed the name is drawn from; the name is stored",
    "_staged_au": "not stored: a position handed to an object before its motion is set (`mass_object`); the position is stored",
    "_staged_velocity_kms": "not stored: a velocity handed to an object before its motion is set (`mass_object`); the velocity is stored",
    "_build": "`BuildRecord` (GEN.1651), the build-time state of a body, one object for what were twenty private attributes: `orbit_birth`, `birth_spin`, `born_moons` and `resonance_e` are kept in `body_history` and `lost_moons` and attached on load (`body_history.attach`, `moon_history.restore_born`); the rest (the air's draws and light, the CO2 sink, the retention cap, the heat at the surface, the envelope and bare core, the late core, the Venus-type spin, the frame axis, the deferred climate, the Kozai outcome, the late accretion) is working state that is not stored: what it produced is",
    "_model_age_lifespan": "not stored: the age and lifespan the star model drew, which `Star._generate_from_population_model` sets once; `age_gy` and `lifespan_gy` are stored",
    "_model_fe_h": "not stored: the metallicity the star model drew; `metallicity_fe_h` is stored",
    "_model_rotation_hours": "not stored: the rotation the star model drew; `rotation_period_hours` is stored",
    # kept, in `body_history` and `lost_moons` (GEN.1497, GEN.1533, GEN.1753, DB.132)
    "_birth_distance_au": "`body_history.birth_distance_au` (a moon's orbit when it was born), attached on load",
    # the air and the tide state of a body, kept for the passes that make it again at the star's final age (the area of GEN.1745)
    "_air_choice": "not stored: the branch of the air the body took (`air_branches`), kept for the passes that make the air again at the final age; the air is stored (`planet_gases`)",
    "_air_context": "not stored: the setting the air was solved in (`air_final`), kept for the life stage; the air is stored",
    "_air_haze": "not stored: the haze of the air kept for the life stage (GEN.1589); the air is stored",
    "_air_stale": "not stored: marks air made on an orbit that has since moved, which the final-age refresh makes again; a loaded system is final",
    "_tide_parent": "not stored: the planet that tides a moon, set for the system by `tide_sources.register` (the parent is stored)",
    "_tide_siblings": "not stored: the other bodies that pull on a body, set for the system by `tide_sources.register` (their orbits are stored)",
    "_tidal_radius_km": "not stored: the radius a giant's envelope is tidally cut at (`apply_envelope`); the radius is stored",
    "_exchange_damped_yr": "not stored: the age the secular exchange wore the eccentricity down (GEN.1459), which `settle_orbit` reads at generation; the eccentricity is stored",
    # working state of the generation of a body
    "_envelope_helium_y": "not stored: the helium share a hydrogen-first loss left (TEST.761); the composition is stored",
    "_formation_fractions": "not stored: the volatile fractions the body formed with (`envelope_share`), read while its envelope is made",
    "_giant_draw": "not stored: the draw a giant's radius was settled with (`giant_draw`); the radius is stored",
    "_belt_exposure_msv_yr": "not stored: the dose a moon takes from its planet's radiation belts, kept for the dose pass; the dose is stored",
    "_motion_done": "not stored: whether `finish_motion` has made the body's position and velocity (GEN.1548); a loaded body has them",
    "_companions": "not stored: the companion stars the light of a body is summed over (`radiation_field`), found again from the system",
    "_members": "not stored: the members of a belt drawn on demand (`belt_ondemand`); a belt read back keeps the members it was saved with",
    "_secular_planets": "not stored: the axes and masses of the planets beside a belt (GEN.1492), which dilute a companion star's forcing of its members; the planets are stored",
}
"""dict: Private (underscore) attribute of a generated object -> why it is not kept. Working state is not a stored output; what it produced is. A new one must say so here."""

def _loaded_by(fields):
    """The attribute names the `fields` load back, and the columns of a field that is only written (its column is named as the attribute, so the number has a column)."""
    class Row(collections.defaultdict):
        def __missing__(self, key):
            return 1.0

    names = set()
    for field in fields:
        if field.load is None:
            names |= set(field.columns)
            continue
        try:
            names |= set(field.load(Row()))
        except Exception:      # noqa: BLE001 -- as `_loaded`
            names |= set(field.columns)
    return names


def _kinds():
    """`{kind: [generated objects]}` for the kinds of `KIND_FIELDS` (made once for the module): the belts and comets of five systems and five of each stand-alone kind."""
    from planetgen.generation.phenomena.asteroid_field import AsteroidField
    from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
    from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet

    if _KINDS:
        return _KINDS
    found = {kind: [] for kind in KIND_FIELDS}
    for seed in range(60):
        config = SystemConfig()
        config.PLANETS = config.MAX_PLANETS = config.MOONS = config.COMETS = True
        with draw.bound(seed):
            system = StarSystem(system_config=config)
        found["asteroid belts"] += [p for p in system.planets if p.body_type == "a"]
        found["comets"] += list(system.comets)
        if seed >= 4:
            break
    for seed in range(5):
        for kind, cls in (("rogue planets", RoguePlanet), ("interstellar comets", InterstellarComet), ("asteroid fields", AsteroidField), ("black holes", BlackHole), ("neutron stars", NeutronStar)):
            with draw.bound(seed):
                found[kind].append(cls(SystemConfig()))
    _KINDS.update(found)
    return _KINDS


_KINDS = {}


@pytest.mark.parametrize("kind", sorted(KIND_FIELDS))
def test_every_public_number_of_a_generated_object_of_each_other_kind_is_loaded_from_a_column_or_has_a_home(kind):
    objects = _kinds()[kind]
    assert objects, f"no generated {kind} to check"
    loaded = _loaded_by(KIND_FIELDS[kind])
    found = collections.defaultdict(set)
    for obj in objects:
        for name in _public(obj):
            found[name].add(getattr(obj, "name", kind))
    homeless = sorted(name for name in found if name not in loaded and name not in KIND_HOMES[kind])
    stale = sorted(name for name in KIND_HOMES[kind] if name in loaded or name not in found)
    assert not homeless, f"public attributes of a generated object of the kind {kind!r} with no column that loads them and no line in KIND_HOMES (DB.134): {homeless}"
    assert not stale, f"KIND_HOMES lists attributes of {kind!r} that its rows load or no generated object has: {stale}"


def test_every_private_attribute_of_a_generated_body_star_or_other_object_says_why_it_is_not_kept():
    found = collections.defaultdict(set)
    for system in _systems():
        for obj in (*system.planets, *(m for p in system.planets if p.body_type != "a" for m in p.moons), *system.stars, *system.comets):
            for name in vars(obj):
                if name.startswith("_"):
                    found[name].add(getattr(obj, "name", type(obj).__name__))
    for objects in _kinds().values():
        for obj in objects:
            for name in vars(obj):
                if name.startswith("_"):
                    found[name].add(getattr(obj, "name", type(obj).__name__))
    homeless = sorted(name for name in found if name not in PRIVATE_HOMES)
    assert not homeless, f"private attributes of a generated object with no line in PRIVATE_HOMES (DB.134): {homeless}; a stored output belongs in a column, working state says here why it is not kept"
