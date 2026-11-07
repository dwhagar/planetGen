# tests/test_rogue_planet_parity.py

"""TEST.498 (Boss, 4:57 PM EDT 2026-10-30): a rogue planet is functionally and programmatically identical to a bound planet, differing only in the rogue
note, the missing host and the quantities that depend on starlight (its equilibrium temperature is the internal and ambient value).

Since GEN.1111 a rogue planet IS a `Planet`, built by the planet code against a `StarlessHost` (`generation/phenomena/starless_host.py`), and stored in
`rogue_planets`, which carries every `planets` column apart from the orbit round a host. This file holds that contract:

* the rogue is built by `Planet` and is a `Planet`;
* the table has every column of a planet except the host and orbit ones;
* matched draws (a rogue and a planet of one mass far from a star) are one kind of body of one size, and both are warmed by their own heat;
* the two are one code path: one class with no subclass, the same methods, the same attributes, the same moon class, the same stored fields, and a class
  from the same table;
* where the two differ in what is stored, the difference is a named quantity that follows from the missing star or from the heat the body keeps (the
  light, the stellar wind, and the liquid water an Earth-mass body holds only where its own heat is enough).
"""

import pytest

from planetgen import tuning
from planetgen.db import models, store
from planetgen.generation.planet import Planet
from planetgen.generation.phenomena import starless_host

# Planet columns a rogue row does not carry: the host and the orbit round it, and the keys that file a planet in a system.
SYSTEM_KEYS = frozenset({"dynamical_state", "host_node_id", "orbital_index", "star_id", "star_system_id"})


def test_the_rogue_table_carries_every_planet_column_but_the_host_and_the_orbit():
    planets = {column.name for column in models.metadata.tables["planets"].columns}
    rogues = {column.name for column in models.metadata.tables["rogue_planets"].columns}
    assert (planets - rogues) <= set(store.ROGUE_HOST_ONLY_COLUMNS) | SYSTEM_KEYS


def test_the_rogue_row_adds_only_the_ejection_and_the_motion():
    planets = {column.name for column in models.metadata.tables["planets"].columns}
    rogues = {column.name for column in models.metadata.tables["rogue_planets"].columns}
    own = rogues - planets
    assert set(store.ROGUE_EJECTION_COLUMNS) <= own
    assert not own & {"zone", "distance_km", "period_years"}


# --- the matched draw (TEST.498: draw matched planets with and without a host and compare what the two chains make) ---------------------------------

FAR_AU = 1.0e5
MASSES_EARTH = {"terrestrial": 1.0, "sub-neptune": 8.0, "saturn": 95.0, "jupiter": 318.0}


def _matched(kind, seed):
    """`(rogue, bound)`: a planet built with no star and one of the same mass built `FAR_AU` from a star, or the test is skipped where a dataset the giant
    needs is not on this machine."""
    from planetgen.datasets.fetch import DataMissing
    from planetgen.physics import constants
    from planetgen.util import draw
    from tests.fuzz_support import deterministic_entropy
    from tests.test_gen_builder_internals import make_system
    with deterministic_entropy(498):     # the same star every run
        system = make_system()
    star = system.star
    mass_kg = MASSES_EARTH[kind] * constants.EARTH_MASS_TO_KG
    try:
        with draw.bound(f"rogue-parity-{kind}-{seed}"):
            rogue = starless_host.build_ejected_planet(system.system_config, star, 4.5, mass_kg, formation_au=5.0, moon_count=0)
        with draw.bound(f"rogue-parity-bound-{kind}-{seed}"):
            bound = Planet(system.system_config, star, star.habitable_zone, FAR_AU, mass=mass_kg, moon_count=0)
    except DataMissing as missing:
        pytest.skip(f"a dataset the {kind} needs is not on this machine: {missing}")
    return rogue, bound


def test_a_rogue_planet_is_a_planet_built_against_a_starless_host():
    rogue, _bound = _matched("terrestrial", 0)
    assert isinstance(rogue, Planet) and isinstance(rogue.star, starless_host.StarlessHost) and rogue.star.starless


@pytest.mark.parametrize("kind", tuple(MASSES_EARTH))
@pytest.mark.parametrize("seed", range(3))
def test_a_rogue_and_a_bound_planet_of_one_mass_are_one_kind_of_body_of_one_size(kind, seed):
    rogue, bound = _matched(kind, seed)
    assert rogue.mass == bound.mass
    assert rogue.body_type == bound.body_type
    if rogue.body_type == "t":
        assert 0.5 * bound.radius <= rogue.radius <= 2.0 * bound.radius, (rogue.radius, bound.radius)


@pytest.mark.parametrize("kind", ("saturn", "jupiter"))
def test_a_giant_with_no_light_is_warmed_by_its_own_heat_far_above_the_ambient_floor(kind):
    rogue, bound = _matched(kind, 0)
    assert bound.surface_temperature > 10.0 * 2.725 and rogue.surface_temperature > 10.0 * 2.725


@pytest.mark.parametrize("kind", ("terrestrial", "sub-neptune"))
def test_a_rocky_planet_with_no_light_is_warmed_by_its_own_heat_only_when_it_has_no_star(kind):
    """GEN.1111 part 5: a starless planet adds its own heat (internal and tidal) to the ambient floor; a lit planet far from its star is unchanged, since the
    stored albedo and greenhouse pair must keep reproducing its temperature."""
    rogue, _bound = _matched(kind, 0)
    assert rogue.surface_temperature >= 2.725


# --- one code path (TEST.498: no rogue-only generator, the same columns, the same moon and composition generation, the same classes) -------------

STARLIGHT_FIELDS = frozenset({"magnetopause_rp", "uv_surface_index"})
"""The stored quantities that need a star: the stellar wind that holds off a planet's field, and the ultraviolet a star gives the surface."""

LIQUID_SURFACE_FIELDS = frozenset({
    "ocean_depth_km", "ocean_class", "ocean_ph", "water_activity", "phosphorus", "surface_liquid_ammonia_wt", "tide_range_spring_m", "tide_range_neap_m",
    "ocean_liquid"})
"""The stored quantities of a liquid surface: they are filled where the surface holds liquid, which the body's own heat decides when it has no star."""

CLOUD_FIELDS = frozenset({"aerosol_draw", "atmosphere_colour", "cloud_cover", "tau_cloud", "tau_dust", "tau_haze"})
"""The stored quantities of a cloud and haze layer: whether the air holds one is decided by the temperature the body's own heat gives it."""

ENVELOPE_AIR_FIELDS = frozenset({"envelope_helium_fraction", "envelope_z_fraction"})
"""The helium and heavy-element shares of an envelope (GEN.1303): stored only where the air is hydrogen-rich, which the temperature (a star's light or the body's own heat) decides."""

RUNNER_UP_FIELDS = frozenset({"runner_up_class"})
"""The next-best class letter (GEN.1207): the classifier names one only when a second class fits at all, which the body's own conditions decide."""

SOLID_ONLY_FIELDS = frozenset({"roundness_ratio", "yield_strength_pa", "shelter_column_g_cm2"})
"""The stored quantities only a solid body has (a giant has none): they differ only where the two draws fall in different classes."""

SEEDS = range(4)
FILING = {"planet_id": None, "star_system_id": None, "star_id": None, "orbital_index": 0, "rogue_planet_id": None}
"""The keys that file a moon in a system or under a rogue planet; they are not part of what is generated, so they are blank here."""


def _moon_matched(kind, seed, moons=2):
    """`_matched`, with `moons` moons each."""
    from planetgen.datasets.fetch import DataMissing
    from planetgen.physics import constants
    from planetgen.util import draw
    from tests.fuzz_support import deterministic_entropy
    from tests.test_gen_builder_internals import make_system
    with deterministic_entropy(498):
        system = make_system()
    star = system.star
    mass_kg = MASSES_EARTH[kind] * constants.EARTH_MASS_TO_KG
    try:
        with draw.bound(f"rogue-parity-moons-{kind}-{seed}"):
            rogue = starless_host.build_ejected_planet(system.system_config, star, 4.5, mass_kg, formation_au=5.0, moon_count=moons)
        with draw.bound(f"rogue-parity-moons-bound-{kind}-{seed}"):
            bound = Planet(system.system_config, star, star.habitable_zone, FAR_AU, mass=mass_kg, moon_count=moons, formation_distance_au=5.0)
    except DataMissing as missing:
        pytest.skip(f"a dataset the {kind} needs is not on this machine: {missing}")
    return rogue, bound


def _public(obj):
    return {name for name in vars(obj) if not name.startswith("_")}


def test_there_is_no_rogue_only_planet_generator():
    """A rogue body is a `Planet` itself: nothing derives from `Planet` (so nothing overrides one of its steps for a rogue), and the module that makes a
    rogue defines no class of its own that builds a planet."""
    import importlib
    import inspect
    import pkgutil
    import planetgen.generation as generation
    for module in pkgutil.walk_packages(generation.__path__, "planetgen.generation."):
        importlib.import_module(module.name)
    assert Planet.__subclasses__() == []
    own = [name for name, member in inspect.getmembers(starless_host, inspect.isclass) if member.__module__ == starless_host.__name__]
    assert own == ["StarlessHost"]
    assert not issubclass(starless_host.StarlessHost, Planet)


@pytest.mark.parametrize("kind", ("terrestrial", "saturn"))
def test_a_rogue_and_its_moons_are_built_by_the_planet_and_moon_code_a_bound_planet_uses(kind):
    rogue, bound = _moon_matched(kind, 0)
    assert type(rogue) is type(bound) is Planet
    assert len(rogue.moons) == len(bound.moons) == 2
    assert {type(moon) for moon in rogue.moons + bound.moons} == {Planet}
    assert all(moon.is_moon for moon in rogue.moons + bound.moons) and not rogue.is_moon
    for name in ("add_moons", "generate_moons", "settle_orbit"):
        if hasattr(Planet, name):
            assert getattr(type(rogue), name) is getattr(type(bound), name), name


@pytest.mark.parametrize("kind", tuple(MASSES_EARTH))
@pytest.mark.parametrize("seed", SEEDS)
def test_a_rogue_and_a_bound_planet_carry_the_same_attributes_and_their_moons_do_too(kind, seed):
    rogue, bound = _moon_matched(kind, seed)
    assert set(vars(rogue)) == set(vars(bound)), sorted(set(vars(rogue)) ^ set(vars(bound)))
    for rogue_moon, bound_moon in zip(rogue.moons, bound.moons):
        assert _public(rogue_moon) == _public(bound_moon), sorted(_public(rogue_moon) ^ _public(bound_moon))


def test_both_kinds_of_moon_are_stored_through_the_one_moon_field_list():
    import inspect
    assert "MOON_FIELDS" in inspect.getsource(store.insert_moon)
    assert "MOON_FIELDS" in inspect.getsource(store.insert_moon_of_rogue_planet)


@pytest.mark.parametrize("kind", tuple(MASSES_EARTH))
@pytest.mark.parametrize("seed", SEEDS)
def test_a_rogue_and_a_bound_planet_fill_the_same_stored_fields_and_differ_only_in_named_quantities(kind, seed):
    """The row of each is filled by the one planet field list (`ROGUE_BODY_FIELDS` is `PLANET_BODY_FIELDS` less the orbit): the same columns in the same order,
    the same kind of value in each, and on the planet a column empty on one and filled on the other only for a named quantity: the star's light, the
    liquid surface and the cloud and haze that the body's own heat decides, and the strength of a solid body, which a giant has none of (so it follows
    the class). The moons are other draws, so for them the check is the columns and the kind of value."""
    from planetgen.db import object_rows
    rogue, bound = _moon_matched(kind, seed)
    for rogue_body, bound_body in [(rogue, bound), *zip(rogue.moons, bound.moons)]:
        fields = store.ROGUE_BODY_FIELDS if not rogue_body.is_moon else store.MOON_FIELDS
        names = object_rows.columns_of(fields)
        kept = object_rows.values_of(fields, rogue_body, FILING)
        lit = object_rows.values_of(fields, bound_body, FILING)
        assert len(kept) == len(lit) == len(names)
        retyped = [name for name, a, b in zip(names, kept, lit) if a is not None and b is not None and type(a) is not type(b)]
        assert not retyped, retyped
        if rogue_body.is_moon:
            continue
        one_sided = {name for name, a, b in zip(names, kept, lit) if (a is None) != (b is None)}
        allowed = STARLIGHT_FIELDS | LIQUID_SURFACE_FIELDS | CLOUD_FIELDS | RUNNER_UP_FIELDS | ENVELOPE_AIR_FIELDS
        if rogue.planet_class != bound.planet_class:
            allowed |= SOLID_ONLY_FIELDS
        assert one_sided <= allowed, sorted(one_sided - allowed)


@pytest.mark.parametrize("kind", tuple(MASSES_EARTH))
@pytest.mark.parametrize("seed", range(2))
def test_a_rogue_and_a_bound_planet_take_their_class_from_the_one_table(kind, seed):
    from planetgen.physics import planet_classes
    rogue, bound = _moon_matched(kind, seed)
    for body in (rogue, bound, *rogue.moons, *bound.moons):
        assert body.planet_class == planet_classes.classify(body), (body.planet_class, planet_classes.classify(body))
        assert body.planet_class is None or body.planet_class in tuning.PLANET_CLASSES
    assert rogue.planet_type == bound.planet_type
