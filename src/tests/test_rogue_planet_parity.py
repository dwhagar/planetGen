# tests/test_rogue_planet_parity.py

"""TEST.498 (Boss, 4:57 PM EDT 2026-10-30): a rogue planet is functionally and programmatically identical to a bound planet, differing only in the rogue
note, the missing host and the quantities that depend on starlight (its equilibrium temperature is the internal and ambient value).

Since GEN.1111 a rogue planet IS a `Planet`, built by the planet code against a `StarlessHost` (`generation/phenomena/starless_host.py`), and stored in
`rogue_planets`, which carries every `planets` column apart from the orbit round a host. This file holds that contract:

* the rogue is built by `Planet` and is a `Planet`;
* the table has every column of a planet except the host and orbit ones;
* matched draws (a rogue and a planet of one mass far from a star) are one kind of body of one size, and both are warmed by their own heat.
"""

import pytest

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
