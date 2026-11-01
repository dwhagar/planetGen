# tests/test_rogue_planet_parity.py

"""TEST.498 (Boss, 4:57 PM EDT 2026-10-30): a rogue planet is functionally and programmatically identical to a bound planet, differing only in
the rogue note, the missing host and the quantities that depend on starlight (its equilibrium temperature is the internal and ambient value).

The rogue planet is a separate generator today (`phenomena/rogue.RoguePlanet`: its own mass, radius, surface and class code, 72 columns under
its own names, `has_moons` a flag and no moon system). GEN.1111 (rogue rows from the architecture's ejections) and the rogue-train items that
follow it make it a `Planet` without a star. This file is the contract they are built to:

* the inventory tests pass now and measure the gap (`MISSING_CEILING` can only fall: a lower number is written into it by whoever closes
  columns);
* the contract tests are strict expected failures until the rogue builds through the planet code; when one starts to pass, its mark comes off
  (the strict mark makes the run fail until it does).
"""

import pytest

from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.planet import Planet

HOST_ONLY = frozenset({
    # where the planet is around its host: no host, none of it
    "zone", "habitable_zone", "distance", "hill_radius", "min_orbit_distance", "period", "orbital_inclination_deg",
    "orbital_ascending_node_deg", "eccentricity", "position_x", "position_y", "position_z", "orbital_speed_kms", "velocity_x_kms",
    "velocity_y_kms", "velocity_z_kms", "min_update_interval_years", "reflex_offset_x", "reflex_offset_y", "reflex_offset_z",
    "host_eccentricity", "born_beyond_snow_line", "parent_shadow_share", "solar_day_hours", "is_moon",
})
"""frozenset: Planet columns a rogue planet need not carry (the host and the orbit round it). Everything else a bound planet stores a rogue
stores too; a starlight quantity (the equilibrium temperature, the sky, the UV) is a column of both, with the internal and ambient value for the
rogue."""

ROGUE_ONLY_PREFIXES = ("rogue_", "galactic_", "ejection_", "parent_system")
"""tuple: The columns only a rogue planet has: its note (`rogue_note`), its motion round the galaxy in place of an orbit, and the system it was
ejected from (GEN.1111)."""

MISSING_CEILING = 202
"""int: Planet columns (net of `HOST_ONLY`) the rogue planet does not carry under the same name: 148 when the test was written (TEST.498), plus the
columns the planet gained while the rogue is still a generator of its own (below). It falls to zero when the rogue is built through the planet code
(GEN.1111 and the rogue train, where every one of them comes with the planet); until then it rises only by the planet columns a lane adds, named in
`PLANET_COLUMNS_ADDED`, never by a column the rogue lost.

TEST.512 (found red on main at 8:48 AM EDT 2026-10-31): the gap was 160, not 148, and no change to the rogue path had opened it. Twelve planet
columns were added by the planet items listed in `PLANET_COLUMNS_ADDED`; the rogue's own columns are the same as when the test was written."""
ORIGINAL_GAP = 134
"""int: The 148 columns of TEST.498, less the 14 old latitude temperature columns (`temp_polar_k` to `temp_terminator_swing_k`) that GEN.1108 part B dropped from
the planet: the gap of the columns that were there when the test was written and still are."""
PLANET_COLUMNS_ADDED = {
    "GEN.474 (intermediate ground)": ("ground_intermediate_fraction",),
    "GEN.768 (hours outside)": ("outside_hours_day", "outside_hours_limit"),
    "GEN.1034 (open-ocean tide)": ("ocean_liquid", "tide_range_neap_m", "tide_range_spring_m"),
    "GEN.1035 (stored surface liquid)": ("surface_liquid", "surface_liquid_depth_km", "surface_liquid_fraction"),
    "GEN.1215 (what the surface liquids hold)": ("surface_liquid_acid_wt", "surface_liquid_acid_layer", "surface_liquid_ammonia_wt"),
    "GEN.1159 (roundness, strength and binding energy)": ("roundness_ratio", "yield_strength_pa", "binding_energy_j"),
    "GEN.1097 (day and night state)": ("day_night_state", "days_per_year", "night_hours"),
    "GEN.983 (belt dose of a moon)": ("dose_belt_msv_yr",),
    "GEN.1108 (seven latitude bands)": tuple(f"band{n}_{name}" for n in range(1, 8) for name in (
        "mean_k", "swing_k", "day_hours", "night_hours", "longest_day_days", "longest_night_days", "state")),
}
"""dict: The planet columns added since TEST.498 was written, by the item that added them. A lane that adds a planet column raises
`MISSING_CEILING` by its count and lists the names here; `test_the_columns_added_since_the_test_was_written_are_all_missing_from_the_rogue` checks
them against the gap, so a name that is not a planet column, or that the rogue now carries, is caught."""
EXTRA_CEILING = 13
"""int: Rogue columns that are neither a planet column nor a rogue-only one (the rogue's own names for what a planet stores, `mass_kg` for
`mass`). It falls to zero with the above."""


def _planet_columns():
    return set(Planet.SERIALIZABLE_FIELDS)


def _rogue_columns():
    return set(RoguePlanet.SERIALIZABLE_FIELDS)


def missing_columns():
    """The columns of a bound planet that a rogue planet lacks, apart from those that need a host."""
    return sorted(_planet_columns() - HOST_ONLY - _rogue_columns())


def extra_columns():
    """The rogue columns that are not a planet column and not rogue-only."""
    return sorted(name for name in _rogue_columns() - _planet_columns() if not name.startswith(ROGUE_ONLY_PREFIXES))


# --- the inventory (passes now; the gap can only shrink) -----------------------------------------------------------------------------

def test_the_host_only_columns_are_all_planet_columns():
    assert HOST_ONLY <= _planet_columns(), sorted(HOST_ONLY - _planet_columns())


def test_the_gap_between_a_rogue_and_a_bound_planet_never_grows():
    missing, extra = missing_columns(), extra_columns()
    assert len(missing) <= MISSING_CEILING, f"{len(missing)} planet columns missing from a rogue: {missing}"
    assert len(extra) <= EXTRA_CEILING, f"{len(extra)} rogue-only names for planet columns: {extra}"


def test_the_columns_added_since_the_test_was_written_are_all_missing_from_the_rogue():
    added = {name for names in PLANET_COLUMNS_ADDED.values() for name in names}
    assert len(added) == MISSING_CEILING - ORIGINAL_GAP, "the ceiling rises by exactly the columns named"
    assert added <= set(missing_columns()), sorted(added - set(missing_columns()))


def test_the_gap_ceilings_are_the_gap_until_the_rogue_is_built_through_the_planet_code():
    """A ceiling a long way above the gap hides a regression; the two numbers are the measured gap, to within 10 and 3."""
    assert MISSING_CEILING - len(missing_columns()) < 10 and EXTRA_CEILING - len(extra_columns()) < 3


# --- the contract (strict expected failures until GEN.1111 and the rogue train build it) ---------------------------------------------

_UNTIL = "the rogue is a RoguePlanet generator of its own; GEN.1111 and the rogue train (TEST.498) build it through the planet code"


@pytest.mark.xfail(strict=True, reason=_UNTIL)
def test_a_rogue_planet_carries_every_column_a_bound_planet_does():
    assert not missing_columns(), missing_columns()


@pytest.mark.xfail(strict=True, reason=_UNTIL)
def test_a_rogue_planet_has_no_columns_of_its_own_but_the_note_the_motion_and_the_ejecting_system():
    assert not extra_columns(), extra_columns()


@pytest.mark.xfail(strict=True, reason=_UNTIL)
def test_a_rogue_planet_is_a_planet_so_it_runs_the_planet_code_and_no_rogue_only_generator():
    assert issubclass(RoguePlanet, Planet)


@pytest.mark.xfail(strict=True, reason=_UNTIL)
def test_a_rogue_planet_holds_a_moon_system_of_planets_like_a_bound_one():
    """Moons are `Planet` objects with `is_moon` (`planetPhysics.generate_moons`); a rogue today stores a flag."""
    assert hasattr(RoguePlanet, "add_moons") and "has_moons" not in _rogue_columns()


@pytest.mark.xfail(strict=True, reason=_UNTIL)
def test_a_rogue_planet_takes_its_class_from_the_same_classifier_as_a_bound_planet():
    """`planetPhysics.apply_class` is the one place a body's class is set; the rogue has its own `rogue_class` wrapper."""
    from planetgen.generation.phenomena import rogue
    assert not hasattr(rogue, "rogue_class")


# --- the matched draw (TEST.498: draw matched planets with and without a host and compare what the two chains make) ---------------------
# Until a starless host exists (GEN.1111) the bound side is a `Planet` built `FAR_AU` from a G2 star with the rogue's own mass: the host's light
# is then the ambient floor, which is the only thing a rogue has in its place.

FAR_AU = 1.0e5
BINS = ("terrestrial", "sub-neptune", "saturn", "jupiter")
"""tuple: The rogue mass bins matched (the brown dwarf is a star's track, not a planet's)."""


def _matched(mass_bin, seed):
    """`(rogue, bound)`: a `RoguePlanet` of the bin and a `Planet` of the same mass far from a star, or the test is skipped where a dataset the
    giant needs is not on this machine."""
    from planetgen.datasets.fetch import DataMissing
    from planetgen.util import draw
    from tests.test_gen_builder_internals import make_system
    system = make_system()
    star = system.star
    try:
        with draw.bound(f"rogue-parity-{mass_bin}-{seed}"):
            rogue = RoguePlanet(system.system_config, mass_bin=mass_bin)
        with draw.bound(f"rogue-parity-bound-{mass_bin}-{seed}"):
            bound = Planet(system.system_config, star, star.habitable_zone, FAR_AU, mass=rogue.mass_kg, moon_count=0)
    except DataMissing as missing:
        pytest.skip(f"a dataset the {mass_bin} needs is not on this machine: {missing}")
    return rogue, bound


@pytest.mark.parametrize("mass_bin", BINS)
@pytest.mark.parametrize("seed", range(3))
def test_a_rogue_and_a_bound_planet_of_one_mass_are_one_kind_of_body_of_one_size(mass_bin, seed):
    from planetgen.physics import constants
    rogue, bound = _matched(mass_bin, seed)
    assert rogue.mass_kg == bound.mass
    assert rogue.planet_type == bound.body_type
    if rogue.planet_type == "t":
        # a rocky body's radius follows its density, drawn from the same range in both chains
        low, high = (value * 1000.0 for value in constants.PLANET_DENSITY["t"])
        radii = [(rogue.mass_kg / (4.0 / 3.0 * 3.141592653589793 * density)) ** (1.0 / 3.0) / 1000.0 for density in (high, low)]
        assert radii[0] * 0.8 <= bound.radius <= radii[1] * 1.2, (bound.radius, radii)


@pytest.mark.xfail(strict=True, reason="the bound chain's equilibrium temperature has no internal heat in it (GEN.1111 slice: the starless host)")
@pytest.mark.parametrize("mass_bin", BINS)
def test_a_rogue_and_a_bound_planet_far_from_light_have_the_same_surface_temperature(mass_bin):
    """The one quantity the rogue and the bound chain should differ in is starlight; with none, the surface is the internal heat's over the ambient
    floor in both. Found (TEST.498, 2026-10-31): the bound chain returns the 2.725 K floor for a 5 Earth-mass body at 1e5 AU and the rogue chain 270 K."""
    rogue, bound = _matched(mass_bin, 0)
    assert abs(bound.surface_temperature - rogue.surface_temperature_k) <= 0.1 * rogue.surface_temperature_k, (
        bound.surface_temperature, rogue.surface_temperature_k)
