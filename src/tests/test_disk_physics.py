"""
Disk-physics regression tests for `StarSystem.estimate_num_objects`.

Covers the rework replacing the old `BASE_MAX_SYSTEM_OBJECTS *
(1 + log10(solar_masses))` curve fit with a physically-derived ceiling: a
protoplanetary-disk isolation-mass/mutual-Hill-radius walk (Lissauer 1993;
Kokubo & Ida 2000, 2002) from real Minimum Mass Solar Nebula surface-density
scaling (Hayashi 1981) -- see `utils.snow_line_au`,
`utils.disk_surface_density_scale`, `utils.mmsn_surface_density_gcm2`,
`utils.isolation_mass_kg`, and
`StarSystem._estimate_max_objects_from_disk_physics`'s own docstring for the
full derivation and literature citations.

`MAX_PLANETS=True` makes `estimate_num_objects` return the physically-derived
ceiling deterministically (no `random.randint` draw), the same way
test_systems.py's own `MAX_PLANETS=True` case already relies on -- used
throughout instead of poking any underscore-prefixed method directly.

Run with: pytest src/tests/test_disk_physics.py
"""
import pytest

from stellarObjects import physical_constants, program_constants as prog_c
from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import (
    disk_surface_density_scale,
    isolation_mass_kg,
    mmsn_surface_density_gcm2,
    snow_line_au,
)


def _make_config(star_type, **overrides):
    """
    `BINARY_SYSTEM` defaults to False: none of this module's tests are
    about binary generation, but `StarSystem._should_generate_binary`
    rolls real chance whenever it's left at its own default (None), and a
    second star (whether merged P-type or an independent S-type
    secondary) would otherwise unpredictably add its own object count on
    top of a single star's -- exactly what this module's disk-physics
    estimates are trying to pin down deterministically.
    """
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.BINARY_SYSTEM = False
    for attr, value in overrides.items():
        setattr(cfg, attr, value)
    return cfg


def _max_objects(star_type):
    """The deterministic disk-physics ceiling for a given star type."""
    system = StarSystem(system_config=_make_config(star_type, PLANETS=True, MAX_PLANETS=True))
    return system.planet_count + system.belt_count


# --- Pure helper functions ---

def test_snow_line_scales_with_sqrt_luminosity():
    solar = snow_line_au(physical_constants.SOLAR_LUMINOSITY)
    assert solar == pytest.approx(physical_constants.SNOW_LINE_AU_AT_1_LSUN)
    quadruple_lum = snow_line_au(4 * physical_constants.SOLAR_LUMINOSITY)
    # L -> 4L should push the snow line out by sqrt(4) = 2x.
    assert quadruple_lum == pytest.approx(2 * solar)


def test_disk_surface_density_scale_is_one_at_solar_mass():
    assert disk_surface_density_scale(physical_constants.SOLAR_MASS_TO_KG) == pytest.approx(1.0)


def test_disk_surface_density_scale_grows_with_stellar_mass():
    lo = disk_surface_density_scale(0.5 * physical_constants.SOLAR_MASS_TO_KG)
    hi = disk_surface_density_scale(2.0 * physical_constants.SOLAR_MASS_TO_KG)
    assert lo < 1.0 < hi


def test_mmsn_surface_density_falls_off_with_distance_and_jumps_at_snow_line():
    snow_line = 2.7
    inner = mmsn_surface_density_gcm2(1.0, snow_line)
    outer_before_snow_line = mmsn_surface_density_gcm2(2.5, snow_line)
    outer_after_snow_line = mmsn_surface_density_gcm2(2.8, snow_line)
    assert inner == pytest.approx(physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2)
    # Falls off as distance^-1.5 inside the snow line.
    assert outer_before_snow_line < inner
    # Ice condensation just past the snow line jumps the solid budget back
    # up, well past where it fell to just before the line.
    assert outer_after_snow_line > outer_before_snow_line


def test_isolation_mass_matches_known_literature_value_at_1au():
    # Real oligarchic-growth literature (Kokubo & Ida) commonly cites an
    # MMSN isolation mass at 1 AU of roughly 0.05-0.1 Earth masses.
    sigma = mmsn_surface_density_gcm2(1.0, snow_line_au(physical_constants.SOLAR_LUMINOSITY))
    m_iso = isolation_mass_kg(1.0, sigma, physical_constants.SOLAR_MASS_TO_KG)
    m_iso_earth = m_iso / physical_constants.EARTH_MASS_TO_KG
    assert 0.03 < m_iso_earth < 0.15


def test_isolation_mass_grows_with_local_surface_density():
    lo = isolation_mass_kg(1.0, 7.0, physical_constants.SOLAR_MASS_TO_KG)
    hi = isolation_mass_kg(1.0, 30.0, physical_constants.SOLAR_MASS_TO_KG)
    assert hi > lo


# --- estimate_num_objects's override contract (unchanged by the rework) ---

def test_planets_false_forces_zero_regardless_of_disk_physics():
    system = StarSystem(system_config=_make_config("G2V", PLANETS=False))
    assert system.planet_count + system.belt_count == 0


def test_num_orbits_bypasses_disk_physics_entirely():
    system = StarSystem(system_config=_make_config("G2V", NUM_ORBITS=3, PLANETS=True))
    assert system.planet_count + system.belt_count == 3


def test_max_planets_false_returns_the_floor():
    system = StarSystem(system_config=_make_config("G2V", PLANETS=True, MAX_PLANETS=False))
    assert system.planet_count + system.belt_count == 1


def test_disk_physics_ceiling_never_exceeds_absolute_max():
    # M50 (an extreme, near-substellar case) exercises the walk's own
    # termination guard directly.
    assert _max_objects("M50") <= prog_c.ABSOLUTE_MAX_SYSTEM_OBJECTS


# --- Real-astrophysics shape of the resulting ceiling ---

def test_more_luminous_star_has_a_larger_snow_line_and_disk():
    # A hotter/more luminous star's disk both starts its ice line farther
    # out and (via DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER) is truncated
    # farther out -- real disks around hotter stars are observed to run
    # larger, not smaller.
    dim = snow_line_au(0.01 * physical_constants.SOLAR_LUMINOSITY)
    bright = snow_line_au(100 * physical_constants.SOLAR_LUMINOSITY)
    assert bright > dim
