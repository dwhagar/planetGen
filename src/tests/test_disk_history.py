# tests/test_disk_history.py

"""GEN.469: the thin disk's star formation history comes from gas continuity (the Kennicutt-Schmidt law, an exponential infall
with the time scale growing outward, GEN.1230), and the disk's young, intermediate and old shares are what that history puts in
their age windows.

The Milky Way's figures here (its local gas, its present star formation rate against its mean, the local B-star share)
are test targets that show the derived history is plausible. None of them is an input of the generator."""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import density, disk_history, mass_chain, star_formation
from planetgen.physics.stellar_evolution import sample_living_star, star_params
from planetgen.util import draw


def _shape(length=2600.0, height=300.0, bulge=1580.0, stellar_mass=None):
    masses = None if stellar_mass is None else (stellar_mass, mass_chain.median_halo_mass(stellar_mass))
    return density.build_galaxy_shape(
        disk_scale_length_pc=length, disk_scale_height_pc=height, bulge_scale_radius_pc=bulge, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4, masses=masses)


MILKY_WAY_MASS = 5.0e10
"""float: The mass of the galaxy the tests call Milky Way-like ([R2681]: 6.08 +- 1.14e10, bulge and
disk); its halo is the median one for it, as the mass chain draws."""
SHAPE = _shape(stellar_mass=MILKY_WAY_MASS)
FLAT = {name: (high - low) / (tuning.STAR_FORMATION_AGE_RANGE_GY[1] - tuning.STAR_FORMATION_AGE_RANGE_GY[0])
        for name, (low, high) in ((n, tuning.THIN_DISC_OLD_WINDOW_GY if n == "old" else tuning.STELLAR_POPULATION_AGE_RANGES_GY[n])
                                  for n in tuning.STELLAR_POPULATION_ARM_AMPLITUDE)}


def test_the_return_fraction_reproduces_madau_and_dickinsons_salpeter_value():
    """Their section 2.1: a Salpeter slope of 2.35 between 0.1 and 100 solar masses returns 0.27 of a generation."""
    assert disk_history.mass_return_fraction((0.1, 100.0), (2.35,)) == pytest.approx(0.27, abs=0.01)
    # The generator's Kroupa mass function has more massive stars per unit mass, so more comes back.
    assert 0.27 < disk_history.return_fraction() < 0.5


def test_the_infall_time_scale_is_the_milky_ways_line_between_its_two_models_limits():
    """[R5095] section 3: `tau = 1.033 R - 1.267` Gy, 7 Gy at the Sun; [R5094]: 2 Gy inside and 8 Gy at 20 kpc."""
    sun = 8.0 / tuning.DISK_INFALL_REFERENCE_SCALE_LENGTH_KPC
    assert disk_history.infall_time_gy(sun) == pytest.approx(7.0, abs=0.05)
    assert disk_history.infall_time_gy(0.1) == tuning.DISK_INFALL_MIN_GY
    assert disk_history.infall_time_gy(10.0) == tuning.DISK_INFALL_MAX_GY
    scales = np.linspace(0.0, 10.0, 50)
    assert np.all(np.diff(disk_history.infall_time_gy(scales)) >= 0.0)


def test_the_inflow_of_each_ring_decays_with_that_rings_time_scale():
    rings = np.array([1.0, 3.0])
    shapes = disk_history.inflow_shapes(rings, 10.0, 800, tuning.DISK_INFALL_REFERENCE_SPEED_KMS)
    assert shapes.shape == (801, 2)
    assert np.all(shapes[0] == 1.0)
    # The outer ring's longer time scale keeps more of its inflow to the present (inside-out formation).
    assert shapes[-1, 1] > shapes[-1, 0]
    for column, scale in enumerate(rings):
        assert shapes[-1, column] == pytest.approx(math.exp(-10.0 / disk_history.infall_time_gy(scale)))


def test_a_disc_slower_than_the_milky_ways_is_fed_ever_more_slowly_and_a_faster_one_the_other_way():
    """[R5103] equation 21: the inflow is multiplied by exp(0.4 t (1 - V / 220)) per Gy, so it is the Milky Way's at 220 km/s."""
    rings = np.array([3.0])
    reference = disk_history.inflow_shapes(rings, 10.0, 800, 220.0)
    slow = disk_history.inflow_shapes(rings, 10.0, 800, 100.0)
    fast = disk_history.inflow_shapes(rings, 10.0, 800, 360.0)
    assert slow[-1, 0] / slow[0, 0] > reference[-1, 0] / reference[0, 0] > fast[-1, 0] / fast[0, 0]
    assert slow[-1, 0] / reference[-1, 0] == pytest.approx(math.exp(0.4 * 10.0 * (1.0 - 100.0 / 220.0)))


def test_the_speed_of_a_disc_is_the_cube_root_of_its_mass_from_the_milky_ways_two_hundred_and_twenty():
    assert disk_history.disc_speed_kms(tuning.DISK_INFALL_REFERENCE_MASS_SOLAR) == pytest.approx(220.0)
    assert disk_history.disc_speed_kms(8.0 * tuning.DISK_INFALL_REFERENCE_MASS_SOLAR) == pytest.approx(440.0)


def test_light_discs_rise_to_the_present_and_heavy_ones_fall():
    """[R5104] section 4.1: a disc of 100 km/s has a rising star formation, one of 360 km/s a falling one; the Milky Way's is near constant."""
    light, milky_way, heavy = (disk_history.history_of(_shape(stellar_mass=mass)).present_over_mean() for mass in (5.0e9, 5.0e10, 2.0e11))
    assert light > 1.0 > milky_way > heavy


def test_each_ring_ends_with_the_stars_it_was_solved_to_have():
    found = disk_history.history(3.5e10, 2600.0, disk_history.disc_speed_kms(MILKY_WAY_MASS))
    built = (1.0 - disk_history.return_fraction()) * np.trapezoid(found.rate, dx=10.0 / 800, axis=1)
    assert built == pytest.approx(found.stars_now, rel=1e-4)
    assert np.all(found.gas_now >= 0.0)


def test_the_gas_left_at_the_suns_radius_is_a_fraction_of_its_stars():
    """Target: the Galaxy's local gas is about a third of its local stellar surface density (13 against 38 solar
    masses per pc2); the derived rings should leave tenths, not hundredths or ones."""
    found = disk_history.history(3.5e10, 2600.0, disk_history.disc_speed_kms(MILKY_WAY_MASS))
    ring = int(np.argmin(np.abs(found.radius_scale_lengths - tuning.REFERENCE_RADIUS_SCALE_LENGTHS)))
    assert 0.1 <= found.gas_now[ring] / found.stars_now[ring] <= 0.5


def test_the_shares_add_up_and_differ_from_the_flat_history():
    shares = disk_history.population_shares(SHAPE)
    assert set(shares) == set(tuning.STELLAR_POPULATION_ARM_AMPLITUDE)
    assert sum(shares.values()) == pytest.approx(1.0)
    assert all(share > 0.0 for share in shares.values())
    # A galaxy this massive formed its stars early: fewer recent ones than a flat history gives.
    assert shares["young"] < FLAT["young"]
    assert shares["intermediate"] < FLAT["intermediate"]
    assert shares["old"] > FLAT["old"]


def test_the_share_function_reads_the_same_numbers():
    for name, share in disk_history.population_shares(SHAPE).items():
        assert star_formation.disk_population_share(name, SHAPE) == share


def test_a_shape_without_masses_takes_the_median_galaxy_of_its_disk():
    bare = _shape()
    stellar = disk_history.stellar_mass_for_scale_length(bare.disk_scale_length_pc)
    named = _shape(stellar_mass=stellar)
    assert disk_history.population_shares(bare) == pytest.approx(disk_history.population_shares(named))


def test_lighter_galaxies_have_more_recent_stars():
    """Downsizing: a spiral of 1e10 solar masses is fed more slowly (the infall of [R5103] equation 21) and spends its thinner gas more slowly, so it forms its
    stars later than one of 5e10. The shares fall with mass from 1e10 to 5e10 (the chain's range is 3e9 to 2e11): 2.3 times between them (GEN.1246; 7.5 times
    with the halo-growth inflow, which kept the Milky Way's own disc forming 2.7 times the stars per year it does)."""
    young = [disk_history.population_shares(_shape(stellar_mass=mass))["young"] for mass in (1.0e10, 3.0e10, 5.0e10)]
    assert young == sorted(young, reverse=True)
    assert young[0] > 2.0 * young[-1]


def test_a_flat_inflow_settles_a_ring_toward_a_constant_rate():
    """Gas continuity checked on its own: a constant inflow fills a ring and its rate climbs, ever more slowly, toward
    the rate that uses the inflow up (the second half of ten Gy changes by under 15 percent)."""
    flat = np.ones((801, 1))
    rate, gas = disk_history._solve_rings(np.array([40.0]), flat, 10.0)
    late = rate[0, 400:]
    assert np.all(np.diff(rate[0]) >= 0.0)
    assert late.max() / late.min() < 1.15


def _observed_present_over_mean():
    """Target: [R2681] ([R2680], read): the Galaxy forms 1.65 +- 0.19 solar masses a year (a
    Kroupa mass function) and its disk holds 5.17 +- 1.11e10, which a life of ten Gy and the stars' return fraction
    make 5.17e10 / (1 - R) / 10 Gy of births a year on average. (The disk's mass includes the thick disk, which formed
    before the thin disk's ten Gy, so the ratio is somewhat higher than this.)"""
    mean_births_per_year = 5.17e10 / (1.0 - disk_history.return_fraction()) / 1.0e10
    return 1.65 / mean_births_per_year


def test_the_present_rate_over_the_mean_is_within_a_factor_of_two_and_a_half_of_the_milky_ways():
    """TEST.533, tightened by GEN.1230. The bound was a factor of three while the inflow followed the stars a halo builds, which for the
    lighter median halo of GEN.1154 (1.57e12 solar masses for 5e10 of stars) rose to the present: the ratio was 0.55 against the observed
    0.20. The infall time scale of the Milky Way's chemical evolution model ([R5095], limited by [R5094]) gives 0.45, 2.3 times the observed
    figure, which also counts the thick disc's stars over a ten-Gy life, so it is low for the thin disc alone."""
    found = disk_history.history(MILKY_WAY_MASS * 0.6, 2600.0, disk_history.disc_speed_kms(MILKY_WAY_MASS))
    observed = _observed_present_over_mean()
    assert 0.5 * observed <= found.present_over_mean() <= 2.5 * observed


def _depletion_time_gy(found):
    """The disk's gas over its star formation now, rings weighted by area, Gy (what `radiogenic` reads as the host's depletion time)."""
    radii = found.radius_scale_lengths
    area = radii * (radii[1] - radii[0])
    return float((found.gas_now * area).sum() / ((found.rate[:, -1] * area).sum() * 1.0e-9) / 1.0e9)


def test_a_milky_way_mass_disc_keeps_its_gas_for_the_measured_depletion_time_and_the_local_gas():
    """GEN.999: with the law steepening where the gas thins (GEN.450) the disc of 3e10 solar masses of stars depletes its gas in 2.5
    to 4.2 Gy (1.0 to 4.2 Gy is the measured range, Sources 2 on GEN.998; the plain Kennicutt law gave 1.8) and holds a local gas of the order of
    the measured 11.9 +- 1.6 solar masses per pc2 of atomic and molecular gas with helium ([R3002: Table 2]: 10.9 + 1.0), at least half of it (the plain law gave 6.0 here)."""
    found = disk_history.history(MILKY_WAY_MASS * 0.6, 2600.0, disk_history.disc_speed_kms(MILKY_WAY_MASS))
    assert 2.5 <= _depletion_time_gy(found) <= 4.2
    solar = int(np.argmin(np.abs(found.radius_scale_lengths * 2600.0 - 8000.0)))
    assert 0.5 * 11.9 <= found.gas_now[solar] <= 1.5 * 11.9


def test_the_steepening_leaves_a_thin_ring_more_gas_than_the_plain_law(monkeypatch):
    """The plain law (no steepening: a saturation of zero) lets a thin ring spend its gas as fast as a dense one; the
    steepened law keeps more of it, and a dense ring (whose gas is far over the saturation) hardly changes."""
    flat = np.ones((201, 2))
    stars = np.array([2.0, 400.0])
    steep = disk_history._solve_rings(stars, flat, 10.0)[1]
    monkeypatch.setattr(tuning, "STAR_FORMATION_SATURATION_GAS_MSUN_PC2", 1.0e-12)
    plain = disk_history._solve_rings(stars, flat, 10.0)[1]
    assert steep[0] > 1.3 * plain[0]
    assert steep[1] == pytest.approx(plain[1], rel=0.1)


def _b_share_between_arms_at_the_sun(shape=SHAPE, draws=100_000, spread_ages=True, seed=469):
    """The share of the stars between the arms at the solar radius that are B dwarfs: each disk population's B share
    (drawn from the generator's own sampler, with the ages spread by the history's rate inside each window unless
    `spread_ages` is false, GEN.546) over the populations' densities there."""
    laws = disk_history.population_age_laws(shape) if spread_ages else {}
    solar = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc
    angles = np.linspace(0.0, 2.0 * math.pi, 720, endpoint=False)
    cosines = [density._arm_cosine(solar * math.cos(a), solar * math.sin(a), shape) for a in angles]
    angle = angles[int(np.argmin(cosines))]
    mix = density.population_densities((solar * math.cos(angle), solar * math.sin(angle), 0.0), shape)
    every = sum(max(value, 0.0) for value in mix.values())
    b_share = {}
    with draw.bound(seed):
        for name in ("young", "intermediate"):
            count = 0
            for _ in range(draws):
                mass, age, state = sample_living_star(population=name, age_law=laws.get(name))
                if state["yerkes_class"] in ("VII", "IV", "III", "II", "I"):
                    continue
                count += star_params(mass, age, state)["type"][0] == "B"
            b_share[name] = count / draws
    return sum(mix[name] * b_share[name] for name in b_share) / every


def test_the_local_b_share_is_what_the_derived_history_and_the_mass_function_give():
    """Recorded: the B-star share between the arms at the solar radius is near 0.02 percent for the Milky Way-like
    galaxy with the ages spread by the history's rate (GEN.546; 0.05 before the lifetimes read the birth [Fe/H],
    GEN.651, 0.04 before the main sequence cooled with age, GEN.652: a B dwarf of 2.5 to 3 solar masses is A-type by the
    end of its main sequence), 0.07 percent with them spread evenly inside the windows and 0.19 percent with a flat
    history, against the observed 0.04."""
    share = _b_share_between_arms_at_the_sun()
    assert 0.00012 <= share <= 0.0004
    assert _b_share_between_arms_at_the_sun(spread_ages=False) > share


def test_the_local_b_share_reaches_the_observed_one():
    """Mamajek's census (docs/design/galaxy/star-types-by-galactic-radius.md) is a test target, never an input (GEN.547)."""
    for seed in (469, 1, 2):
        assert _b_share_between_arms_at_the_sun(seed=seed) <= 1.5 * 0.0004
