# tests/test_planet_architecture.py

"""
GEN.728, GEN.729, GEN.730, GEN.828, GEN.829, GEN.920: the planet architecture (`planetgen.physics.planet_architecture`): Pu and Wu's
survival law and the thinning it drives, the inner system and the outer chain, the gas and the migration of giants, and the sky of a
few thousand seeded stars against the observed figures of `planet_occurrence` (the falsifiers of the registry rows of section 6.36 of the
gap filling plan).
"""

import math

import pytest

from planetgen import tuning
from planetgen.physics import constants, formation, kepler_multiplicity, planet_architecture as pa, planet_mass, planet_occurrence as po
from planetgen.util import draw

SUN = constants.SOLAR_MASS_TO_KG
EARTH = constants.EARTH_MASS_TO_KG


def sun_like():
    return {"star_mass_kg": SUN, "snow_line_au": 2.7, "metal_ratio": 1.0, "fe_h": 0.0, "outer_au": 48.0}


def period_d(star_mass_kg, a_au):
    return pa.orbital_period_s(a_au, star_mass_kg) / 86400.0


# --- Pu and Wu's survival law ----------------------------------------------------------------------------------------------


def test_the_survival_law_gives_half_at_the_published_critical_spacing():
    for log_tau in (4.0, 7.0, 10.5):
        k50 = 0.7 * log_tau + 2.87
        assert pa.survival_probability(k50, log_tau) == pytest.approx(0.5, abs=0.01)


def test_a_wider_pair_survives_better_and_a_longer_time_worse():
    assert pa.survival_probability(12.0, 9.0) > pa.survival_probability(10.0, 9.0) > pa.survival_probability(8.0, 9.0)
    assert pa.survival_probability(10.0, 7.0) > pa.survival_probability(10.0, 9.0) > pa.survival_probability(10.0, 11.0)
    assert pa.survival_probability(1.0, 0.0) == 0.0                     # under the law's own minimum spacing
    assert pa.survival_probability(30.0, 9.0) > 0.99


def test_the_drawn_instability_time_follows_the_survival_law():
    draw.set_run_seed(4)
    k, log_tau, n = 10.0, 10.0, 20000
    lasted = sum(1 for _ in range(n) if pa.instability_log_tau(k, 1.0 - draw.random()) > log_tau)
    assert lasted / n == pytest.approx(pa.survival_probability(k, log_tau), abs=0.012)


def test_a_pair_closer_than_the_law_allows_goes_at_once_and_none_outlasts_the_cap():
    assert pa.instability_log_tau(1.0, 0.5) == 0.0
    assert pa.instability_log_tau(500.0, 0.5) == pa.MAX_LOG_TAU


# --- spacing and the thinning ----------------------------------------------------------------------------------------------


def test_the_spacing_of_two_planets_is_their_gap_in_mutual_hill_radii():
    m = 3.0 * EARTH
    r_h = pa.mutual_hill_radius_au(1.0, m, 1.2, m, SUN)
    assert r_h == pytest.approx(1.1 * (2.0 * m / (3.0 * SUN)) ** (1.0 / 3.0))
    assert pa.spacing_k(1.0, m, 1.2, m, SUN) == pytest.approx(0.2 / r_h)


def test_a_well_spaced_chain_is_left_alone_and_a_tight_one_is_thinned():
    draw.set_run_seed(8)
    wide = [[0.1 * 1.8 ** i, 1.0 * EARTH] for i in range(6)]
    assert [b[:2] for b in pa.sculpt(wide, SUN, 1e9)] == wide
    tight = [[0.1 * 1.04 ** i, 1.0 * EARTH] for i in range(12)]                 # touching: they all fall together
    thinned = pa.sculpt(tight, SUN, 1e9)
    assert 1 <= len(thinned) < len(tight)
    assert sum(b[2] for b in thinned) == len(tight)                      # every oligarch is in some planet
    distances = [b[0] for b in thinned]
    assert distances == sorted(distances)


def test_a_collision_keeps_the_angular_momentum_and_the_retained_share_of_the_mass():
    first, second = [1.0, 2.0 * EARTH, 1], [1.4, 1.0 * EARTH, 1]
    merged = pa._merge(first, second)
    assert merged[1] == pytest.approx(3.0 * EARTH * tuning.COLLISION_MASS_RETENTION) and merged[2] == 2
    assert math.sqrt(merged[0]) == pytest.approx((2.0 * math.sqrt(1.0) + 1.0 * math.sqrt(1.4)) / 3.0)


def test_what_is_left_after_the_thinning_is_spaced_wider_than_the_oligarchs_were():
    draw.set_run_seed(15)
    chain = pa.oligarch_chain(SUN, 2.7, 1.0, 1.0, 6.0, 0.05, 1.0)
    before = [pa.spacing_k(a[0], a[1], b[0], b[1], SUN) for a, b in zip(chain, chain[1:])]
    after_bodies = pa.sculpt(chain, SUN, 1e9)
    after = [pa.spacing_k(a[0], a[1], b[0], b[1], SUN) for a, b in zip(after_bodies, after_bodies[1:])]
    assert sorted(after)[len(after) // 2] > sorted(before)[len(before) // 2]
    assert sum(b[1] for b in after_bodies) <= sum(b[1] for b in chain) + 1.0


# --- the inner system ------------------------------------------------------------------------------------------------------


def test_the_share_of_stars_with_an_inner_system_follows_the_observed_anchors():
    shares = {m: pa.inner_system_fraction(m * SUN) for m in (0.08, 0.15, 0.6, 1.0, 1.4, 3.0)}
    assert shares[0.6] == pytest.approx(0.75) and shares[1.0] == pytest.approx(0.30)
    assert shares[0.08] == shares[0.15] and shares[1.4] == shares[3.0]
    assert shares[0.15] > shares[0.6] > shares[1.0]
    assert 0.30 < pa.inner_system_fraction(0.8 * SUN) < 0.75
    assert po.compare(pa.inner_system_fraction(SUN), po.BANDS["inner_fraction.G"]).ok
    assert po.compare(pa.inner_system_fraction(0.5 * SUN), po.BANDS["inner_fraction.M"]).ok


def test_an_inner_system_starts_at_the_inner_edge_and_neighbours_are_alike():
    draw.set_run_seed(2)
    found = 0
    for _ in range(60):
        chain = pa.inner_chain(0.5 * SUN, 1.0, 0.0, None)
        if len(chain) < 2:
            continue
        found += 1
        assert period_d(0.5 * SUN, chain[0][0]) < 30.0
        masses = [m for _a, m in chain]
        assert max(masses) / min(masses) < 10.0 ** (6 * tuning.INNER_NEIGHBOUR_SIGMA_DEX + len(chain) * tuning.INNER_OUTWARD_MASS_DEX_PER_PLANET)
        assert [a for a, _m in chain] == sorted(a for a, _m in chain)
    assert found > 10


def test_a_heavier_disk_builds_heavier_inner_planets():
    def median_mass(scatter):
        draw.set_run_seed(6)
        masses = [m for _ in range(200) for _a, m in pa.inner_chain(SUN, scatter, 0.0, None)]
        return sorted(masses)[len(masses) // 2]
    assert median_mass(5.0) > median_mass(1.0) > median_mass(0.2)


def test_systems_differ_in_scale_and_the_planets_of_one_system_stay_alike():
    """GEN.1335: with no spread between systems every Sun-like inner system started at the same 3.6 Earth masses and the Earth-size planets the surveys count were missing.
    The mean of a system's log mass now spreads by about `INNER_SYSTEM_SCALE_SIGMA_DEX`, and the planets of one system stay within the neighbour spread."""
    draw.set_run_seed(31)
    means, within = [], []
    for _ in range(400):
        chain = pa.inner_chain(SUN, 1.0, 0.0, None, force=True)
        logs = [math.log10(m / EARTH) for _a, m in chain]
        if len(chain) >= 3:
            means.append(sum(logs) / len(logs))
            within.append(max(logs) - min(logs))
    spread = (sum((x - sum(means) / len(means)) ** 2 for x in means) / len(means)) ** 0.5
    assert len(means) > 200
    assert tuning.INNER_SYSTEM_SCALE_SIGMA_DEX * 0.8 < spread < math.hypot(tuning.INNER_SYSTEM_SCALE_SIGMA_DEX, tuning.INNER_NEIGHBOUR_SIGMA_DEX) * 1.3, spread
    assert sorted(within)[len(within) // 2] < 6 * tuning.INNER_NEIGHBOUR_SIGMA_DEX
    assert sum(means) / len(means) == pytest.approx(math.log10(tuning.INNER_PLANET_MASS_EARTH), abs=0.08)       # the median of the scale is unchanged


def test_a_share_of_the_inner_systems_of_a_k_dwarf_hold_an_earth_size_planet():
    """GEN.1335 against [R1459] (16.5 +- 3.6 percent of FGK stars have a planet of 0.8 to 1.25 Earth radii inside 85 days): of the inner systems of a K dwarf
    about a quarter hold a bare core of that size (0.19 with every system at the same scale). GEN.1406 lowered the median scale to 2.5 Earth masses (the rocky peak of the radius
    histogram at 1.3 Earth radii, [R1478: Abstract]): 0.38 of the inner systems of a K dwarf and 0.33 of a Sun's."""
    mass_solar, hosts, earth = 0.7, 0, 0
    draw.set_run_seed(78)
    for _ in range(600):
        z = 0.2 * draw.gauss(0.0, 1.0)
        arch = pa.build(mass_solar * SUN, 1.2, 10.0 ** z, z, 30.0, inner=True)
        short = [x for x in arch.planets if not x.runaway and period_d(mass_solar * SUN, x.distance_au) < 85.0]
        hosts += bool(short)
        earth += any(0.8 <= _radius_earth(x.mass_kg / EARTH) <= 1.25 for x in short)
    assert hosts > 500
    assert 0.30 <= earth / hosts <= 0.46, earth / hosts       # 0.19 with no spread between systems, 0.25 at a median scale of 3.6 Earth masses


def test_a_share_of_the_inner_systems_of_a_sun_like_star_hold_an_earth_size_planet():
    """GEN.1406: with the Sun's inner systems at a median of 3.6 Earth masses only 0.21 of them held a bare core of 0.8 to 1.25 Earth radii inside 85 days, which left 5.5 percent of
    G2V stars one against the 16.5 +- 3.6 percent [R1459] finds for FGK stars; at 2.5 Earth masses the share is a third of the hosts."""
    hosts, earth = 0, 0
    draw.set_run_seed(79)
    for _ in range(800):
        z = 0.2 * draw.gauss(0.0, 1.0)
        arch = pa.build(SUN, 1.2, 10.0 ** z, z, 30.0, inner=True)
        short = [x for x in arch.planets if not x.runaway and period_d(SUN, x.distance_au) < 85.0]
        hosts += bool(short)
        earth += any(0.8 <= _radius_earth(x.mass_kg / EARTH) <= 1.25 for x in short)
    assert hosts > 700
    assert 0.27 <= earth / hosts <= 0.42, earth / hosts


# --- giants ----------------------------------------------------------------------------------------------------------------


def test_the_chance_of_a_hot_jupiter_rises_steeply_with_metallicity():
    low, solar, high = (pa.migration_probability(x) for x in (-0.4, 0.0, 0.4))
    assert low < solar < high
    assert solar == pytest.approx(tuning.GIANT_MIGRATION_SHARE)
    assert high / low == pytest.approx(10.0 ** (0.8 * tuning.GIANT_MIGRATION_METALLICITY_EXTRA_EXPONENT))
    assert pa.migration_probability(5.0) == 1.0


def test_a_hot_jupiter_has_nothing_inside_its_orbit_and_a_moved_giant_clears_its_neighbours():
    hot = pa.Seed(0.05, 300 * EARTH, 30 * EARTH, True, 4, True, 0.05)
    inside = pa.Seed(0.03, 2 * EARTH, 2 * EARTH, False, 1, False, 0.03)
    outside = pa.Seed(1.5, 2 * EARTH, 2 * EARTH, False, 1, False, 1.5)
    assert pa._clear_around_giants([inside, hot, outside], SUN) == [hot, outside]
    moved = pa.Seed(2.0, 300 * EARTH, 30 * EARTH, True, 4, True, 2.0)
    near = pa.Seed(2.2, 5 * EARTH, 5 * EARTH, False, 1, False, 2.2)
    far = pa.Seed(9.0, 5 * EARTH, 5 * EARTH, False, 1, False, 9.0)
    assert pa._clear_around_giants([moved, near, far], SUN) == [moved, far]
    still = pa.Seed(2.0, 300 * EARTH, 30 * EARTH, True, 4, False, 2.0)
    assert pa._clear_around_giants([still, far], SUN) == [still, far]      # a giant that did not move clears nothing
    assert pa._clear_around_giants([still, near, far], SUN) == [still, far]      # but one that landed on a neighbour throws it out


def test_two_giants_that_land_on_one_another_leave_the_heavier():
    heavy = pa.Seed(1.57, 300 * EARTH, 30 * EARTH, True, 4, False, 4.0)
    light = pa.Seed(1.59, 80 * EARTH, 20 * EARTH, True, 3, False, 6.0)
    apart = pa.Seed(8.0, 80 * EARTH, 20 * EARTH, True, 3, False, 8.0)
    assert pa._clear_around_giants([heavy, light, apart], SUN) == [heavy, apart]
    assert pa._clear_around_giants([light, heavy, apart], SUN) == [heavy, apart]


# --- the whole architecture ------------------------------------------------------------------------------------------------


def test_a_star_gets_its_planets_in_order_with_the_masses_a_planet_has_and_a_seeded_run_repeats():
    kwargs = sun_like()
    draw.set_run_seed(31)
    first = [pa.build(**kwargs) for _ in range(40)]
    draw.set_run_seed(31)
    second = [pa.build(**kwargs) for _ in range(40)]
    assert first == second
    for architecture in first:
        distances = [s.distance_au for s in architecture.planets]
        assert distances == sorted(distances)
        for s in architecture.planets:
            assert s.merged >= 1 and s.mass_kg >= 0.9 * tuning.PLANET_MIN_MASS_EARTH * EARTH      # the soft cap shaves a little off a core
            assert s.core_kg <= s.mass_kg or s.runaway


def test_the_architecture_of_a_dim_disk_is_empty_and_of_a_massive_one_is_not():
    draw.set_run_seed(13)
    counts = [len(pa.build(**sun_like()).planets) for _ in range(300)]
    assert min(counts) == 0 and max(counts) >= 6


def _sky(mass_solar, fe_h, n, seed, spread=0.0):
    """`n` stars of `mass_solar`; each has [Fe/H] `fe_h` plus a normal of width `spread` (the survey stars are a mix, and the hot Jupiters
    go as a steep power of the metallicity, so their rate is that of the mix, not of its middle)."""
    draw.set_run_seed(seed)
    metallicities = [fe_h + spread * draw.gauss(0.0, 1.0) for _ in range(n)]
    mass_kg = mass_solar * SUN
    luminosity_w = 3.828e26 * (mass_solar ** 4.0 if mass_solar > 0.43 else 0.23 * mass_solar ** 2.3)
    snow = formation.snow_line_au(luminosity_w)
    return mass_kg, snow, [pa.build(mass_kg, snow, 10.0 ** z, z, min(18.0 * snow, 100.0)) for z in metallicities]


def _radius_earth(mass_earth):
    """[R0758] mean mass-radius law (the three branches), Earth radii: the sky's size without the envelope loss."""
    if mass_earth < 2.04:
        return 1.008 * mass_earth ** 0.279
    if mass_earth < 131.6:
        return 0.808 * mass_earth ** 0.589
    return 17.78 * mass_earth ** -0.044


def _inner(architecture, mass_kg):
    return [s for s in architecture.planets if pa.orbital_period_s(s.distance_au, mass_kg) / 86400.0 < 400.0
            and 1.0 <= _radius_earth(s.mass_kg / EARTH) <= 20.0]


def test_a_sun_like_sky_has_the_inner_systems_the_surveys_count():
    mass_kg, _snow, sky = _sky(1.0, 0.0, 1500, 101)
    hosts = [_inner(a, mass_kg) for a in sky]
    share = sum(1 for h in hosts if h) / len(sky)
    assert po.compare(share, po.BANDS["inner_fraction.G"], po.share_sigma(share, len(sky))).ok, share
    mean = sum(len(h) for h in hosts if h) / max(1, sum(1 for h in hosts if h))
    band = po.BANDS["inner_multiplicity.G"]
    assert po.compare(mean, band, 0.1).ok, mean


def test_the_spacing_of_neighbours_peaks_near_twenty_hill_radii_and_few_are_closer_than_ten():
    mass_kg, _snow, sky = _sky(1.0, 0.0, 800, 202)
    spacings = []
    for architecture in sky:
        p = architecture.planets
        spacings += [pa.spacing_k(a.distance_au, a.mass_kg, b.distance_au, b.mass_kg, mass_kg) for a, b in zip(p, p[1:])]
    assert len(spacings) > 800
    wide = sum(1 for k in spacings if k >= 10.0) / len(spacings)
    assert po.compare(wide, po.BANDS["spacing_k_ge_10"], po.share_sigma(wide, len(spacings))).ok, wide
    median = sorted(spacings)[len(spacings) // 2]
    assert 12.0 < median < 28.0, median


def test_giants_rise_with_the_mass_of_the_star_and_with_its_metallicity():
    def giants(mass_solar, fe_h, seed):
        _mass, _snow, sky = _sky(mass_solar, fe_h, 1500, seed)
        return sum(1 for a in sky if any(s.runaway for s in a.planets)) / len(sky)

    assert giants(0.3, 0.0, 1) < giants(1.0, 0.0, 2)
    assert giants(1.0, -0.4, 3) * 2.0 < giants(1.0, 0.4, 3)


def test_an_m_dwarf_has_giants_within_2_5_au_as_often_as_the_survey_finds():
    # GEN.1336 / [R2172]: 3.5 percent of 0.5 solar mass stars at solar metallicity, in proportion to the star's mass; the model's snow line lies close in round
    # an M dwarf, so before the gas of a lighter star's disk fell with its mass the rate was about 7 percent.
    _mass, _snow, sky = _sky(0.5, 0.0, 3000, 808)
    share = sum(1 for a in sky if any(s.runaway and s.distance_au < 2.5 for s in a.planets)) / len(sky)
    assert 0.02 < share < 0.055, share
    _mass, _snow, sky = _sky(1.0, 0.0, 3000, 809)
    sun = sum(1 for a in sky if any(s.runaway and s.distance_au < 2.5 for s in a.planets)) / len(sky)
    assert share < sun, (share, sun)


def test_cold_giants_of_a_sun_like_star_are_as_common_as_the_survey_finds_them():
    mass_kg, _snow, sky = _sky(1.0, 0.0, 2500, 303)
    share = sum(1 for a in sky if any(s.runaway and 2.0 <= s.distance_au <= 8.0 for s in a.planets)) / len(sky)
    assert po.compare(share, po.BANDS["giants_2_8au"], po.share_sigma(share, len(sky))).ok, share


def test_hot_jupiters_are_rare_and_rarer_round_a_red_dwarf():
    def rate(mass_solar, n, seed):
        _mass, _snow, sky = _sky(mass_solar, -0.03, n, seed, spread=0.2)
        return sum(1 for a in sky if any(s.migrated and s.distance_au < tuning.GIANT_HOT_ORBIT_AU for s in a.planets)) / len(sky)

    # GEN.1336: the red dwarf's rate is the pool of three seed blocks (it passed on seed 405 alone before the gas of a lighter star's disk fell with its mass)
    sun, red = rate(1.0, 6000, 404), sum(rate(0.55, 2000, seed) for seed in (405, 406, 407)) / 3.0
    assert po.compare(sun, po.BANDS["hot_jupiter_rate.G"], po.share_sigma(sun, 6000)).ok or \
        po.compare(sun, po.BANDS["hot_jupiter_rate.G.rv"], po.share_sigma(sun, 6000)).ok, sun
    assert po.compare(red, po.BANDS["hot_jupiter_rate.M"], po.share_sigma(red, 6000)).ok, red


def test_an_m_dwarf_has_inner_systems_more_often_than_a_sun_like_star():
    mass_kg, _snow, red = _sky(0.3, 0.0, 800, 505)
    _mass, _snow, sun = _sky(1.0, 0.0, 800, 506)
    red_share = sum(1 for a in red if _inner(a, mass_kg)) / len(red)
    sun_share = sum(1 for a in sun if _inner(a, SUN)) / len(sun)
    assert red_share > 1.8 * sun_share


def test_the_multiplicity_kepler_would_have_seen_matches_the_survey():
    mass_kg, _snow, sky = _sky(1.0, 0.0, 1200, 606)
    host_radius = 1.0
    systems = [{"host_mass_solar": 1.0, "host_radius_solar": host_radius,
                "planets": [{"period_d": period_d(mass_kg, s.distance_au), "radius_earth": _radius_earth(s.mass_kg / EARTH)}
                            for s in a.planets if 0.5 <= _radius_earth(s.mass_kg / EARTH) <= 20.0]} for a in sky]
    draw.set_run_seed(607)
    orientations = 25
    counts = kepler_multiplicity.transit_counts(systems, orientations=orientations)
    assert counts.get(1, 0) > counts.get(2, 0) > counts.get(3, 0)
    # the orientations of one system are not independent skies: the sample is the systems, so the counts are put per orientation
    chi2, dof = po.multiplicity_chi2({k: v / orientations for k, v in counts.items()})
    assert chi2 < 5.0 * dof + 10.0, (chi2, dof, counts)
