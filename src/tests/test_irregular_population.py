# tests/test_irregular_population.py

"""GEN.1077: the capture windows, the ground-down count and the stored irregulars of a planet (`planetgen.physics.irregular_population`)."""

import pytest

from planetgen.physics import belt_ledger, constants, irregular_zone, irregular_population as ip
from planetgen.util import draw

SUN = constants.SOLAR_MASS_TO_KG
EARTH = constants.EARTH_MASS_TO_KG
AGE_YR = 4.5e9
ENCOUNTERS = ip.REFERENCE_ENCOUNTERS
GIANTS = {"Jupiter": (5.203, 1.898e27, True), "Saturn": (9.537, 5.683e26, True), "Uranus": (19.19, 8.68e25, False), "Neptune": (30.07, 1.024e26, False)}


def _planet(distance_au, mass_kg, gas_giant, encounters=ENCOUNTERS, reservoir=ip.REFERENCE_DISC_MASS_EARTH):
    hill = irregular_zone.hill_radius_au(distance_au, mass_kg, SUN)
    pool = ip.make_pool(hill, encounters, reservoir, gas_giant, 3.0 if gas_giant else 0.0)
    return hill, pool, ip.intrinsic_probability(mass_kg, ip.zone_radius_au(hill, distance_au))


def _ground(name, radius, **kwargs):
    distance, mass, gas = GIANTS[name]
    _, pool, probability = _planet(distance, mass, gas, **kwargs)
    return ip.ground_count_above(pool, radius, probability, AGE_YR)


def test_the_source_law_has_the_papers_counts_and_inverts():
    assert ip.source_count_above(40.0) == pytest.approx(3.5e7)
    assert ip.source_count_above(40.0, 2.0) == pytest.approx(7.0e7)
    below = ip.source_count_above(10.0) / ip.source_count_above(20.0)
    assert below == pytest.approx(2.0 ** 1.8)
    above = ip.source_count_above(140.0) / ip.source_count_above(280.0)
    assert above == pytest.approx(2.0 ** 5.0)
    for radius in (2.0, 40.0, 70.0, 150.0, 600.0):
        assert ip.source_radius_at_count(ip.source_count_above(radius, 0.7), 0.7) == pytest.approx(radius)
    with pytest.raises(ip.IrregularPopulationError):
        ip.source_count_above(-1.0)


def test_mass_and_radius_round_trip_and_the_floor_is_ten_kilometres():
    assert ip.radius_km(ip.mass_kg(37.0)) == pytest.approx(37.0)
    assert ip.radius_km(ip.STORE_FLOOR_KG) == pytest.approx(10.6, abs=0.1)


def test_the_hill_volume_scales_the_yield_up_to_jupiters():
    assert ip.geometry_factor(0.35) == 1.0 and ip.geometry_factor(0.77) == 1.0
    assert ip.geometry_factor(0.175) == pytest.approx(0.125)
    earth = ip.geometry_factor(irregular_zone.hill_radius_au(1.0, EARTH, SUN))
    assert earth == pytest.approx(2.0e-5, rel=0.3)           # the design's figure
    assert ip.geometry_factor(irregular_zone.hill_radius_au(0.05, 1.898e27, SUN)) < 0.001


def test_jupiter_like_planet_gets_about_one_big_body_and_a_few_tens_above_five_km():
    big = ip.make_pool(0.355, ENCOUNTERS, 30.0).count_above(40.0)
    assert big == pytest.approx(1.1, abs=0.2)                # NDD14: about one above R = 40 km
    ground_big = _ground("Jupiter", 40.0)
    assert 0.2 < ground_big < 1.0
    small = _ground("Jupiter", 5.0)
    assert 8.0 < small < 40.0                                # the observed 8 to 11 sits at the low end; the model is within a factor 4
    assert small < ip.make_pool(0.355, ENCOUNTERS, 30.0, True, 3.0).count_above(5.0)


def test_the_four_giants_hold_the_same_order_of_bodies():
    counts = {name: _ground(name, 10.6) for name in GIANTS}
    assert max(counts.values()) < 3.0 * min(counts.values())
    observed = {"Jupiter": 8, "Saturn": 5, "Uranus": 6, "Neptune": 8}        # above R = 10.6 km (row CM6)
    for name, count in counts.items():
        assert observed[name] / 2.5 < count < observed[name] * 2.5
    for name in GIANTS:
        assert 0.1 < _ground(name, 40.0) < 2.0


def test_grinding_loses_little_where_the_zone_is_big_and_is_a_lower_bound():
    distance, mass, gas = GIANTS["Jupiter"]
    hill, pool, probability = _planet(distance, mass, gas)
    # the intrinsic probability is below the 3.8e-15 to 1.6e-14 Bottke et al. compute (row CM9) by 1.4 to 6
    for name, (d, m, g) in GIANTS.items():
        h = irregular_zone.hill_radius_au(d, m, SUN)
        p = ip.intrinsic_probability(m, ip.zone_radius_au(h, d))
        assert 1.0e-16 < p < 1.6e-14
    tc = ip.collision_time_yr(pool, 5.0, probability)
    assert 1.0e9 < tc < 2.0e10
    assert ip.surviving_fraction(0.0, tc) == 1.0
    assert ip.surviving_fraction(tc, tc) == pytest.approx(0.5)
    assert ip.ground_count_above(pool, 5.0, probability, 0.0) == pytest.approx(pool.count_above(5.0))
    older = ip.ground_count_above(pool, 5.0, probability, 2 * AGE_YR)
    assert older < ip.ground_count_above(pool, 5.0, probability, AGE_YR) < pool.count_above(5.0)


def test_a_denser_zone_grinds_more_and_the_late_limit_is_independent_of_what_was_captured():
    distance, mass, gas = GIANTS["Saturn"]
    hill = irregular_zone.hill_radius_au(distance, mass, SUN)
    p = ip.intrinsic_probability(mass, ip.zone_radius_au(hill, distance))
    light, heavy = (ip.make_pool(hill, n, 30.0) for n in (30.0, 3000.0))
    late = 1.0e13
    assert ip.ground_count_above(heavy, 5.0, p, late) / ip.ground_count_above(light, 5.0, p, late) < 2.0     # 100 times the captures, a ratio under 2 once saturated
    assert heavy.count_above(5.0) / light.count_above(5.0) == pytest.approx(100.0)
    smaller = ip.intrinsic_probability(mass, 0.5 * ip.zone_radius_au(hill, distance))
    assert smaller == pytest.approx(2.0 ** 3.5 * p)                                      # speed falls as r^-1/2 and the volume as r^3


def test_terrestrial_planets_and_hot_jupiters_hold_nothing():
    hill = irregular_zone.hill_radius_au(1.0, EARTH, SUN)
    pool = ip.make_pool(hill, ENCOUNTERS, 30.0)
    assert pool.count_above(5.0) < 1.0e-2
    p = ip.intrinsic_probability(EARTH, ip.zone_radius_au(hill, 1.0))
    for seed in range(5):
        found = ip.draw_irregulars(pool, p, AGE_YR, draw.Stream(seed))
        assert found.radii_km == () and found.above_floor == 0
    hot = irregular_zone.hill_radius_au(0.05, 1.898e27, SUN)
    pool = ip.make_pool(hot, ENCOUNTERS, 30.0, True, 3.0)
    assert pool.count_above(5.0) < 0.1 * ip.make_pool(0.355, ENCOUNTERS, 30.0, True, 3.0).count_above(5.0)
    no_disc = ip.make_pool(0.355, ENCOUNTERS, 0.0, True, 3.0)
    assert no_disc.count_above(5.0) == 0.0 and no_disc.count_above(0.5) == 0.0
    no_encounters = ip.make_pool(0.355, 0.0, 30.0)
    assert no_encounters.count_above(5.0) == 0.0


def test_the_gas_window_needs_an_envelope_a_lasting_disc_and_takes_a_size_range():
    none = ip.make_pool(0.45, 0.0, 30.0, gas_giant=False, gas_disc_myr=5.0)
    assert none.count_above(5.0) == 0.0                                           # an ice giant with no instability holds none
    short = ip.make_pool(0.45, 0.0, 30.0, gas_giant=True, gas_disc_myr=0.5)
    assert short.count_above(5.0) == 0.0
    gas = ip.make_pool(0.45, 0.0, 30.0, gas_giant=True, gas_disc_myr=3.0)
    full = ip.make_pool(0.45, ip.REFERENCE_ENCOUNTERS, 30.0)
    assert gas.count_above(5.0) == pytest.approx(ip.GAS_WINDOW_STRENGTH * full.count_above(5.0), rel=0.02)
    assert gas.count_above(600.0) == 0.0                                          # bodies over 500 km barely feel the drag
    assert gas.count_above(0.5) == pytest.approx(gas.count_above(1.0))            # smaller ones spiral in: nothing between 0.5 and 1 km
    assert gas.gas_share(5.0) == 1.0 and gas.gas_share(600.0) == 0.0
    both = ip.make_pool(0.45, ip.REFERENCE_ENCOUNTERS, 30.0, True, 3.0)
    assert both.gas_share(5.0) == pytest.approx(ip.GAS_WINDOW_STRENGTH / (1.0 + ip.GAS_WINDOW_STRENGTH))
    assert both.gas_share(600.0) == 0.0


def test_the_reservoir_scales_every_count():
    light, heavy = ip.make_pool(0.45, 100.0, 15.0), ip.make_pool(0.45, 100.0, 60.0)
    assert heavy.count_above(8.0) == pytest.approx(4.0 * light.count_above(8.0))


def test_the_exchange_is_rare_massive_and_absent_at_gas_giants():
    assert ip.exchange_chance(1.0, 0.77, gas_giant=True) == 0.0
    chance = ip.exchange_chance(1.0, 0.77, gas_giant=False)
    assert chance == pytest.approx(ip.EXCHANGE_CHANCE) and 0.002 <= chance <= 0.035
    assert ip.exchange_chance(0.2, 0.77, False) == pytest.approx(0.2 * chance)
    assert ip.exchange_chance(1.0, 0.1, False) < chance
    assert ip.draw_exchange(0.9, 0.5, chance) is None
    assert ip.draw_exchange(0.001, 0.0, chance) == ip.EXCHANGE_MIN_RADIUS_KM
    assert ip.draw_exchange(0.001, 0.9, chance) > ip.draw_exchange(0.001, 0.5, chance)
    with pytest.raises(ip.IrregularPopulationError):
        ip.exchange_chance(1.5, 0.5, False)


def test_temporary_captures_are_a_flag_for_planets_with_a_zone():
    assert ip.holds_temporary_captures(0.001, 6.371e3)                  # Earth: 150 thousand km
    assert not ip.holds_temporary_captures(0.0001, 7.0e4)               # a hot Jupiter's 4 radii inside the limit
    assert ip.holds_temporary_captures(0.14, 7.0e4)


def test_the_draw_follows_the_rank_rule_and_is_the_same_stream_whatever_the_planet():
    distance, mass, gas = GIANTS["Jupiter"]
    _, pool, probability = _planet(distance, mass, gas)
    a = ip.draw_irregulars(pool, probability, AGE_YR, draw.Stream(3))
    b = ip.draw_irregulars(pool, probability, AGE_YR, draw.Stream(3))
    assert a == b
    assert list(a.radii_km) == sorted(a.radii_km, reverse=True)
    assert a.masses_kg == tuple(ip.mass_kg(r) for r in a.radii_km)
    assert a.above_floor == sum(1 for m in a.masses_kg if m >= ip.STORE_FLOOR_KG)
    assert len(a.radii_km) >= min(a.above_floor, ip.STORE_MAX_BODIES) and len(a.radii_km) <= ip.STORE_MAX_BODIES
    assert len(a.gas_window) == len(a.radii_km)
    assert a.covered >= belt_ledger.LEDGER_COVERAGE or len(a.radii_km) == ip.STORE_MAX_BODIES
    # rank k always uses the k-th exponential: a lighter reservoir gives the same stream, a smaller body at each rank
    light = ip.make_pool(irregular_zone.hill_radius_au(distance, mass, SUN), ENCOUNTERS, 10.0, True, 3.0)
    c = ip.draw_irregulars(light, probability, AGE_YR, draw.Stream(3))
    assert c.radii_km[0] < a.radii_km[0]
    exponentials, windows = ip.draw_unit_numbers(draw.Stream(3))
    assert len(exponentials) == len(windows) == ip.STORE_MAX_BODIES
    with pytest.raises(ip.IrregularPopulationError):
        ip.irregulars_from(pool, probability, AGE_YR, exponentials[:5], windows)


def test_poisson_top_end_some_jupiters_have_none_of_the_big_ones_and_some_several():
    distance, mass, gas = GIANTS["Jupiter"]
    _, pool, probability = _planet(distance, mass, gas)
    big = []
    for seed in range(40):
        exponentials, windows = ip.draw_unit_numbers(draw.Stream(seed))
        gap, n = 0.0, 0
        for e in exponentials:
            gap += e
            if gap > ip.ground_count_above(pool, 40.0, probability, AGE_YR):
                break
            n += 1
        big.append(n)
    assert min(big) == 0 and max(big) >= 2
    assert 0.2 < sum(big) / len(big) < 0.9                # expected 0.4 above 40 km


def test_the_stored_list_of_the_four_giants_is_about_the_observed_size():
    for name, (distance, mass, gas) in GIANTS.items():
        _, pool, probability = _planet(distance, mass, gas)
        counts = [ip.draw_irregulars(pool, probability, AGE_YR, draw.Stream(s)).above_floor for s in range(8)]
        assert 2.0 <= sum(counts) / len(counts) <= 14.0          # observed 8, 5, 6, 8 above 1e16 kg
        found = ip.draw_irregulars(pool, probability, AGE_YR, draw.Stream(1))
        assert found.remainder_kg < 0.2 * sum(found.masses_kg) + 1.0e17


def test_mass_below_a_radius_is_the_integral_of_the_ground_law():
    distance, mass, gas = GIANTS["Saturn"]
    _, pool, probability = _planet(distance, mass, gas)
    total = ip.mass_below_kg(pool, probability, AGE_YR, 500.0)
    half = ip.mass_below_kg(pool, probability, AGE_YR, 50.0)
    assert 0.0 < half < total
    assert ip.mass_below_kg(pool, probability, AGE_YR, ip.MODEL_MIN_RADIUS_KM) == 0.0
    # a population of cumulative index under 3 below its break carries its mass in its largest bodies: the stored share is most of it
    assert half / total < 0.5


def test_radius_at_count_inverts_the_ground_count_and_returns_none_past_the_population():
    distance, mass, gas = GIANTS["Uranus"]
    _, pool, probability = _planet(distance, mass, gas)
    for count in (0.5, 3.0, 12.0):
        r = ip.radius_at_count(pool, probability, AGE_YR, count)
        assert ip.ground_count_above(pool, r, probability, AGE_YR) == pytest.approx(count, rel=1e-3)
    assert ip.radius_at_count(pool, probability, AGE_YR, 1.0e9) is None
