# tests/test_chemically_peculiar.py

"""GEN.1027: the Ap and Am stars (`planetgen/physics/chemically_peculiar.py`)."""

import pytest

from planetgen.physics import chemically_peculiar as cp


def test_the_shares_of_each_range_are_the_stand_in_shares():
    us = [(k + 0.5) / 10000.0 for k in range(10000)]
    in_both = [cp.kind_of(8500.0, u) for u in us]
    assert in_both.count(cp.AP) / len(us) == pytest.approx(cp.AP_SHARE, abs=1.0e-3)
    assert in_both.count(cp.AM) / len(us) == pytest.approx(cp.AM_SHARE, abs=1.0e-3)
    assert cp.AP not in [cp.kind_of(7200.0, u) for u in us]    # cooler than the Ap range
    assert cp.AM not in [cp.kind_of(11000.0, u) for u in us]    # hotter than the Am range
    assert {cp.kind_of(5000.0, u) for u in us} == {None} and {cp.kind_of(30000.0, u) for u in us} == {None}


def test_a_star_is_never_both():
    for temperature in (7000.0, 7500.0, 8500.0, 9500.0, 13000.0):
        kinds = [cp.kind_of(temperature, k / 1000.0) for k in range(1000)]
        assert set(kinds) <= {None, cp.AP, cp.AM}


def test_an_ap_star_rotates_in_days_inside_its_range_and_an_am_star_slower_than_100_km_s():
    periods = [cp.rotation_period_hours(cp.AP, 2.0, (k + 0.5) / 1000.0) for k in range(1000)]
    assert min(periods) >= 0.5 * 24.0 - 1.0e-9 and max(periods) <= 30.0 * 24.0 + 1.0e-9
    ordered = sorted(periods)
    assert ordered[500] / 24.0 == pytest.approx(cp.AP_PERIOD_MEDIAN_DAYS, rel=0.2)
    speed_low = 2.0 * 3.141592653589793 * 1.8 * cp.SOLAR_RADIUS_KM / (cp.rotation_period_hours(cp.AM, 1.8, 0.0) * 3600.0)
    speed_high = 2.0 * 3.141592653589793 * 1.8 * cp.SOLAR_RADIUS_KM / (cp.rotation_period_hours(cp.AM, 1.8, 1.0) * 3600.0)
    assert speed_low == pytest.approx(cp.AM_SPEED_RANGE_KMS[0]) and speed_high == pytest.approx(cp.AM_SPEED_RANGE_KMS[1])


def test_the_ap_field_is_kilogauss_scale_and_an_am_star_has_none():
    low, high = cp.AP_POLAR_FIELD_GAUSS
    assert cp.polar_field_gauss(cp.AP, 0.0) == pytest.approx(low) and cp.polar_field_gauss(cp.AP, 1.0) == pytest.approx(high)
    assert low < cp.polar_field_gauss(cp.AP, 0.5) < high
    assert cp.polar_field_gauss(cp.AM, 0.5) is None
    with pytest.raises(ValueError):
        cp.polar_field_gauss("Fm", 0.5)
    with pytest.raises(ValueError):
        cp.rotation_period_hours("Fm", 1.0, 0.5)


def _stars(temperature_k, count, yerkes_class="V"):
    return [cp.of_star(yerkes_class, 3.0e30 * (1.0 + k * 1.0e-3), 1.8 * cp.SOLAR_RADIUS_KM, temperature_k, 0.3 + k * 1.0e-3) for k in range(count)]


def test_of_star_gives_the_shares_and_the_same_answer_every_time():
    found = _stars(8500.0, 4000)
    kinds = [p.kind for p in found if p]
    assert kinds.count(cp.AP) / len(found) == pytest.approx(cp.AP_SHARE, abs=0.02)
    assert kinds.count(cp.AM) / len(found) == pytest.approx(cp.AM_SHARE, abs=0.03)
    assert found == _stars(8500.0, 4000)    # a keyed stream, not a draw of the run


def test_only_a_main_sequence_dwarf_in_range_is_peculiar():
    assert not any(_stars(8500.0, 300, "III")) and not any(_stars(8500.0, 300, "VI")) and not any(_stars(8500.0, 300, "D"))
    assert not any(_stars(5800.0, 300)) and not any(_stars(25000.0, 300))
    assert cp.of_star("V", None, 1.0e6, 8500.0, 0.3) is None


def test_the_period_and_field_of_a_found_star_are_those_of_its_kind():
    for found in _stars(8500.0, 800):
        if found is None:
            continue
        if found.kind == cp.AP:
            assert 0.5 * 24.0 - 1.0e-9 <= found.period_hours <= 30.0 * 24.0 + 1.0e-9
            assert cp.AP_POLAR_FIELD_GAUSS[0] <= found.field_gauss <= cp.AP_POLAR_FIELD_GAUSS[1]
        else:
            speed = 2.0 * 3.141592653589793 * 1.8 * cp.SOLAR_RADIUS_KM / (found.period_hours * 3600.0)
            assert cp.AM_SPEED_RANGE_KMS[0] - 1.0e-9 <= speed <= cp.AM_SPEED_RANGE_KMS[1] + 1.0e-9
            assert found.field_gauss is None


def test_the_spin_cascade_gives_a_peculiar_dwarf_its_kinds_period_and_leaves_the_draws_of_the_run_alone():
    from planetgen.physics import spin
    from planetgen.util import draw
    mass_kg, radius_km, age_gy = 3.0e30, 1.8 * cp.SOLAR_RADIUS_KM, 0.4
    peculiar = [(m, p) for m, p in ((mass_kg * (1.0 + k * 1.0e-3), _of(mass_kg * (1.0 + k * 1.0e-3), radius_km, age_gy)) for k in range(200)) if p]
    assert peculiar
    for mass, found in peculiar[:20]:
        with draw.bound(1):
            period = spin.star_rotation_period_hours(mass, radius_km, 8500.0, age_gy, "V")
            after = draw.random()
        with draw.bound(1):
            assert draw.random() == after    # no draw of the run was made for it
        assert period == pytest.approx(max(found.period_hours, spin.period_hours_from_speed(radius_km, spin.BREAKUP_FRACTION * spin.breakup_speed_kms(mass, radius_km))))
    ordinary = next(mass_kg * (1.0 + k * 1.0e-3) for k in range(200) if not _of(mass_kg * (1.0 + k * 1.0e-3), radius_km, age_gy))
    with draw.bound(1):
        spin.star_rotation_period_hours(ordinary, radius_km, 8500.0, age_gy, "V")
        used = draw.random()
    with draw.bound(1):
        assert draw.random() != used    # an ordinary dwarf still draws its speed


def _of(mass_kg, radius_km, age_gy):
    return cp.of_star("V", mass_kg, radius_km, 8500.0, age_gy)
