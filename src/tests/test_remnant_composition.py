"""
GEN.1271: what a supernova remnant is made of (`physics/remnant_composition.py`): the swept-up gas, the ejecta of a core-collapse
and of a Type Ia explosion, the progenitor's mass, and the Cassiopeia A check of the GEN.1268 note (section 6).
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation.phenomena.supernova_remnant import sedov_taylor_radius_ly
from planetgen.physics import imf, remnant_composition as rc, remnant_progenitors

LY_PER_PC = 3.261563777


def test_the_tables_are_the_sourced_ones():
    """[R5096: Table 3]: W7 is a 1.37 solar mass explosion of iron-group and intermediate-mass elements with no hydrogen or helium."""
    table = rc.table_elements()
    assert sum(row["w7"] for row in table.values()) == pytest.approx(1.37, abs=0.005)
    assert table["Fe"]["w7"] == pytest.approx(0.749, abs=0.002)
    assert table["O"]["type_ii_10_50"] == pytest.approx(1.80, abs=0.01)
    assert table["Fe"]["type_ii_10_50"] == pytest.approx(0.0907, abs=0.001)
    ia = rc.type_ia_masses()
    assert ia["H"] == 0.0 and ia["He"] == 0.0 and sum(ia.values()) == pytest.approx(1.37, abs=0.005)


def test_swept_up_mass_is_the_sphere_times_the_density_times_1_4_hydrogen_masses():
    # 1 parsec of radius at 1 hydrogen nucleus per cm^3: 4/3 pi pc^3 times 2.343e-24 g cm^-3 = 0.145 solar masses.
    assert rc.swept_up_mass_sol(LY_PER_PC, 1.0) == pytest.approx(0.1450, rel=0.01)
    assert rc.swept_up_mass_sol(2 * LY_PER_PC, 1.0) == pytest.approx(8 * rc.swept_up_mass_sol(LY_PER_PC, 1.0))
    assert rc.swept_up_mass_sol(LY_PER_PC, 3.0) == pytest.approx(3 * rc.swept_up_mass_sol(LY_PER_PC, 1.0))
    assert rc.swept_up_mass_sol(0.0, 5.0) == 0.0


def test_cassiopeia_a_swept_up_mass_within_30_percent():
    """The repo's Sedov law at 350 years in about 3 particles per cm^3 gives 8 solar masses; the XMM-Newton budget finds 7.9 [R5097]."""
    radius = sedov_taylor_radius_ly(350.0, 3.0)
    assert 2.0 < radius / LY_PER_PC < 3.2, "Cas A is about 2.5 pc across its radius"
    assert rc.swept_up_mass_sol(radius, 3.0) == pytest.approx(7.9, rel=0.30)


def test_cassiopeia_a_ejecta_within_35_percent_for_a_stripped_15_solar_mass_star():
    """The budget's 2.2 solar masses of ejecta [R5097] is the helium star of a Type IIb (a stripped progenitor of 15 solar masses, whose helium core
    is 4.352 in [R4251: Table 2]) less its neutron star. This is an independent check of the table, which was not tuned to it."""
    ejecta = rc.helium_core_mass_sol(15.0) - 1.5
    assert ejecta == pytest.approx(2.2, rel=0.35)


def test_the_yield_scales_average_to_one_over_the_imf_window_and_grow_with_mass():
    import numpy as np
    masses = [10.0 + 40.0 * (i + 0.5) / 2000 for i in range(2000)]
    weights = [float(w) for w in imf.density(np.array(masses))]
    for symbol in ("O", "Si", "Fe"):
        mean = sum(rc.yield_scale(symbol, m) * w for m, w in zip(masses, weights)) / sum(weights)
        assert mean == pytest.approx(1.0, abs=0.005), symbol
    assert rc.yield_scale("O", 25.0) > 2.5 * rc.yield_scale("O", 15.0)
    assert rc.yield_scale("Fe", 25.0) == rc.yield_scale("Fe", 12.0) == 1.0
    assert rc.yield_scale("Si", 25.0) > rc.yield_scale("Si", 12.0)


def test_a_core_collapse_ejecta_is_oxygen_rich_and_a_type_ia_is_iron_rich():
    cc = rc.remnant_composition("cc", "core-collapse", "neutron_star", 8.0, 3.0)
    ia = rc.remnant_composition("ia", "Type Ia", None, 8.0, 3.0)
    metals = lambda r: {k: v for k, v in r.ejecta.items() if k not in ("H", "He")}  # noqa: E731
    assert cc.ejecta["O"] > 3 * cc.ejecta["Fe"] or cc.progenitor_mass_sol < 14, (cc.progenitor_mass_sol, cc.ejecta)
    heavy = rc.processed_masses(25.0)
    assert heavy["O"] > 5 * heavy["Fe"] and heavy["O"] > heavy["Si"] > 0
    assert ia.ejecta["Fe"] > 0.5 * sum(metals(ia).values())
    assert ia.ejecta["Fe"] > 4 * ia.ejecta["Si"]
    assert ia.ejecta["O"] < 0.2 * ia.ejecta["Fe"]
    assert ia.ejecta["H"] == 0.0 and ia.ejecta["He"] == 0.0
    assert ia.compact_mass_sol == 0.0 and ia.stripped
    assert ia.ejecta_mass_sol == pytest.approx(1.37, abs=0.005)


def test_every_part_sums_and_the_fractions_sum_to_one():
    for kind, compact in (("core-collapse", "neutron_star"), ("core-collapse", "black_hole"), ("core-collapse", None), ("Type Ia", None)):
        for key in ("a", "b", "c", "d"):
            r = rc.remnant_composition(key, kind, compact, 10.0, 2.0)
            assert sum(r.ejecta.values()) == pytest.approx(r.ejecta_mass_sol, rel=1e-9)
            assert sum(r.swept.values()) == pytest.approx(r.swept_mass_sol, rel=1e-9)
            assert all(v >= 0 for v in r.ejecta.values()) and all(v >= 0 for v in r.swept.values())
            ejecta, swept, total = r.fractions()
            assert sum(ejecta.values()) == pytest.approx(1.0) and sum(swept.values()) == pytest.approx(1.0) and sum(total.values()) == pytest.approx(1.0)
            assert 0.0 < r.ejecta_share < 1.0
            assert r.ejecta_mass_sol > 0 and r.total_mass_sol == pytest.approx(r.ejecta_mass_sol + r.swept_mass_sol)


def test_a_remnant_reads_the_same_every_time_and_by_its_own_key():
    a = rc.remnant_composition("object-1", "core-collapse", "neutron_star", 8.0, 3.0)
    assert rc.remnant_composition("object-1", "core-collapse", "neutron_star", 8.0, 3.0) == a
    masses = {rc.remnant_composition(f"object-{i}", "core-collapse", "neutron_star", 8.0, 3.0).progenitor_mass_sol for i in range(20)}
    assert len(masses) == 20


def test_the_progenitor_lies_in_the_windows_of_its_stored_core():
    low, high = tuning.CORE_COLLAPSE_MASS_RANGE_SOL
    ns_windows = remnant_progenitors.kind_windows(remnant_progenitors.NEUTRON_STAR)
    bh_windows = remnant_progenitors.kind_windows(remnant_progenitors.BLACK_HOLE)
    inside = lambda m, windows: any(a <= m <= b for a, b in windows)  # noqa: E731
    for i in range(40):
        assert inside(rc.progenitor_mass_sol(f"n{i}", "core-collapse", "neutron_star"), ns_windows)
        assert inside(rc.progenitor_mass_sol(f"b{i}", "core-collapse", "black_hole"), bh_windows)
        assert low <= rc.progenitor_mass_sol(f"x{i}", "core-collapse", None) <= high
        assert 3.0 <= rc.progenitor_mass_sol(f"i{i}", "Type Ia", None) <= 8.0


def test_a_progenitor_of_40_solar_masses_or_more_is_always_stripped_and_a_third_of_the_rest_are():
    assert rc.is_stripped(41.0, 0.99) and rc.is_stripped(80.0, 0.99) and not rc.is_stripped(39.0, 0.99)
    assert not rc.is_stripped(20.0, 0.9) and rc.is_stripped(20.0, 0.1)
    stripped = [rc.remnant_composition(f"s{i}", "core-collapse", None, 8.0, 3.0) for i in range(300)]
    light = [r for r in stripped if r.progenitor_mass_sol < 40.0]
    share = sum(r.stripped for r in light) / len(light)
    assert 0.2 < share < 0.45, share
    for r in stripped:
        if r.stripped:
            assert r.ejecta["H"] == 0.0, "a stripped progenitor has no hydrogen to throw out"
        else:
            assert r.ejecta["H"] > 0.0
        if r.progenitor_mass_sol >= 40.0:
            assert r.stripped


def test_a_metal_poor_place_gives_a_metal_poor_swept_up_gas_and_a_given_gas_is_used_as_it_is():
    rich = rc.remnant_composition("m", "Type Ia", None, 10.0, 2.0, fe_h=0.0)
    poor = rc.remnant_composition("m", "Type Ia", None, 10.0, 2.0, fe_h=-1.0)
    z = lambda r: 1.0 - (r.swept["H"] + r.swept["He"]) / r.swept_mass_sol  # noqa: E731
    assert z(poor) < 0.3 * z(rich)
    given = rc.remnant_composition("m", "Type Ia", None, 10.0, 2.0, ambient={"H": 0.75, "He": 0.24, "O": 0.01})
    assert given.swept["O"] / given.swept_mass_sol == pytest.approx(0.01)
    assert given.swept["H"] / given.swept_mass_sol == pytest.approx(0.75)
    assert given.ejecta == rc.remnant_composition("m", "Type Ia", None, 10.0, 2.0).ejecta, "the gas does not change the ejecta"


def test_ejecta_are_the_final_mass_less_the_compact_core_and_a_black_hole_core_is_heavier():
    for i in range(60):
        r = rc.remnant_composition(f"e{i}", "core-collapse", "neutron_star", 8.0, 3.0)
        final = rc.helium_core_mass_sol(r.progenitor_mass_sol) if r.stripped else r.progenitor_mass_sol - rc.wind_loss_sol(r.progenitor_mass_sol)
        assert r.ejecta_mass_sol + r.compact_mass_sol == pytest.approx(min(final, r.progenitor_mass_sol), rel=1e-6)
        assert 1.0 < r.compact_mass_sol < 2.6
    holes = [rc.remnant_composition(f"h{i}", "core-collapse", "black_hole", 8.0, 3.0) for i in range(40)]
    assert all(h.compact_mass_sol > 2.5 for h in holes)
    assert all(h.ejecta_mass_sol > 0 for h in holes)


def test_the_helium_star_and_the_wind_are_the_rows_of_the_presupernova_table():
    """[R4251: Table 2]: the helium core and the final mass of each model, and the rows reproduce exactly."""
    for initial, final, core in rc._presupernova_rows():
        assert rc.helium_core_mass_sol(initial) == pytest.approx(core)
        assert rc.wind_loss_sol(initial) == pytest.approx(initial - final)
    assert rc.wind_loss_sol(15.0) == pytest.approx(2.397, abs=1e-3) and rc.wind_loss_sol(25.0) == pytest.approx(9.151, abs=1e-3)
    assert rc.helium_core_mass_sol(150.0) == pytest.approx(6.16)             # held at the heaviest model
    assert rc.wind_loss_sol(150.0) == pytest.approx(150.0 - 6.16)


def test_between_and_under_the_rows_the_table_is_continuous_and_never_gives_a_core_over_the_final_mass():
    masses = [8.0 + 0.25 * i for i in range(0, 460)]
    cores = [rc.helium_core_mass_sol(m) for m in masses]
    for m, core in zip(masses, cores):
        assert 0.0 < core <= m - rc.wind_loss_sol(m) + 1e-9, m
    assert max(abs(b - a) for a, b in zip(cores, cores[1:])) < 0.7
    assert 0.0 < rc.wind_loss_sol(8.0) < 0.5 and rc.helium_core_mass_sol(8.0) < rc.helium_core_mass_sol(10.0)


def test_a_kind_it_does_not_know_is_an_error():
    with pytest.raises(rc.RemnantCompositionError):
        rc.remnant_composition("k", "pair-instability", None, 8.0, 3.0)
    with pytest.raises(rc.RemnantCompositionError):
        rc.remnant_composition("k", "core-collapse", "white_dwarf", 8.0, 3.0)
