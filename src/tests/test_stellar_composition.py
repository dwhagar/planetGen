# tests/test_stellar_composition.py

"""
GEN.440: a star's composition derived from its stored [Fe/H], age, mass and key (`planetgen/physics/stellar_composition.py`).

The relations themselves are tested in test_stellar_abundances.py (GEN.439); here the key, the draws, the kinds and the phases.
"""

import concurrent.futures
import math

import pytest

from planetgen import tuning
from planetgen.physics import stellar_abundances as sa
from planetgen.physics import stellar_composition as sc
from planetgen.physics.stellar_evolution import evolve_star, main_sequence_lifetime_gy
from planetgen.util import draw


def total(composition):
    return math.fsum(composition.elements.values())


# --- The key and the draws ---------------------------------------------------------------------------------------------

def test_the_key_prefers_the_group_then_the_bright_seed_then_the_uid_then_the_quantities():
    assert sc.abundance_key(uid=5, mass=1.0, age_gy=2.0, distance_ly=3.0, bright_seed=7, group="9") == "9"
    assert sc.abundance_key(uid=5, mass=1.0, age_gy=2.0, distance_ly=3.0, bright_seed=7) == "bright:7"
    assert sc.abundance_key(uid=5, mass=1.0, age_gy=2.0, distance_ly=3.0) == "5"
    assert sc.abundance_key(mass=1.0, age_gy=2.0, distance_ly=3.0) == "1.0:2.0:3.0"
    assert sc.cluster_key(11) == "cluster:11"
    with pytest.raises(sc.StarMakeupError):
        sc.abundance_key(mass=1.0)


def test_the_draws_are_a_function_of_the_key_alone_and_all_are_always_drawn():
    one, two = sc.key_draws("42"), sc.key_draws("42")
    assert one == two and one != sc.key_draws("43")
    assert 0.0 <= one.u_branch < 1.0 and len(one.normals) == 4
    # The fixed order: changing only what is drawn last cannot move what is drawn first.
    stream = draw.Stream("abundance:42")
    assert stream.random() == one.u_branch and tuple(stream.gauss(0.0, 1.0) for _ in range(4)) == one.normals


def test_the_draws_are_a_stream_of_their_own():
    before = draw.Stream("metallicity:bright:7").random()
    sc.key_draws("bright:7")
    assert draw.Stream("metallicity:bright:7").random() == before, "the abundance stream does not touch the [Fe/H] stream"


def test_the_same_star_gets_the_same_composition_at_one_and_four_workers():
    stars = [(-1.5 + 0.1 * n, 1.0 + n, 0.5 + 0.05 * n, f"bright:{n}") for n in range(24)]
    serial = [sc.star_makeup(sc.MAIN_SEQUENCE, fe_h, age, mass, key) for fe_h, age, mass, key in stars]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        parallel = list(pool.map(_main_sequence_makeup, stars))
    assert parallel == serial


def _main_sequence_makeup(star):
    fe_h, age, mass, key = star
    return sc.star_makeup(sc.MAIN_SEQUENCE, fe_h, age, mass, key)


# --- The composition ---------------------------------------------------------------------------------------------------

def test_a_main_sequence_star_follows_the_relations_with_its_own_draws():
    makeup = sc.star_makeup(sc.MAIN_SEQUENCE, -0.45, 9.5, 0.9, "uid-1")
    draws = sc.key_draws("uid-1")
    branch = sa.branch_from_uniform(-0.45, 9.5, draws.u_branch)
    expected = sa.composition(-0.45, sa.abundance_offsets(-0.45, branch, draws.normals, draws.u_enhanced, draws.u_carbon))
    assert makeup.branch == branch and makeup.birth == expected and makeup.surface == makeup.birth
    assert total(makeup.birth) == pytest.approx(1.0, abs=1e-12)
    assert makeup.draws == draws


def test_a_young_star_is_thin_disc_whatever_its_number_says():
    for key in ("a", "b", "c", "d"):
        assert sc.star_makeup(sc.MAIN_SEQUENCE, -1.0, 2.0, 1.0, key).branch == sa.BRANCH_THIN


def test_editing_the_metallicity_re_derives_the_composition_with_the_branch_flipping_at_most_once():
    for key in ("bright:1", "bright:2", "bright:3", "bright:4", "bright:5", "bright:6"):
        branches = [sc.star_makeup(sc.MAIN_SEQUENCE, 0.5 - 0.01 * step, 11.0, 0.8, key).branch for step in range(200)]
        flips = sum(a != b for a, b in zip(branches, branches[1:]))
        assert flips <= 1 and (flips == 0 or branches[0] == sa.BRANCH_THIN), key


def test_the_members_of_a_group_share_the_branch_and_the_offsets_and_take_their_own_phase():
    key = sc.abundance_key(uid=1, group="primary-9")
    dwarf = sc.star_makeup(sc.MAIN_SEQUENCE, -0.3, 6.0, 1.2, key)
    giant = sc.star_makeup(sc.RED_GIANT, -0.3, 6.0, 1.2, key)
    assert dwarf.branch == giant.branch and dwarf.birth == giant.birth
    assert dwarf.surface == dwarf.birth and giant.surface != giant.birth
    assert giant.surface.y > giant.birth.y and giant.surface.elements["N"] > giant.birth.elements["N"]
    assert giant.surface.elements["C"] < giant.birth.elements["C"]


def test_the_stars_of_a_cluster_share_their_residuals_and_only_the_means_move():
    key = sc.cluster_key(33)
    one = sc.star_makeup(sc.MAIN_SEQUENCE, -1.5, 12.0, 0.8, key)
    other = sc.star_makeup(sc.MAIN_SEQUENCE, -1.6, 12.0, 0.7, key)
    assert one.draws == other.draws and one.branch == other.branch == sa.BRANCH_HIGH_ALPHA


def test_a_red_giant_has_the_dredge_up_and_a_supergiant_keeps_its_birth_composition_with_open_items():
    giant = sc.star_makeup(sc.RED_GIANT, -0.1, 4.1, 1.5, "g")
    assert giant.detail["dredge_up"].delta_y == pytest.approx(0.0165, abs=3e-3)
    assert total(giant.surface) == pytest.approx(1.0, abs=1e-12)
    supergiant = sc.star_makeup(sc.SUPERGIANT, 0.0, 0.02, 15.0, "s")
    assert supergiant.surface == supergiant.birth and supergiant.open_items


def test_the_agb_star_can_be_a_carbon_star_above_the_minimum_mass():
    assert sc.carbon_star_minimum_mass_sol(0.0125) == pytest.approx(1.425)
    assert sc.carbon_star_minimum_mass_sol(0.001) == 1.4 and sc.carbon_star_minimum_mass_sol(0.05) == 1.5
    heavy = sc.star_makeup(sc.ASYMPTOTIC_GIANT, -0.1, 0.4, 2.0, "agb")
    light = sc.star_makeup(sc.ASYMPTOTIC_GIANT, -0.1, 0.4, 1.0, "agb")
    assert heavy.detail["can_be_carbon_star"] and not light.detail["can_be_carbon_star"]


def test_a_wolf_rayet_star_is_hydrogen_free_with_carbon_inside_the_range():
    makeup = sc.star_makeup(sc.WOLF_RAYET, 0.0, 0.004, 20.0, "w", wolf_rayet="WC")
    fractions = makeup.detail["surface_fractions"]
    assert 0.20 <= fractions["C"] <= 0.55 and math.fsum(fractions.values()) == pytest.approx(1.0)
    assert fractions["C"] == pytest.approx(0.20 + 0.35 * makeup.draws.u_wolf_rayet)
    early = sc.wolf_rayet_surface("WN-early", 0.5)
    assert early["H"] == 0.0 and early["N"] == pytest.approx(0.01)
    with pytest.raises(sc.StarMakeupError):
        sc.star_makeup(sc.WOLF_RAYET, 0.0, 0.004, 20.0, "w")
    with pytest.raises(sc.StarMakeupError):
        sc.wolf_rayet_surface("WX", 0.5)


def test_a_brown_dwarf_reports_its_limits():
    makeup = sc.star_makeup(sc.BROWN_DWARF, 0.0, 1.0, 0.05, "bd")
    detail = makeup.detail
    assert detail["mass_jupiter"] == pytest.approx(52.4, abs=0.1)
    assert detail["burned_deuterium"] and not detail["burned_lithium"] and not detail["burns_hydrogen"]
    assert detail["radius_solar_at_1_gy"] == pytest.approx(0.0960)
    heavy = sc.star_makeup(sc.BROWN_DWARF, 0.0, 1.0, 0.065, "bd")
    assert heavy.detail["burned_lithium"] and not heavy.detail["burns_hydrogen"]
    light = sc.star_makeup(sc.BROWN_DWARF, 0.0, 1.0, 0.02, "bd")
    assert light.detail["radius_solar_at_1_gy"] is None and not light.detail["burned_lithium"]
    assert light.detail["burned_deuterium"] is (light.detail["mass_jupiter"] > light.detail["deuterium_limit_mj"])


def test_a_white_dwarf_has_a_core_an_atmosphere_and_a_radius():
    makeup = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.61, "wd")
    assert makeup.detail["core"] == "carbon-oxygen" and makeup.detail["radius_solar"] == pytest.approx(0.01258, rel=0.02)
    assert makeup.detail["atmosphere"] in ("DA", "non-DA")
    shares = [sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, f"wd{n}").detail["atmosphere"] == "DA" for n in range(600)]
    assert sum(shares) / 600 == pytest.approx(tuning.STELLAR_WD_DA_SHARE, abs=0.05)
    assert sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 1.2, "x").detail["core"] == "oxygen-neon"
    assert sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 1.6, "x").detail["radius_solar"] is None


def test_a_neutron_star_and_a_black_hole_store_no_metallicity_and_have_none():
    star = sc.star_makeup(sc.NEUTRON_STAR, None, 0.1, 1.4, "ns")
    assert star.birth is None and star.surface is None and star.branch is None
    assert star.detail["outer_crust"][0] == "Fe-56" and star.detail["outer_crust"][-1] == "Kr-118" and star.open_items
    hole = sc.star_makeup(sc.BLACK_HOLE, None, 0.1, 10.0, "bh", spin=0.9)
    assert hole.birth is None and hole.detail["hole"].isco_gm_c2 == pytest.approx(2.321, abs=2e-3)
    with pytest.raises(sc.StarMakeupError):
        sc.star_makeup(sc.MAIN_SEQUENCE, None, 1.0, 1.0, "x")
    with pytest.raises(sc.StarMakeupError):
        sc.star_makeup("quark star", 0.0, 1.0, 1.0, "x")


def test_the_kind_comes_from_the_stored_class():
    assert sc.kind_of("V", 1.0) == sc.MAIN_SEQUENCE and sc.kind_of("VI", 0.5) == sc.MAIN_SEQUENCE
    assert sc.kind_of("IV", 1.5) == sc.SUBGIANT and sc.kind_of("III", 1.5) == sc.RED_GIANT
    assert sc.kind_of("III", 2.0, phase="asymptotic-giant") == sc.ASYMPTOTIC_GIANT and sc.kind_of("II", 5.0) == sc.RED_GIANT
    assert sc.kind_of("IB", 20.0) == sc.SUPERGIANT and sc.kind_of("0", 80.0) == sc.SUPERGIANT
    assert sc.kind_of("VII", 0.6) == sc.WHITE_DWARF and sc.kind_of("D", 0.06) == sc.WHITE_DWARF
    assert sc.kind_of("protostar", 1.0) == sc.PROTOSTAR and sc.kind_of("PMS", 1.0) == sc.PRE_MAIN_SEQUENCE
    assert sc.kind_of("V", 0.05) == sc.BROWN_DWARF, "below the hydrogen limit whatever the class"
    with pytest.raises(sc.StarMakeupError):
        sc.kind_of("Z", 1.0)


# --- The sweep ---------------------------------------------------------------------------------------------------------

def test_a_sweep_of_every_kind_of_star_the_generator_makes_adds_up_and_follows_the_relations():
    stream = draw.Stream("gen440-sweep")
    seen = set()
    masses = [0.04, 0.06, 0.09, 0.3, 0.8, 1.0, 1.5, 2.5, 5.0, 9.0, 15.0, 25.0, 40.0, 80.0]
    ages = [0.0002, 0.005, 0.05, 0.5, 2.0, 5.0, 9.0, 11.0, 13.0]
    count = 0
    for mass in masses:
        lifetime = main_sequence_lifetime_gy(mass)
        for age in ages + [lifetime * share for share in (0.3, 0.95, 1.03, 1.1, 1.2, 1.4, 2.0)]:
            fe_h = -2.0 + 2.4 * stream.random()
            state = evolve_star(mass, age, rng=stream, fe_h=fe_h)
            key = sc.abundance_key(mass=mass, age_gy=age, distance_ly=26000.0)
            if state is None:
                kind = sc.NEUTRON_STAR if mass < 25.0 else sc.BLACK_HOLE
                makeup = sc.star_makeup(kind, None, age, 1.4 if kind == sc.NEUTRON_STAR else mass, key)
                assert makeup.birth is None
            else:
                kind = sc.kind_of(state["yerkes_class"], state["mass_sol"], state.get("phase"))
                makeup = sc.star_makeup(kind, fe_h, age, state["mass_sol"], key)
                draws = sc.key_draws(key)
                branch = sa.branch_from_uniform(fe_h, age, draws.u_branch)
                expected = sa.composition(fe_h, sa.abundance_offsets(fe_h, branch, draws.normals, draws.u_enhanced,
                                                                     draws.u_carbon))
                assert makeup.branch == branch and makeup.birth == expected, (mass, age)
                for composition in (makeup.birth, makeup.surface):
                    assert total(composition) == pytest.approx(1.0, abs=1e-12), (kind, mass, age)
                    assert composition.x + composition.y + composition.z == pytest.approx(1.0, abs=1e-12)
                    assert all(value >= 0.0 for value in composition.elements.values())
                if kind == sc.WHITE_DWARF:
                    assert makeup.detail["core"] in ("helium", "carbon-oxygen", "oxygen-neon")
            seen.add(kind)
            count += 1
    expected_kinds = {sc.MAIN_SEQUENCE, sc.SUBGIANT, sc.RED_GIANT, sc.SUPERGIANT, sc.WHITE_DWARF, sc.NEUTRON_STAR,
                      sc.BLACK_HOLE, sc.BROWN_DWARF}
    assert expected_kinds <= seen, f"the sweep of {count} stars never made {sorted(expected_kinds - seen)}"
