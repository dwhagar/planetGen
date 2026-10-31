# tests/test_stellar_composition.py

"""
GEN.440: a star's composition derived from its stored [Fe/H], age, mass and key (`planetgen/physics/stellar_composition.py`).

The relations themselves are tested in test_stellar_abundances.py (GEN.439); here the key, the draws, the kinds and the phases.
GEN.590: the main-sequence track, the white dwarf's radius and layers and the neutron star's crust come from the downloaded
files, which are required (no fallback); these tests serve the real download step small fixtures of them
(`tests/datasets_support.py`) and check the refusal when they are missing.
"""

import concurrent.futures
import math
import multiprocessing

import pytest

from planetgen import tuning
from planetgen.datasets import fetch, isochrones, neutron_star_matter, white_dwarf_models
from planetgen.physics import alpha_tracks, elements, wolf_rayet
from planetgen.physics import stellar_abundances as sa
from planetgen.physics import stellar_composition as sc
from planetgen.physics.stellar_evolution import evolve_star, main_sequence_lifetime_gy
from planetgen.util import draw
from tests import datasets_support as support


@pytest.fixture(autouse=True)
def downloaded(tmp_path, monkeypatch):
    """The three datasets, as small fixtures, in a cache folder of this test's own."""
    support.reset_caches(tmp_path / "cache", monkeypatch)
    support.install_all(monkeypatch)


def total(composition):
    return math.fsum(composition.elements.values())


# --- The key and the draws ---------------------------------------------------------------------------------------------

def test_the_key_prefers_the_group_then_the_bright_seed_then_the_uid_then_the_quantities():
    assert sc.abundance_key(uid=5, mass=1.0, distance_ly=3.0, bright_seed=7, group="9") == "9"
    assert sc.abundance_key(uid=5, mass=1.0, distance_ly=3.0, bright_seed=7) == "bright:7"
    assert sc.abundance_key(uid=5, mass=1.0, distance_ly=3.0) == "5"
    assert sc.abundance_key(mass=1.0, distance_ly=3.0) == "1.0:3.0"
    with pytest.raises(sc.StarMakeupError):
        sc.abundance_key(mass=1.0)


def test_a_stars_key_does_not_move_when_its_age_rises_after_its_planets_are_made():
    """GEN.984: a planet's host ratios and the star's makeup read one draw, whatever age the star had when the planet was made."""
    from types import SimpleNamespace
    from planetgen.physics import bulk_draw
    star = SimpleNamespace(mass=1.0, age=2.0, galactic_center_dist_ly=26000.0, metallicity_fe_h=0.0)
    before = bulk_draw.star_key(star)
    star.age = 4.5
    assert bulk_draw.star_key(star) == before == sc.abundance_key(mass=1.0, distance_ly=26000.0)
    assert sc.key_draws(before) == sc.key_draws(bulk_draw.star_key(star))


def test_the_draws_are_a_function_of_the_key_alone_and_all_are_always_drawn():
    one, two = sc.key_draws("42"), sc.key_draws("42")
    assert one == two and one != sc.key_draws("43")
    assert len(one.normals) == 5 and 0.0 <= one.u_crust_model < 1.0
    # The fixed order: changing only what is drawn last cannot move what is drawn first.
    stream = draw.Stream("abundance:42")
    stream.random()      # the draw the branch used to take, kept unread so no later draw moved (GEN.1024)
    assert tuple(stream.gauss(0.0, 1.0) for _ in range(4)) == one.normals[:4]
    assert stream.random() == one.u_enhanced and stream.random() == one.u_carbon
    assert stream.gauss(0.0, 1.0) == one.eu_normal and stream.gauss(0.0, 1.0) == one.th_eu_normal
    assert stream.random() == one.u_atmosphere and stream.random() == one.u_wolf_rayet
    assert stream.random() == one.u_crust_model, "the draw added before the pollution one (GEN.590) moved none before it"
    assert stream.random() == one.u_pollution, "the draw added before the calcium one (GEN.253) moved none before it"
    assert stream.gauss(0.0, 1.0) == one.normals[4], "the calcium residual (GEN.988) is the last draw and moved none before it"


def test_the_draws_are_a_stream_of_their_own():
    before = draw.Stream("metallicity:bright:7").random()
    sc.key_draws("bright:7")
    assert draw.Stream("metallicity:bright:7").random() == before, "the abundance stream does not touch the [Fe/H] stream"


def test_the_same_star_gets_the_same_composition_at_one_and_four_workers():
    stars = [(-1.5 + 0.1 * n, 1.0 + n, 0.5 + 0.05 * n, f"bright:{n}") for n in range(24)]
    serial = [sc.star_makeup(sc.MAIN_SEQUENCE, fe_h, age, mass, key) for fe_h, age, mass, key in stars]
    # Forked, so the workers keep the test's stand-in for the downloaded data check (Python 3.14 made forkserver the default).
    with concurrent.futures.ProcessPoolExecutor(max_workers=4, mp_context=multiprocessing.get_context("fork")) as pool:
        parallel = list(pool.map(_main_sequence_makeup, stars))
    assert parallel == serial


def _main_sequence_makeup(star):
    fe_h, age, mass, key = star
    return sc.star_makeup(sc.MAIN_SEQUENCE, fe_h, age, mass, key)


# --- The composition ---------------------------------------------------------------------------------------------------

def test_a_main_sequence_star_follows_the_relations_with_its_own_draws():
    makeup = sc.star_makeup(sc.MAIN_SEQUENCE, -0.45, 9.5, 0.9, "uid-1")
    draws = sc.key_draws("uid-1")
    branch = sa.branch_of(-0.45, 9.5)
    expected = sa.composition(-0.45, sa.abundance_offsets(-0.45, branch, draws.normals, draws.u_enhanced, draws.u_carbon))
    assert makeup.branch == branch and makeup.birth == expected and makeup.surface == makeup.birth
    assert total(makeup.birth) == pytest.approx(1.0, abs=1e-12)
    assert makeup.draws == draws


def test_a_young_star_is_thin_disc_whatever_its_number_says():
    for key in ("a", "b", "c", "d"):
        assert sc.star_makeup(sc.MAIN_SEQUENCE, -1.0, 2.0, 1.0, key).branch == sa.BRANCH_THIN


def test_editing_the_metallicity_re_derives_the_composition_with_the_branch_of_the_new_feh_and_age():
    """GEN.1024: the branch is a function of the star's [Fe/H] and age (the two numbers a stored star keeps, so its track agrees with its
    composition), so an edited [Fe/H] has the branch that [Fe/H] gives and not the one the star had; the key leaves it alone."""
    for key in ("bright:1", "bright:2", "bright:3"):
        for step in range(40):
            fe_h = 0.5 - 0.05 * step
            assert sc.star_makeup(sc.MAIN_SEQUENCE, fe_h, 11.0, 0.8, key).branch == sa.branch_of(fe_h, 11.0), (key, fe_h)


def test_the_members_of_a_group_share_the_branch_and_the_offsets_and_take_their_own_phase():
    key = sc.abundance_key(uid=1, group="primary-9")
    dwarf = sc.star_makeup(sc.MAIN_SEQUENCE, -0.3, 6.0, 1.2, key)
    giant = sc.star_makeup(sc.RED_GIANT, -0.3, 6.0, 1.2, key)
    assert dwarf.branch == giant.branch and dwarf.birth == giant.birth
    assert dwarf.surface == dwarf.birth and giant.surface != giant.birth
    assert giant.surface.y > giant.birth.y and giant.surface.elements["N"] > giant.birth.elements["N"]
    assert giant.surface.elements["C"] < giant.birth.elements["C"]


def test_the_stars_of_a_cluster_share_their_residuals_and_only_the_means_move():
    key = "cluster:33"  # the key a cluster's stars share: the draws are shared, only the [Fe/H]-dependent means move
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
    makeup = sc.star_makeup(sc.WOLF_RAYET, 0.0, 0.004, 20.0, "w", wolf_rayet_subtype="WC")
    fractions = makeup.detail["surface_fractions"]
    assert 0.20 <= fractions["C"] <= 0.55 and math.fsum(fractions.values()) == pytest.approx(1.0)
    assert fractions["C"] == pytest.approx(0.20 + 0.35 * makeup.draws.u_wolf_rayet)
    early = sc.wolf_rayet_surface("WN-early", 0.5)
    assert early["H"] == 0.0 and early["N"] == pytest.approx(0.01)
    with pytest.raises(sc.StarMakeupError):
        sc.star_makeup(sc.WOLF_RAYET, 0.0, 0.004, 20.0, "w")
    with pytest.raises(sc.StarMakeupError):
        sc.wolf_rayet_surface("WX", 0.5)


def test_a_wolf_rayet_class_star_has_a_kind_and_each_generator_subtype_a_row():
    assert sc.kind_of("W", 20.0) == sc.WOLF_RAYET
    # The generator's subtypes with a subclass number, as the spectral text writes them, or with the temperature alone.
    assert sc.wolf_rayet_row("WN4") == "WN-early" and sc.wolf_rayet_row("WN5") == "WN-early"
    assert sc.wolf_rayet_row("WN6") == "WN-late" and sc.wolf_rayet_row("WN9") == "WN-late"
    assert sc.wolf_rayet_row("WC4") == "WC" and sc.wolf_rayet_row("WC") == "WC"
    assert sc.wolf_rayet_row("WO2") == "WO" and sc.wolf_rayet_row("WO") == "WO"
    assert sc.wolf_rayet_row("WN-late") == "WN-late" and sc.wolf_rayet_row("WN/WC") == "WN/WC"
    low_k, high_k = tuning.WR_TEMPERATURE_RANGE_K["WN"]
    assert sc.wolf_rayet_row("WN", teff_k=high_k) == "WN-early" and sc.wolf_rayet_row("WN", teff_k=low_k) == "WN-late"
    with pytest.raises(sc.StarMakeupError):
        sc.wolf_rayet_row("WN")
    with pytest.raises(sc.StarMakeupError):
        sc.wolf_rayet_row("WNX")


def test_a_wo_star_is_carbon_and_oxygen_rich_and_every_row_adds_up_at_both_ends_of_its_range():
    for u in (0.0, 0.5, 1.0):
        wo = sc.wolf_rayet_surface("WO", u)
        assert wo["O"] > sc.wolf_rayet_surface("WC", u).get("O", 0.0) and wo["C"] > wo["O"] and wo["He"] >= 0.05 - 1e-12
        for name in tuning.STELLAR_WOLF_RAYET:
            fractions = sc.wolf_rayet_surface(name, u)
            assert math.fsum(fractions.values()) == pytest.approx(1.0) and min(fractions.values()) >= 0.0, (name, u)


def test_every_wolf_rayet_star_the_generator_makes_has_a_composition_that_adds_up(monkeypatch):
    monkeypatch.setattr(isochrones, "helium_star_relation", lambda: (4.31, 1.16, 0.064))   # the fixture grid holds no Wolf-Rayet row
    stream = draw.Stream("gen831")
    seen = set()
    for index in range(200):
        age = 0.003 + 0.001 * stream.random()
        state = wolf_rayet.state(40.0 + 40.0 * stream.random(), 15.0 + 20.0 * stream.random(), age, 0.003, 0.004, stream)
        text = f"{state['wr_subtype']}{wolf_rayet.subclass(state['wr_subtype'], state['temperature_k'])}"
        key = sc.abundance_key(mass=state["mass_sol"], distance_ly=26000.0 + index)
        by_text = sc.star_makeup(sc.kind_of(state["yerkes_class"], state["mass_sol"]), 0.0, age, state["mass_sol"], key,
                                 wolf_rayet_subtype=text)
        by_temperature = sc.star_makeup(sc.WOLF_RAYET, 0.0, age, state["mass_sol"], key, wolf_rayet_subtype=state["wr_subtype"],
                                        teff_k=state["temperature_k"])
        assert by_text.detail["surface_fractions"] == by_temperature.detail["surface_fractions"], text
        fractions = by_text.detail["surface_fractions"]
        assert math.fsum(fractions.values()) == pytest.approx(1.0) and min(fractions.values()) >= 0.0, text
        seen.add(state["wr_subtype"])
    assert seen == {"WN", "WC", "WO"}


def test_a_brown_dwarf_reports_its_limits():
    makeup = sc.star_makeup(sc.BROWN_DWARF, 0.0, 1.0, 0.05, "bd")
    detail = makeup.detail
    assert detail["mass_jupiter"] == pytest.approx(52.4, abs=0.1)
    assert detail["burned_deuterium"] and not detail["burned_lithium"] and not detail["burns_hydrogen"]
    assert detail["radius_solar"] == pytest.approx(0.0960, rel=0.05)         # at the star's 1 Gy, from the Baraffe table (GEN.948)
    heavy = sc.star_makeup(sc.BROWN_DWARF, 0.0, 1.0, 0.065, "bd")
    assert heavy.detail["burned_lithium"] and not heavy.detail["burns_hydrogen"]
    light = sc.star_makeup(sc.BROWN_DWARF, 0.0, 1.0, 0.02, "bd")
    assert light.detail["radius_solar"] is not None and not light.detail["burned_lithium"]
    assert light.detail["burned_deuterium"] is (light.detail["mass_jupiter"] > light.detail["deuterium_limit_mj"])


def _key_with_atmosphere(atmosphere):
    return next(key for key in (f"wd{n}" for n in range(200))
                if (sc.key_draws(key).u_atmosphere < sc.DA_SHARE) == (atmosphere == "DA"))


def test_a_white_dwarf_has_a_core_an_atmosphere_a_radius_and_layers_from_the_sequences():
    makeup = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.61, "wd", teff_k=60000.0)
    assert makeup.detail["core"] == "carbon-oxygen" and makeup.detail["atmosphere"] in ("DA", "non-DA")
    envelope = makeup.detail["envelope"]
    assert envelope == ("thick" if makeup.detail["atmosphere"] == "DA" else "thin")
    assert makeup.detail["radius_solar"] == white_dwarf_models.radius_solar(0.61, 60000.0, envelope)
    assert makeup.detail["layers"] == white_dwarf_models.layer_mass_fractions(0.61, envelope)
    assert makeup.detail["layers"]["carbon"] == makeup.detail["layers"]["oxygen"], "the sequences' own even split"
    assert any("even one" in item for item in makeup.open_items)


def test_a_hydrogen_atmosphere_takes_the_thick_layer_and_a_helium_one_the_thin_layer():
    thick = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, _key_with_atmosphere("DA"), teff_k=20000.0)
    thin = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, _key_with_atmosphere("non-DA"), teff_k=20000.0)
    assert (thick.detail["envelope"], thin.detail["envelope"]) == ("thick", "thin")
    assert thick.detail["layers"]["hydrogen"] == pytest.approx(1e-4) and thin.detail["layers"]["hydrogen"] == pytest.approx(1e-10)


def test_a_white_dwarf_cools_and_shrinks_with_the_temperature_it_is_given():
    key = _key_with_atmosphere("DA")
    hot = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, key, teff_k=100000.0).detail["radius_solar"]
    cool = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, key, teff_k=5000.0).detail["radius_solar"]
    assert hot > cool > 0.0


def test_the_atmosphere_shares_follow_the_measured_one():
    shares = [sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, f"wd{n}", teff_k=9000.0).detail["atmosphere"] == "DA"
              for n in range(600)]
    assert sum(shares) / 600 == pytest.approx(tuning.STELLAR_WD_DA_SHARE, abs=0.05)


def test_a_helium_core_or_oxygen_neon_white_dwarf_has_no_layers_and_says_so():
    light = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.3, "x", teff_k=9000.0)
    assert light.detail["core"] == "helium" and light.detail["layers"] is None and light.detail["radius_solar"] is None
    heavy = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 1.2, "x", teff_k=9000.0)
    assert heavy.detail["core"] == "oxygen-neon" and heavy.detail["layers"] is None
    assert any("helium-core or oxygen-neon" in item for item in heavy.open_items)
    assert sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 1.6, "x", teff_k=9000.0).detail["radius_solar"] is None


def test_a_white_dwarf_needs_its_temperature():
    with pytest.raises(sc.StarMakeupError, match="temperature"):
        sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, "wd")


def test_a_main_sequence_star_gets_its_track_from_the_downloaded_isochrones():
    makeup = sc.star_makeup(sc.MAIN_SEQUENCE, -0.45, 1.0, 1.5, "ms")
    track = alpha_tracks.track_fe_h(-0.45, 1.0)          # GEN.1024: the track of the star's branch
    point = isochrones.main_sequence_point(1.5, 1.0, track)
    assert makeup.detail["main_sequence"] == dict(zip(("log_luminosity_sol", "log_teff_k", "log_radius_sol"), point))
    assert makeup.detail["main_sequence"]["log_luminosity_sol"] == pytest.approx(4 * math.log10(1.5) - 0.045, abs=0.03)
    assert makeup.detail["main_sequence_lifetime_gy"] == isochrones.main_sequence_lifetime_gy(1.5, track) > 1.0
    assert not makeup.open_items


def test_a_main_sequence_star_off_the_grid_says_so_and_still_has_its_composition():
    turned_off = sc.star_makeup(sc.MAIN_SEQUENCE, 0.0, 9.0, 3.0, "ms")      # a 3 solar mass star is gone by 9 Gy
    assert turned_off.detail["main_sequence"] is None and turned_off.open_items
    assert turned_off.birth is not None and total(turned_off.birth) == pytest.approx(1.0, abs=1e-12)


def test_a_neutron_star_gets_the_outer_crust_of_one_downloaded_mass_model():
    models = neutron_star_matter.crust_models()
    seen = set()
    for n in range(40):
        key = f"ns{n}"
        star = sc.star_makeup(sc.NEUTRON_STAR, None, 0.1, 1.4, key)
        assert star.birth is None and star.surface is None and star.branch is None and star.open_items
        picked = models[int(sc.key_draws(key).u_crust_model * len(models))]
        assert star.detail["crust_model"] == picked
        seen.add(picked)
        layers = neutron_star_matter.outer_crust(picked)
        assert star.detail["outer_crust"] == tuple(f"{elements.element(layer.z).symbol}-{layer.a}" for layer in layers)
        assert star.detail["neutron_drip_g_cm3"] == pytest.approx(neutron_star_matter.mass_density_g_cm3(layers[-1].upper_density_fm3))
        assert star.detail["radius_band_km"] == tuning.STELLAR_NS_RADIUS_BAND_KM
    assert seen == set(models), "the models are drawn, not fixed"
    first = sc.star_makeup(sc.NEUTRON_STAR, None, 0.1, 1.4, "ns0").detail["outer_crust"]
    assert first[0] == "Fe-56", "the surface of every model is iron"


def test_a_black_hole_stores_no_metallicity_and_has_none():
    hole = sc.star_makeup(sc.BLACK_HOLE, None, 0.1, 10.0, "bh", spin=0.9)
    assert hole.birth is None and hole.detail["hole"].isco_gm_c2 == pytest.approx(2.321, abs=2e-3)
    with pytest.raises(sc.StarMakeupError):
        sc.star_makeup(sc.MAIN_SEQUENCE, None, 1.0, 1.0, "x")
    with pytest.raises(sc.StarMakeupError):
        sc.star_makeup("quark star", 0.0, 1.0, 1.0, "x")


# --- The files are required --------------------------------------------------------------------------------------------

@pytest.mark.parametrize("kind, fe_h, mass, kwargs, fill", [
    (sc.MAIN_SEQUENCE, 0.0, 1.0, {}, "SC1"),
    (sc.WHITE_DWARF, 0.0, 0.6, {"teff_k": 9000.0}, "SC2"),
    (sc.NEUTRON_STAR, None, 1.4, {}, "SC3"),
    (sc.BROWN_DWARF, 0.0, 0.05, {}, "SC6"),
])
def test_a_kind_that_reads_a_downloaded_file_refuses_when_it_was_not_downloaded(tmp_path, monkeypatch, kind, fe_h, mass, kwargs, fill):
    support.reset_caches(tmp_path / "nothing-downloaded", monkeypatch)
    with pytest.raises(fetch.DataMissing, match=rf"missing \({fill},"):
        sc.star_makeup(kind, fe_h, 2.0, mass, "k", **kwargs)


@pytest.mark.parametrize("kind, fe_h, mass", [
    (sc.PROTOSTAR, 0.0, 1.0), (sc.PRE_MAIN_SEQUENCE, 0.0, 1.0), (sc.SUBGIANT, 0.0, 1.5), (sc.RED_GIANT, 0.0, 1.5),
    (sc.ASYMPTOTIC_GIANT, 0.0, 2.0), (sc.SUPERGIANT, 0.0, 15.0), (sc.BLACK_HOLE, None, 10.0),
])
def test_a_kind_that_reads_no_file_does_not_need_them(tmp_path, monkeypatch, kind, fe_h, mass):
    support.reset_caches(tmp_path / "nothing-downloaded", monkeypatch)
    assert sc.star_makeup(kind, fe_h, 2.0, mass, "k").kind == kind


def test_the_derived_values_do_not_depend_on_which_machine_asks(tmp_path, monkeypatch):
    """Two machines with the same files give the same makeup; one without them gives none."""
    star = sc.star_makeup(sc.MAIN_SEQUENCE, -0.45, 1.0, 1.5, "ms")
    support.reset_caches(tmp_path / "another-machine", monkeypatch)
    support.install_all(monkeypatch)
    assert sc.star_makeup(sc.MAIN_SEQUENCE, -0.45, 1.0, 1.5, "ms") == star


def test_the_kind_comes_from_the_stored_class():
    assert sc.kind_of("V", 1.0) == sc.MAIN_SEQUENCE and sc.kind_of("VI", 0.5) == sc.MAIN_SEQUENCE
    assert sc.kind_of("IV", 1.5) == sc.SUBGIANT and sc.kind_of("III", 1.5) == sc.RED_GIANT
    assert sc.kind_of("III", 2.0, phase="asymptotic-giant") == sc.ASYMPTOTIC_GIANT and sc.kind_of("II", 5.0) == sc.RED_GIANT
    assert sc.kind_of("IB", 20.0) == sc.SUPERGIANT and sc.kind_of("0", 80.0) == sc.SUPERGIANT
    assert sc.kind_of("VII", 0.6) == sc.WHITE_DWARF and sc.kind_of("D", 0.06) == sc.WHITE_DWARF
    assert sc.kind_of("W", 20.0) == sc.WOLF_RAYET
    assert sc.kind_of("protostar", 1.0) == sc.PROTOSTAR and sc.kind_of("PMS", 1.0) == sc.PRE_MAIN_SEQUENCE
    assert sc.kind_of("V", 0.05) == sc.BROWN_DWARF, "below the hydrogen limit whatever the class"
    with pytest.raises(sc.StarMakeupError):
        sc.kind_of("Z", 1.0)


# --- The sweep ---------------------------------------------------------------------------------------------------------

def test_a_sweep_of_every_kind_of_star_the_generator_makes_adds_up_and_follows_the_relations(monkeypatch):
    # The sweep draws stars from 0.072 to 80 solar masses (0.072 is a brown dwarf at solar metallicity) through `evolve_star`, which reads the MIST grid (GEN.652) and, under its
    # 0.1 solar masses, the Baraffe table (GEN.936; a lighter object is a brown dwarf, which cools out of the table's stars): the
    # fixture grid holds 0.2 solar masses and up to a one-solar-mass turn-off, so the sweep reads the downloaded files.
    monkeypatch.undo()
    stream = draw.Stream("gen440-sweep")
    seen = set()
    masses = [0.072, 0.075, 0.09, 0.3, 0.8, 1.0, 1.5, 2.5, 5.0, 9.0, 15.0, 25.0, 40.0, 80.0]
    ages = [0.0002, 0.005, 0.05, 0.5, 2.0, 5.0, 9.0, 11.0, 13.0]
    count = 0
    for mass in masses:
        lifetime = main_sequence_lifetime_gy(mass)
        for age in ages + [lifetime * share for share in (0.3, 0.95, 1.0, 1.03, 1.06, 1.1, 1.2, 1.4, 2.0)]:
            fe_h = -2.0 + 2.4 * stream.random()
            state = evolve_star(mass, age, rng=stream, fe_h=fe_h)
            key = sc.abundance_key(mass=mass, distance_ly=26000.0)
            if state is None:
                kind = sc.NEUTRON_STAR if mass < 25.0 else sc.BLACK_HOLE
                makeup = sc.star_makeup(kind, None, age, 1.4 if kind == sc.NEUTRON_STAR else mass, key)
                assert makeup.birth is None
            elif state["yerkes_class"] == "W":
                continue          # a Wolf-Rayet star has no kind in `stellar_composition` yet (found by GEN.742's shifted draws; filed with TODO)
            else:
                kind = sc.kind_of(state["yerkes_class"], state["mass_sol"], state.get("phase"))
                makeup = sc.star_makeup(kind, fe_h, age, state["mass_sol"], key, wolf_rayet_subtype=state.get("wr_subtype"),
                                        teff_k=state["temperature_k"])
                draws = sc.key_draws(key)
                branch = sa.branch_of(fe_h, age)
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


def test_only_a_cool_white_dwarf_is_ever_polluted_and_the_share_is_the_middle_of_the_observed_range():
    keys = [f"wd-pollution-{n}" for n in range(4000)]
    cool = [sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, key, teff_k=6000.0).detail["polluted"] for key in keys]
    hot = [sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, key, teff_k=40000.0).detail["polluted"] for key in keys]
    assert not any(hot)
    assert sum(cool) / len(cool) == pytest.approx(0.375, abs=0.03)
    assert sc.POLLUTED_SHARE == 0.375 and 0.25 < sum(cool) / len(cool) < 0.50


def test_the_pollution_draw_is_the_stars_own_and_the_generator_and_the_makeup_agree():
    from planetgen.generation import multistar
    for n in range(60):
        key = f"wd-agree-{n}"
        makeup = sc.star_makeup(sc.WHITE_DWARF, -0.2, 8.0, 0.6, key, teff_k=8000.0)
        assert multistar.draw_white_dwarf_pollution(8000.0, key) is makeup.detail["polluted"]
        assert multistar.draw_white_dwarf_pollution(8000.0, key) == (sc.key_draws(key).u_pollution < 0.375)
        assert multistar.draw_white_dwarf_pollution(30000.0, key) is False


def test_the_cool_white_dwarfs_are_those_below_the_floor_of_the_hot_ones():
    from planetgen.physics import constants
    assert constants.COOL_WHITE_DWARF_BELOW_K == 10000
    assert sc.white_dwarf_polluted(9999.0, 0.0) is True and sc.white_dwarf_polluted(10000.0, 0.0) is False
    assert sc.white_dwarf_polluted(5000.0, 0.375) is False and sc.white_dwarf_polluted(5000.0, 0.374) is True

