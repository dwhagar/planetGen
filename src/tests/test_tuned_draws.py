# tests/test_tuned_draws.py

"""
The tuned draws (GEN.509 to GEN.525): every target of the tuning-target table that names the body it is drawn for and needs no selection
model is checked against the shipped recipe (`tuning_harness.check_target`: the row drawn at 2000 to 3000 seeds, then compared within 2
sigma), the knobs that were retuned are pinned, and a retune away from the data fails.
"""

import copy

import pytest

from planetgen.generation import tuning_harness as th
from planetgen.physics.atmosphere_recipes import RECIPES, Recipes

CHECKED = [t.id for t in th.TARGETS if t.role == "target" and t.subject is not None and t.selection_model is None]


def retuned(row_id, mutate):
    document = copy.deepcopy(RECIPES.document)
    mutate(next(e for e in document["entries"] if e["id"] == row_id))
    return Recipes(document, register=False)


@pytest.mark.parametrize("target_id", CHECKED)
def test_the_shipped_recipe_reproduces_the_published_statistic(target_id):
    target = th.TARGETS.get(target_id)
    result = th.check_target(target, n=3000)
    assert result.extrapolated == (), "the subject lies inside the validity range of the data"
    assert result.passed, (target_id, result.simulated, target.value, result.z)


def test_every_target_with_a_subject_is_checked_here_or_through_its_selection_model():
    assert CHECKED, "no target is checked"
    for target in th.TARGETS:
        if target.subject is not None:
            assert target.selection_model is None or target.selection_model in th.SELECTION_MODEL_ITEMS


# --- giants (GEN.516) --------------------------------------------------------------------------------------------------------


def test_the_giant_metallicity_follows_the_later_larger_fit_and_not_the_first_relation():
    spec = RECIPES.entry("f5").spec["dist"]
    assert spec["mu"] == "log10(2.6/mass_jupiter + 3.3)" and spec["sigma"] == 0.26
    old = retuned("f5", lambda e: e["dist"].update(mu="log10(9.7*mass_jupiter**-0.45)"))
    assert not th.check_target(th.TARGETS.get("f5-z-1p0"), recipes=old, n=3000).passed, "the first relation is 1.6 times the fit at one M_J"
    assert th.check_target(th.TARGETS.get("f5-z-0p1"), recipes=old, n=3000).passed, "and the two meet at 0.1 M_J"


def test_the_giant_core_width_comes_from_the_population_and_not_from_one_planet():
    assert RECIPES.entry("f7").spec["dist"]["core"]["sigma"] == 0.21
    narrow = retuned("f7", lambda e: e["dist"]["core"].update(sigma=0.08))
    assert not th.check_target(th.TARGETS.get("f7-low"), recipes=narrow, n=3000).passed, "12 of 69 planets are below 0.35"
    wide_tail = retuned("f7", lambda e: e["constants"].update(f_tail=[0.3, 0.34, "log"]))
    assert not th.check_target(th.TARGETS.get("f7-tail"), recipes=wide_tail, n=3000).passed, "the old tail bracket is outside the data"


# --- the harness's statistic fields ------------------------------------------------------------------------------------------


def test_a_log_scale_a_threshold_and_a_subject_are_read_from_the_target():
    log_sd = th.TARGETS.get("f5-width")
    assert log_sd.scale == "log10" and log_sd.subject == {"body_class": "giant", "quantities": {"mass_jupiter": 1.0}}
    assert th.compare(log_sd, [1.0, 10.0, 100.0, 10.0]).simulated == pytest.approx(0.8165, abs=1e-3), "sample sd of 0, 1, 2, 1"
    below = th.TARGETS.get("f7-low")
    assert below.below and below.threshold == 0.35
    result = th.compare(below, [0.1, 0.2, 0.5, 0.6], count=None)
    assert result.simulated == 0.5
    above = th.TARGETS.get("f7-tail")
    assert not above.below and th.compare(above, [0.5, 0.9, 1.2, 0.3]).simulated == 0.5
    with pytest.raises(th.TargetError, match="no logarithm"):
        th.compare(log_sd, [1.0, -2.0])
    assert th.subject_of(log_sd).body_class == "giant" and th.subject_of(th.TARGETS.get("a4-mean")) is not None


def test_the_table_refuses_a_malformed_scale_threshold_or_subject():
    base = dict(th.TARGETS.document["targets"][0])
    for key, value, match in (("scale", "ln", "scale"), ("threshold", 0.5, "threshold"), ("below", True, "below"),
                              ("subject", {"quantities": {}}, "subject"), ("subject", "x", "subject")):
        document = {"version": 1, "targets": [dict(base, **{key: value})]}
        with pytest.raises(th.TargetError, match=match):
            th.TargetTable(document)
    proportion = {"id": "p", "row": "t:r", "statistic": "share", "kind": "proportion", "counts": [1, 10], "source": "s",
                  "label": "DATA", "validity": {"q": [0, 1]}}
    th.TargetTable({"version": 1, "targets": [dict(proportion, threshold=0.5, below=True)]})
    with pytest.raises(th.TargetError, match="needs a threshold"):
        th.TargetTable({"version": 1, "targets": [dict(proportion, below=True)]})
    with pytest.raises(th.TargetError, match="names no subject"):
        th.check_target(th.TargetTable({"version": 1, "targets": [proportion]}).get("p"))


# --- radioisotopes (GEN.514) -------------------------------------------------------------------------------------------------


def test_the_thorium_scale_of_non_thin_disc_hosts_is_the_thick_disc_offset_on_the_solar_twin_scale_of_row_o1():
    spec = RECIPES.entry("o2").spec["dist"]
    assert spec["kind"] == "normal" and spec["log10"] and spec["mu"] == 0.061 and spec["sigma"] == 0.19
    assert RECIPES.entry("o1").spec["dist"]["mu"] == 0.08 and 0.08 - 0.019 == pytest.approx(0.061)
    shifted = retuned("o2", lambda e: e["dist"].update(mu=-0.063))
    assert not th.check_target(th.TARGETS.get("o2-mean"), recipes=shifted, n=3000).passed, "the Hypatia scale is 0.12 dex under row o1"


def test_the_thorium_to_uranium_ratio_grows_with_the_planets_age_and_follows_the_published_formula():
    from planetgen.generation.possibility_table import Subject
    dist = lambda age, f: RECIPES._dist_of(RECIPES.entry("o4"), Subject("terrestrial", {"age_gy": age, "th_u_f": f}))   # noqa: E731
    low, high = dist(4.5, 1.16254)["lo"], dist(4.5, 1.71724)["hi"]
    # the report's table (PR 1.751, mass Th/U at 4.5 Gyr) is 3.19 for T = 3 Gyr and 4.72 for T = 13 Gyr; the bounds use PR 1.667 and 2.0
    assert low == pytest.approx(3.19 * 1.667 / 1.751, rel=0.005) and high == pytest.approx(4.72 * 2.0 / 1.751, rel=0.005)
    assert dist(6.0, 1.16254)["lo"] / dist(4.5, 1.16254)["lo"] == pytest.approx(1.1719, rel=0.002), "x1.17 over 1.5 Gyr (decay constants 0.1058 per Gyr)"
    assert dist(0.0, 1.16254)["lo"] == pytest.approx(0.97483 * 1.667 * 1.16254, rel=1e-6)
    assert RECIPES.entry("o4").state == "partly verified" and "disc's star formation history" in RECIPES.entry("o4").limit
    assert "T = 3" not in str(RECIPES.entry("o4").spec["dist"]), "the uniform birth epoch is gone"


def test_the_enriched_aluminium_share_bracket_is_the_published_one():
    constant = RECIPES.entry("o6").constants["c_o6"]
    assert (constant[0], constant[1]) == (0.04, 0.66), "at least 4 percent of systems (Curry 2022) to the white-dwarf upper estimate"


# --- retention and lava worlds (GEN.511) ---------------------------------------------------------------------------------------


def test_the_shoreline_is_the_published_logistic():
    from planetgen.generation.possibility_table import Subject

    def bare(insolation, v_esc_kms=11.186, luminosity=1.0):
        s = Subject("terrestrial", {"insolation": insolation, "v_esc_kms": v_esc_kms, "luminosity_solar": luminosity})
        return dict(RECIPES.probabilities("d", s)).get("d1", 0.0)
    # log10 f0 = 2.62 at the Earth's escape velocity and a solar luminosity (Berta-Thompson 2025): half the planets keep an atmosphere
    assert bare(10 ** 2.62) == pytest.approx(0.5, abs=1e-6)
    assert bare(10 ** (2.62 + 0.247)) == pytest.approx(1.0 / (1.0 + 10 ** -1.0), abs=1e-6), "one width above the line"
    assert bare(10 ** 2.62 * 10 ** 5.9, v_esc_kms=11.186 * 10.0) == pytest.approx(0.5, abs=1e-6), "p = 5.9 in the escape velocity"
    assert bare(10 ** 2.62 * 10 ** 1.17, luminosity=10.0) == pytest.approx(0.5, abs=1e-6), "q = 1.17 in the luminosity"
    assert bare(1.0) < 1e-9, "the Earth's flux is far below the line"


def test_the_retained_share_may_be_zero_and_the_lava_world_share_is_the_beta_interval():
    j1 = RECIPES.groups["j"].constants["c_j1"]
    assert j1[0] == 0.0 and j1[1] == 1.0 and j1[2] == "uniform", "0 of 2 bare M-dwarf planets cannot be forbidden"
    i6 = RECIPES.entry("i6").constants["c_i6"]
    assert (i6[0], i6[1]) == (0.29, 0.99), "the 95 percent interval of Beta(3, 1), 2 volatile-bearing of 2 molten planets"


# --- aerosols and the titanium cold trap (GEN.519) ---------------------------------------------------------------------------


def test_titanium_oxide_is_seen_above_a_drawn_temperature_between_the_last_non_detection_and_the_detection():
    from planetgen.generation.possibility_table import Subject
    from planetgen.util import draw

    def share(teq_k, n=400):
        s = Subject("giant", {"teq_k": teq_k, "zone": "h", "tday_k": 2500.0, "drawn.h": "h3"})
        k = 0
        for i in range(n):
            with draw.bound(draw.Stream(f"h4-test:{teq_k}:{i}")):
                k += RECIPES.draw_binary("h4", s)[0]
        return k / n

    assert RECIPES.entry("h4").constants["c_h4"] == (2350.0, 2641.0, "uniform"), "no TiO at 2350 K, TiO at 2641 K (N=5)"
    assert share(2349.0) == 0.0 and share(2641.0) == 1.0
    assert share(2350.0 + 291.0 / 2) == pytest.approx(0.5, abs=0.08), "a uniform threshold: half the worlds at the middle"
    assert share(2350.0 + 291.0 / 4) == pytest.approx(0.25, abs=0.08)


# --- comets (GEN.520) --------------------------------------------------------------------------------------------------------


def test_the_coma_class_shares_are_the_printed_counts_and_not_the_rounded_words():
    probabilities = dict(RECIPES.probabilities("k", th.subject_of(th.TARGETS.get("k1-share"))))
    assert probabilities == pytest.approx({"k1": 13 / 25, "k2": 9 / 25, "k3": 3 / 25}), "the three printed counts of 25 comets"


def test_the_carbon_dioxide_and_monoxide_ratios_follow_the_comet_database_and_a_wrong_centre_fails():
    assert RECIPES.entry("k4").spec["dist"]["sigma"] == 0.51 and RECIPES.entry("k5").spec["dist"]["sigma"] == 0.65
    assert "range" not in RECIPES.entry("k5").spec["dist"], "the 0.3 to 26 percent range is the sample's extremes, not a bound"
    high = retuned("k5", lambda e: e["dist"].update(mu="log10(0.1)"))
    assert not th.check_target(th.TARGETS.get("k5-mean"), recipes=high, n=3000).passed, "CO/H2O of 10 percent is 0.6 dex above 2.5"
    narrow = retuned("k4", lambda e: e["dist"].update(sigma=0.2))
    assert not th.check_target(th.TARGETS.get("k4-sd"), recipes=narrow, n=3000).passed


def test_the_species_with_a_population_follow_it_and_the_single_comet_species_stay_flagged():
    dists = RECIPES.entry("k7").spec["dists"]
    assert dists["CH4"]["mu"] == "log10(0.0068)" and dists["CH4"]["sigma"] == 0.33
    assert {name: d["sigma"] for name, d in dists.items() if name in ("N2", "O2", "Xe")} == {"N2": "comet_species_sigma()", "O2": "comet_species_sigma()", "Xe": "comet_species_sigma()"}
    assert dists["N2"]["mu"].startswith("log10(comet_species_median('N2', n_over_c, ar_over_c, s_over_c))")
    assert "one comet" in RECIPES.entry("k7").limit and RECIPES.entry("k7").replacement == "GEN.1021"


# --- host-star rotation track (GEN.523) --------------------------------------------------------------------------------------


def test_the_saturation_time_follows_the_mass_and_agrees_with_tu_2015_at_one_solar_mass_and_the_percentile_is_drawn_uniform():
    from planetgen.generation.possibility_table import Subject
    from planetgen.physics import saturation_time

    def t_sat(u, mass=1.0):
        return RECIPES.derive("p6", Subject("terrestrial", {"drawn.p": "p5", "host_mass_solar": mass}), value=u)["t_sat_myr"]
    assert (t_sat(0.1), t_sat(0.5), t_sat(0.9)) == pytest.approx((5.7, 23.0, 226.0), rel=0.1), "Tu 2015: the 10th, 50th and 90th percentile tracks"
    assert t_sat(0.02) < 5.7 < 23.0 < t_sat(0.98), "a monotonic run through the tail percentiles"
    assert t_sat(0.3, 0.3) == pytest.approx(saturation_time.saturation_time_myr(0.3, 30.0)) and t_sat(0.5, 0.3) > 50 * t_sat(0.5, 1.0), "Johnstone 2021: low-mass stars stay saturated longer"
    spec = RECIPES.entry("p6").spec["dist"]
    assert (spec["kind"], spec["lo"], spec["hi"]) == ("uniform", 0.02, 0.98), "no paper weights the tracks"


def test_the_m_dwarf_saturation_time_passes_through_the_three_measured_masses():
    from planetgen.generation.possibility_table import Subject

    def tau(mass):
        return RECIPES.derive("p4", Subject("terrestrial", {"host_type": "M", "host_mass_solar": mass, "age_gy": 1.0}))["tau_sat_gy"]
    assert (tau(0.55), tau(0.25), tau(0.15)) == pytest.approx((0.162, 0.6, 2.3)), "Jackson 2012 and Medina 2022"
    assert tau(0.55) < tau(0.4) < tau(0.25) < tau(0.2) < tau(0.15), "later saturation end for lower masses"


# --- the power of the targets (GEN.665) ---------------------------------------------------------------------------------------

WEAK_TARGETS = {"k5-mean": "the CO/CO2 ratio of 3 comets is a 0.16 dex mean: a factor-2 error in either direction passes",
                "k7-H2S-sd": "half the published spread (0.195 against 0.39, sigma 0.10) passes; the doubled spread is refused"}


def test_a_target_refuses_a_value_wrong_by_a_factor_of_two_unless_it_is_listed_as_weak():
    weak = {}
    for target_id in CHECKED:
        target = th.TARGETS.get(target_id)
        refused = th.power(target, th.check_target(target, n=2000))
        assert refused, target_id
        if not all(refused.values()):
            weak[target_id] = refused
    assert set(weak) == set(WEAK_TARGETS), "a target lost or gained its power: say so in WEAK_TARGETS and the tuning notes"
    assert not any(weak["k5-mean"].values()) and any(weak["k7-H2S-sd"].values())


def test_the_power_of_a_target_is_the_distance_of_the_wrong_value_in_its_combined_sigma():
    target = th.TARGETS.get("f5-z-1p0")
    result = th.check_target(target, n=2000)
    assert th.wrong_values(target) == [11.8, 2.95]
    assert all(th.power(target, result).values())
    proportion = th.TARGETS.get("f7-tail")
    assert th.wrong_values(proportion) == [pytest.approx(0.0869565 + 0.25)], "the lower side would be negative"
