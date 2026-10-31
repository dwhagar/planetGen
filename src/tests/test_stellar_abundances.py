# tests/test_stellar_abundances.py

"""
GEN.439: the stellar composition relations (`planetgen/physics/stellar_abundances.py`, docs/design/composition/stellar-composition.md).

The Sun's composition is a target here and never a preset: it lives in this file and the module takes only the reference
scale. The expected numbers are those the design note computed and checked against its sources (sections 3, 10 and 11).
"""

import math
import random

import pytest

from planetgen import tuning
from planetgen.generation import possibility_table as pt
from planetgen.physics import elements
from planetgen.physics import stellar_abundances as sa

SUN = {"x": 0.7121, "y": 0.2725, "z": 0.0154}
"""dict: The proto-solar X, Y and Z of Asplund et al. 2021 (a test target, section 11)."""

ZERO = {symbol: 0.0 for symbol in sa.METAL_SYMBOLS}
NORMALS_ZERO = [0.0, 0.0, 0.0, 0.0, 0.0]


def mean_composition(fe_h, branch, neon="O"):
    return sa.composition(fe_h, sa.abundance_offsets(fe_h, branch, NORMALS_ZERO, neon_follows=neon))


# --- The scale --------------------------------------------------------------------------------------------------------

def test_the_scale_reproduces_itself():
    star = sa.composition(0.0, ZERO)
    assert (star.x, star.y, star.z) == pytest.approx((SUN["x"], SUN["y"], SUN["z"]), abs=1e-9)
    assert math.fsum(star.elements.values()) == pytest.approx(1.0, abs=1e-12)


def test_the_helium_per_metal_ratio_is_derived_and_inside_the_observed_one():
    assert sa.HELIUM_PER_METAL == pytest.approx(1.674, abs=1e-3)
    assert abs(sa.HELIUM_PER_METAL - 2.1) < 0.9, "Casagrande et al. 2007: 2.1 +- 0.9"
    assert sa.PROTO_Y == pytest.approx(SUN["y"], abs=1e-9)


def test_the_mass_fractions_follow_from_the_log_eps_and_the_atomic_weights():
    part_a = {"C": 2.8356e-3, "N": 7.7519e-4, "O": 6.4144e-3, "Ne": 1.8966e-3, "Mg": 7.0592e-4, "Si": 7.4394e-4,
              "Fe": 1.3184e-3, "Ca": 6.5459e-5}
    for symbol, expected in part_a.items():
        assert sa.PROTO_MASS_FRACTION[symbol] == pytest.approx(expected, rel=2e-3), symbol
    assert sa.UNLISTED_METAL_FRACTION == pytest.approx(4.9e-5, rel=0.05)
    assert sa.UNLISTED_METAL_FRACTION / sa.PROTO_Z < 0.004, "the 14 metals are 99.7 percent of Z"


def test_the_thin_branch_at_solar_metallicity_gives_the_sun_within_three_percent():
    star = mean_composition(0.0, sa.BRANCH_THIN)
    assert star.x == pytest.approx(0.7129, abs=5e-4) and star.y == pytest.approx(0.2720, abs=5e-4)
    assert star.z == pytest.approx(0.01511, abs=2e-4)
    assert abs(star.z - SUN["z"]) / SUN["z"] < 0.03


# Z against 10^[Fe/H] times the solar Z, from the means alone (part A 2.5; GEN.464 reads this table).
Z_OVER_SOLAR = [
    # F, branch, neon follows O, neon follows Mg
    (0.2, "thin", 0.906, 0.897), (0.0, "thin", 0.982, 0.969), (0.0, "high-alpha", 1.298, 1.274),
    (-0.5, "thin", 1.428, 1.382), (-0.5, "high-alpha", 2.043, 1.962), (-1.0, "high-alpha", 2.982, 2.760),
    (-2.0, "high-alpha", 3.027, 2.798)]


@pytest.mark.parametrize("fe_h, branch, with_oxygen, with_magnesium", Z_OVER_SOLAR)
def test_the_metal_fraction_against_the_solar_scaling_matches_the_design_table(fe_h, branch, with_oxygen, with_magnesium):
    assert mean_composition(fe_h, branch, "O").z_over_solar == pytest.approx(with_oxygen, abs=0.01)
    assert mean_composition(fe_h, branch, "Mg").z_over_solar == pytest.approx(with_magnesium, abs=0.01)


def test_the_two_neon_choices_differ_by_at_most_eight_percent_in_z():
    for fe_h in (-2.0, -1.0, -0.5, 0.0):
        a, b = (mean_composition(fe_h, sa.BRANCH_HIGH_ALPHA, neon).z for neon in ("O", "Mg"))
        assert abs(a - b) / max(a, b) <= 0.08 + 1e-9


@pytest.mark.parametrize("branch, low, high", [("thin", -0.9, 0.5), ("high-alpha", -1.6, 0.1)])
def test_the_fractions_add_up_and_the_metals_rise_with_iron_on_each_branch(branch, low, high):
    previous = 0.0
    for step in range(0, 41):
        fe_h = low + (high - low) * step / 40
        star = mean_composition(fe_h, branch)
        assert math.fsum(star.elements.values()) == pytest.approx(1.0, abs=1e-12)
        assert star.x + star.y + star.z == pytest.approx(1.0, abs=1e-12)
        assert star.y >= sa.PRIMORDIAL_Y - 1e-12 and star.x > 0 and star.z > previous, (branch, fe_h)
        previous = star.z


def test_the_scale_keeps_helium_from_primordial_up_to_the_helium_per_metal_line():
    star = mean_composition(-0.4, sa.BRANCH_THIN)
    assert star.y == pytest.approx(sa.PRIMORDIAL_Y + sa.HELIUM_PER_METAL * star.z, abs=1e-9)
    metal_free = sa.composition(-30.0, ZERO)
    assert metal_free.z < 1e-20 and metal_free.x == pytest.approx(1.0 - sa.PRIMORDIAL_Y, abs=1e-12)


# --- The worked examples of the design note ----------------------------------------------------------------------------

def test_example_a_a_thin_main_sequence_star():
    fe_h = -0.20
    offsets = sa.abundance_offsets(fe_h, "thin", [-0.15, 1.30, 0.005, 0.42, 0.0])
    assert (offsets["C"], offsets["O"], offsets["Mg"], offsets["Si"]) == pytest.approx((-0.012, 0.194, 0.003, 0.017), abs=1e-3)
    assert offsets["Ca"] == pytest.approx(0.072, abs=1e-3) and offsets["Ne"] == offsets["O"]
    star = sa.composition(fe_h, offsets)
    assert (star.x, star.y, star.z) == pytest.approx((0.7190, 0.2682, 0.01282), abs=1e-4)
    assert star.z_over_solar == pytest.approx(1.319, abs=2e-3)
    assert star.elements["O"] == pytest.approx(6.38e-3, rel=5e-3) and star.elements["C"] == pytest.approx(1.76e-3, rel=1e-2)
    assert sa.mass_ratio(star, "Fe", "Si") == pytest.approx(1.703, abs=2e-3)
    assert sa.mass_ratio(star, "Mg", "Si") == pytest.approx(0.917, abs=2e-3)
    assert sa.mass_ratio(star, "C", "O") == pytest.approx(0.2756, abs=2e-3)


def test_the_ratios_the_planets_inherit_at_the_zero_point_are_the_scales():
    star = sa.composition(0.0, ZERO)
    assert sa.mass_ratio(star, "Fe", "Si") == pytest.approx(1.772, abs=3e-3)
    assert sa.mass_ratio(star, "Mg", "Si") == pytest.approx(0.949, abs=3e-3)
    assert sa.mass_ratio(star, "C", "O") == pytest.approx(0.442, abs=3e-3)
    assert sa.number_ratio(star, "Fe", "Si") == pytest.approx(0.891, abs=1e-3)
    assert sa.number_ratio(ZERO, "Fe", "Si") == sa.number_ratio(star, "Fe", "Si"), "the offsets dict gives the same"


def test_the_ratios_do_not_depend_on_the_overall_metallicity():
    one = sa.composition(-1.0, ZERO)
    other = sa.composition(0.3, ZERO)
    assert sa.mass_ratio(one, "Mg", "Si") == sa.mass_ratio(other, "Mg", "Si")


def test_high_alpha_stars_have_less_iron_per_silicon_and_more_magnesium():
    thin = mean_composition(-0.5, "thin")
    high = mean_composition(-0.5, "high-alpha")
    assert sa.mass_ratio(high, "Fe", "Si") < sa.mass_ratio(thin, "Fe", "Si")
    assert sa.mass_ratio(high, "Mg", "Si") > sa.mass_ratio(thin, "Mg", "Si")
    assert sa.mass_ratio(high, "C", "O") < sa.mass_ratio(thin, "C", "O")


# --- The branch --------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("fe_h, expected", [(-0.35, 0.5), (-0.5, 0.757), (-0.8, 0.968), (-1.6, 0.9999)])
def test_the_high_alpha_probability(fe_h, expected):
    assert sa.high_alpha_probability(fe_h) == pytest.approx(expected, abs=1e-3)


def test_no_star_younger_than_the_age_gate_is_high_alpha():
    assert tuning.STELLAR_THICK_START_GY == 8.0
    for fe_h in (-2.0, -0.5, 0.0):
        assert sa.branch_probability(fe_h, 7.999) == 0.0
        assert sa.branch_from_uniform(fe_h, 7.999, 0.0) == sa.BRANCH_THIN
        assert sa.branch_probability(fe_h, 8.0) == sa.high_alpha_probability(fe_h)


def test_the_branch_flips_at_most_once_as_metallicity_falls_for_a_fixed_number():
    for u in (0.05, 0.3, 0.5, 0.9, 0.999):
        flips, previous = 0, None
        for step in range(0, 301):
            fe_h = 0.6 - step * 0.01
            branch = sa.branch_from_uniform(fe_h, 11.0, u)
            if previous is not None and branch != previous:
                flips += 1
            previous = branch
        assert flips <= 1, u


def test_the_branch_table_agrees_with_the_branch_probability():
    for fe_h, age in ((-0.5, 9.0), (-1.6, 12.0), (0.2, 12.0), (-0.5, 2.0), (-0.2, 7.5)):
        found = sa.branch_probabilities(fe_h, age)
        assert math.fsum(found.values()) == pytest.approx(1.0, abs=1e-12)
        assert found["high-alpha"] == pytest.approx(sa.branch_probability(fe_h, age), abs=1e-12)
    table = sa.BRANCH_TABLE
    assert table in pt.REGISTRY.tables() and [row.name for row in table.rows][0] == "high-alpha sequence"
    assert table.state_counts()[pt.PARTLY_VERIFIED] == 2 and table.state_counts()[pt.OPEN_STATE] == 0
    subject = pt.Subject("star", {"fe_h": -0.5, "age_gy": 9.0})
    assert [name for name, _p in table.draw(subject).probabilities] == ["high-alpha sequence", "thin sequence"]


# --- The offsets -------------------------------------------------------------------------------------------------------

def test_the_means_hold_at_the_edge_of_the_fitted_range():
    assert sa.mean_xfe("Mg", 1.0, "thin") == sa.mean_xfe("Mg", 0.5, "thin")
    assert sa.mean_xfe("O", -3.0, "high-alpha") == sa.mean_xfe("O", -1.6, "high-alpha")
    assert sa.mean_xfe("O", -1.6, "high-alpha") == pytest.approx(0.606 - 0.495 * 0, abs=1e-9), "the plateau below the knee"
    with pytest.raises(sa.StellarAbundanceError):
        sa.mean_xfe("Zr", 0.0, "thin")
    with pytest.raises(sa.StellarAbundanceError):
        sa.mean_xfe("Mg", 0.0, "thick")


def test_the_residuals_have_the_stated_scatter_and_correlation():
    generator = random.Random(439)
    draws = [sa.residual_offsets("thin", [generator.gauss(0, 1) for _ in range(5)]) for _ in range(40000)]
    sd = {symbol: math.sqrt(sum(d[symbol] ** 2 for d in draws) / len(draws)) for symbol in sa.RESIDUAL_ELEMENTS}
    for symbol, expected in tuning.STELLAR_XFE_SIGMA["thin"].items():
        if symbol in sd:
            assert sd[symbol] == pytest.approx(expected, rel=0.03), symbol
    correlation = sum(d["Mg"] * d["Si"] for d in draws) / len(draws) / (sd["Mg"] * sd["Si"])
    assert correlation == pytest.approx(0.68, abs=0.03)
    ratio_sd = math.sqrt(sum((d["Mg"] - d["Si"]) ** 2 for d in draws) / len(draws))
    assert ratio_sd == pytest.approx(0.054, rel=0.20), "the Mg/Si ratio scatter is narrower than the sum of either element's"
    with pytest.raises(sa.StellarAbundanceError):
        sa.residual_offsets("thin", [0.0, 0.0])


def test_the_residual_transformation_is_the_cholesky_factor_of_the_correlation():
    for i, row in enumerate(tuning.STELLAR_XFE_CORRELATION):
        for j, value in enumerate(row):
            assert sum(sa._CHOLESKY[i][k] * sa._CHOLESKY[j][k] for k in range(5)) == pytest.approx(value, abs=1e-12)


def test_carbon_enhanced_metal_poor_stars():
    assert sa.carbon_enhanced(-1.9, 0.0, 0.5) is None, "only below [Fe/H] -2"
    assert sa.carbon_enhanced(-3.0, 0.30, 0.5) is None, "a share of 0.25"
    assert sa.carbon_enhanced(-3.0, 0.069, 0.396) == pytest.approx(0.7 + 1.3 * 0.396)
    assert sa.carbon_enhanced(-4.2, 0.99, 0.0) == pytest.approx(0.7), "most stars below -4 (out of the stored range)"
    offsets = sa.abundance_offsets(-3.0, "high-alpha", NORMALS_ZERO, u_enhanced=0.069, u_carbon=0.396)
    assert offsets["C"] == pytest.approx(1.2148, abs=1e-3)
    star = sa.composition(-3.0, offsets)
    assert sa.mass_ratio(star, "C", "O") > 1.0, "a carbon-rich rocky planet is possible there"
    plain = sa.abundance_offsets(-3.0, "high-alpha", NORMALS_ZERO)
    assert plain["C"] == pytest.approx(sa.mean_xfe("C", -3.0, "high-alpha"))


def test_the_other_elements_follow_the_stated_ones():
    offsets = sa.abundance_offsets(-0.3, "thin", [0.3, -0.2, 0.1, 0.4, 0.0])
    assert offsets["S"] == offsets["Ar"] == offsets["Ca"] != sa.mean_xfe("Ca", -0.3, "thin"), "they carry calcium's residual (GEN.988)"
    assert offsets["Ne"] == offsets["O"]
    assert sa.abundance_offsets(-0.3, "thin", [0.3, -0.2, 0.1, 0.4, 0.0], neon_follows="Mg")["Ne"] == offsets["Mg"]
    assert all(offsets[symbol] == 0.0 for symbol in ("N", "Na", "Al", "K", "Fe", "Ni"))
    with pytest.raises(sa.StellarAbundanceError):
        sa.abundance_offsets(0.0, "thin", NORMALS_ZERO, neon_follows="C")


def test_the_same_inputs_give_the_same_composition():
    one = sa.composition(-0.7, sa.abundance_offsets(-0.7, "high-alpha", [0.1, 0.2, -0.3, 0.4, 0.0]))
    two = sa.composition(-0.7, sa.abundance_offsets(-0.7, "high-alpha", [0.1, 0.2, -0.3, 0.4, 0.0]))
    assert one == two


def test_a_random_sweep_of_stars_always_adds_up():
    generator = random.Random(7)
    for _ in range(300):
        fe_h = generator.uniform(-3.5, 0.6)
        age = generator.uniform(0.0, 13.0)
        branch = sa.branch_from_uniform(fe_h, age, generator.random())
        offsets = sa.abundance_offsets(fe_h, branch, [generator.gauss(0, 1) for _ in range(5)], generator.random(),
                                       generator.random())
        star = sa.composition(fe_h, offsets)
        assert math.fsum(star.elements.values()) == pytest.approx(1.0, abs=1e-12)
        assert all(value >= 0 for value in star.elements.values()) and star.y >= sa.PRIMORDIAL_Y - 1e-12


def test_the_solar_neighbourhood_mix_reproduces_the_research_ranges():
    generator = random.Random(2026)
    fe_si, mg_si, c_o, high = [], [], [], 0
    for _ in range(6000):
        fe_h = min(max(generator.gauss(-0.10, 0.22), -3.5), 0.6)
        branch = sa.branch_from_uniform(fe_h, 11.0, generator.random()) if generator.random() < 1.0 else "thin"
        high += branch == sa.BRANCH_HIGH_ALPHA
        star = sa.composition(fe_h, sa.abundance_offsets(fe_h, branch, [generator.gauss(0, 1) for _ in range(5)]))
        fe_si.append(sa.mass_ratio(star, "Fe", "Si"))
        mg_si.append(sa.mass_ratio(star, "Mg", "Si"))
        c_o.append(sa.mass_ratio(star, "C", "O"))
    median = lambda values: sorted(values)[len(values) // 2]  # noqa: E731
    # The medians of Brewer and Fischer 2016 (Fe/Si 1.81, Mg/Si 0.93, C/O 0.385), which the thin Mg and Si were refitted to (GEN.989).
    assert median(fe_si) == pytest.approx(1.81, abs=0.2) and median(mg_si) == pytest.approx(0.93, abs=0.1)
    assert median(c_o) == pytest.approx(0.385, abs=0.05)
    assert high / 6000 == pytest.approx(0.216, abs=0.03), "the high-alpha share with the age gate off, as part A drew it"


# --- Dredge-up and the radioactive birth abundances --------------------------------------------------------------------

def test_the_first_dredge_up_of_a_one_and_a_half_solar_mass_star():
    dredge = sa.first_dredge_up(1.5, 0.01257)
    assert dredge.weight == pytest.approx(0.845, abs=1e-3) and dredge.delta_y == pytest.approx(0.0165, abs=2e-4)
    assert dredge.carbon_factor == pytest.approx(0.7257, abs=1e-3) and dredge.nitrogen_factor == pytest.approx(2.02, abs=5e-3)
    assert not dredge.mass_held and not dredge.z_held


def test_the_first_dredge_up_at_the_tabulated_corner():
    dredge = sa.first_dredge_up(1.0, 0.02)
    assert dredge.delta_y == pytest.approx(0.3025 - 0.28, abs=1e-9)
    assert dredge.carbon_factor == pytest.approx(0.1530 / 0.1733, abs=1e-9)
    assert dredge.nitrogen_factor == pytest.approx(0.0771 / 0.0531, abs=1e-9)
    # The design note's check: 1 solar mass at Z 0.02 gives +0.0225 helium, C/Z x0.88, N/Z x1.45.
    assert (dredge.delta_y, dredge.carbon_factor, dredge.nitrogen_factor) == pytest.approx((0.0225, 0.88, 1.45), rel=0.1)


def test_the_dredge_up_holds_the_table_edge_outside_its_range_and_says_so():
    assert sa.first_dredge_up(0.8, 0.02).mass_held and sa.first_dredge_up(3.0, 0.02).mass_held
    assert sa.first_dredge_up(1.5, 0.05).z_held and sa.first_dredge_up(1.5, 1e-5).z_held
    assert sa.first_dredge_up(0.8, 0.02).delta_y == sa.first_dredge_up(1.0, 0.02).delta_y
    assert sa.first_dredge_up(2.5, 0.001).delta_y == sa.first_dredge_up(2.0, 0.001).delta_y, "no 2.5 row at low Z"


def test_example_f_the_red_giant_surface():
    fe_h = -0.10
    birth = sa.composition(fe_h, sa.abundance_offsets(fe_h, "thin", NORMALS_ZERO))
    giant = sa.after_first_dredge_up(birth, 1.5)
    dredge = sa.first_dredge_up(1.5, birth.z)
    assert giant.y == pytest.approx(birth.y + dredge.delta_y)
    assert giant.z == pytest.approx(birth.z, rel=0.01), "the net change of C and N is 0.2 percent of Z"
    assert giant.x + giant.y + giant.z == pytest.approx(1.0, abs=1e-12)
    assert math.fsum(giant.elements.values()) == pytest.approx(1.0, abs=1e-12)
    assert giant.x == pytest.approx(0.7032, abs=1e-3) and giant.y == pytest.approx(0.2842, abs=1e-3)
    assert giant.elements["C"] / birth.elements["C"] == pytest.approx(dredge.carbon_factor)
    assert giant.elements["N"] / birth.elements["N"] == pytest.approx(dredge.nitrogen_factor)
    assert giant.elements["Mg"] == birth.elements["Mg"]


def test_the_birth_thorium_and_uranium_of_the_solar_scale():
    assert sa.solar_birth_th_over_eu() == pytest.approx(0.405, rel=0.01)
    assert sa.solar_birth_u_over_th() == pytest.approx(0.572, rel=0.015)


def test_the_radioactive_birth_abundances_follow_europium():
    star = sa.composition(0.0, ZERO)
    found = sa.birth_radiogenic(star, "thin")
    eu, th, uranium = (found["abundances"][symbol] for symbol in ("Eu", "Th", "U"))
    assert th / eu == pytest.approx(sa.solar_birth_th_over_eu() * 10.0 ** 0.014)
    assert uranium / th == pytest.approx(sa.solar_birth_u_over_th())
    assert found["eu_fe"] == pytest.approx(sa.mean_xfe("Eu", 0.0, "thin"))
    heavy = sa.birth_radiogenic(star, "thin", eu_normal=1.0, th_eu_normal=1.0)
    assert heavy["abundances"]["Th"] > th and heavy["th_eu"] == pytest.approx(0.014 + 0.045)
    metal_poor = sa.birth_radiogenic(sa.composition(-1.0, ZERO), "high-alpha")
    assert metal_poor["abundances"]["Eu"] < eu, "less europium at lower iron"
    atoms_per_kg = found["mass_fractions"]["Th"] / elements.element("Th").molar_mass
    assert atoms_per_kg > 0 and found["mass_fractions"]["U"] < found["mass_fractions"]["Th"] * 2


# --- Remnants ----------------------------------------------------------------------------------------------------------

def test_the_white_dwarf_core_by_mass():
    assert sa.white_dwarf_core(0.3) == "helium" and sa.white_dwarf_core(0.6) == "carbon-oxygen"
    assert sa.white_dwarf_core(1.2) == "oxygen-neon" and sa.white_dwarf_core(1.05) == "carbon-oxygen"
    assert tuning.STELLAR_WD_DA_SHARE == 0.80 and tuning.STELLAR_WD_POLLUTED_SHARE == (0.25, 0.50)


def test_the_neutron_star_constants():
    low, high = tuning.STELLAR_NS_RADIUS_BAND_KM
    assert (low, high) == (11.6, 13.1) and low < 12.45 < high, "Miller 2021's R(1.4)"
    assert tuning.STELLAR_NS_MASS_FLOOR_SOL < tuning.STELLAR_NS_MAXIMUM_NONROTATING_MASS_SOL[0]
    assert (tuning.STELLAR_NS_MASS_HEAVIEST_SOL, tuning.STELLAR_NS_MASS_HEAVIEST_SIGMA_SOL) == (2.35, 0.17)
    # Mean density of example K: 1.36 solar masses within 12.45 km is 1.22 times saturation.
    mean_density = 3 * 1.36 * 1.98841e33 / (4 * math.pi * (12.45e5) ** 3)
    assert mean_density / sa.NUCLEAR_SATURATION_DENSITY_G_CM3[0] == pytest.approx(1.22, rel=0.03)


def test_black_hole_quantities():
    hole = sa.black_hole(10.0)
    assert hole.schwarzschild_radius_km == pytest.approx(29.53, rel=1e-3)
    assert hole.hawking_temperature_k == pytest.approx(6.17e-9, rel=2e-3)
    assert hole.horizon_radius_km == pytest.approx(hole.schwarzschild_radius_km) and hole.isco_gm_c2 == pytest.approx(6.0)
    assert hole.efficiency == pytest.approx(0.057, abs=1e-3)
    spinning = sa.black_hole(10.0, 0.9)
    assert spinning.isco_gm_c2 == pytest.approx(2.321, abs=2e-3) and spinning.efficiency == pytest.approx(0.156, abs=1e-3)
    example = sa.black_hole(12.0, 0.7)
    assert example.schwarzschild_radius_km == pytest.approx(35.44, abs=0.01) and example.horizon_radius_km == pytest.approx(30.37, abs=0.01)
    assert example.isco_gm_c2 == pytest.approx(3.393, abs=1e-3) and example.efficiency == pytest.approx(0.1036, abs=5e-4)
    with pytest.raises(sa.StellarAbundanceError):
        sa.black_hole(10.0, 1.0)


def test_brown_dwarf_limits_and_radius():
    assert sa.deuterium_burning_limit_mj(0.0189, 0.25) == pytest.approx(13.55)
    assert sa.deuterium_burning_limit_mj(0.0154) == pytest.approx(13.2, abs=0.3)
    assert sa.deuterium_burning_limit_mj(0.0154) > sa.deuterium_burning_limit_mj(0.05), "metal-rich burns deuterium sooner"
    assert 11.0 <= sa.deuterium_burning_limit_mj(1e-8) <= 16.3 and 11.0 <= sa.deuterium_burning_limit_mj(0.2) <= 16.3
    assert sa.lithium_burning_limit_sol() * 1.98841e30 / 1.89813e27 == pytest.approx(63, abs=1.0)
    assert sa.hydrogen_burning_limit_sol(0.0) == pytest.approx(0.092) and sa.hydrogen_burning_limit_sol(0.0154) == pytest.approx(0.074)
    assert sa.hydrogen_burning_limit_sol(0.0154) < sa.hydrogen_burning_limit_sol(0.005) < sa.hydrogen_burning_limit_sol(0.0)
    assert sa.brown_dwarf_radius_solar(0.07) == pytest.approx(0.0944, rel=0.05)           # at 1 Gyr, the default
    assert sa.brown_dwarf_radius_solar(0.04, 0.1) == pytest.approx(0.13)                   # GEN.948: any age of the table
    assert sa.brown_dwarf_radius_solar(0.07, 0.1) > sa.brown_dwarf_radius_solar(0.07, 1.0) >= sa.brown_dwarf_radius_solar(0.07, 10.0)
    for mass in (0.005, 2.0):
        with pytest.raises(sa.StellarAbundanceError):
            sa.brown_dwarf_radius_solar(mass)
    # A fifty-two Jupiter-mass object (example I) burned its deuterium, kept its lithium and burns no hydrogen.
    mass = 0.050
    assert mass * 1.98841e30 / 1.89813e27 > sa.deuterium_burning_limit_mj(0.0154)
    assert mass < sa.lithium_burning_limit_sol() < sa.hydrogen_burning_limit_sol(0.0154)


# --- The relations list ------------------------------------------------------------------------------------------------

def test_every_relation_names_a_source_and_a_state_and_a_limit_when_it_is_not_verified():
    assert sa.RELATIONS and len({relation.name for relation in sa.RELATIONS}) == len(sa.RELATIONS)
    for relation in sa.RELATIONS:
        assert relation.source.strip() and relation.state in pt.STATES, relation.name
        if relation.state != sa.VERIFIED:
            assert relation.limit.strip(), f"{relation.name} is {relation.state} and names no limit"
    states = {relation.state for relation in sa.RELATIONS}
    assert states == set(pt.STATES), "the list shows each of the four states"


def test_no_generator_table_holds_the_suns_composition_as_a_default():
    # The scale's zero point is a reference, so no star kind, branch or fit defaults to the Sun's own numbers.
    assert not hasattr(sa, "SUN") and not hasattr(sa, "SOLAR_COMPOSITION")
    for fits in tuning.STELLAR_XFE_FITS.values():
        assert "Fe" not in fits, "iron is the reference of [X/Fe]"


def test_alpha_enhancement_is_the_mean_of_the_four_alpha_offsets_and_the_high_alpha_branch_has_more():
    offsets = {symbol: 0.0 for symbol in sa.METAL_SYMBOLS}
    offsets.update(O=0.1, Mg=0.3, Si=0.2, Ca=0.0, Fe=0.5)    # iron and the others do not count
    assert sa.alpha_enhancement(offsets) == pytest.approx(0.15)
    thin = sa.composition(-0.8, sa.abundance_offsets(-0.8, sa.BRANCH_THIN, (0.0,) * 5))
    high = sa.composition(-0.8, sa.abundance_offsets(-0.8, sa.BRANCH_HIGH_ALPHA, (0.0,) * 5))
    assert sa.alpha_enhancement(thin) == sa.alpha_enhancement(thin.offsets)
    assert 0.0 < sa.alpha_enhancement(thin) < 0.25 < sa.alpha_enhancement(high) < 0.45
