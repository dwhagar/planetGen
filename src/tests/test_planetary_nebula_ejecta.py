"""GEN.1270: the progenitor and ejecta of a planetary nebula (`planetgen/physics/planetary_nebula_ejecta.py`): the progenitor to ejecta
relation, the class windows, determinism and the conservation of mass."""

import math

import pytest

from planetgen.physics import agb_carbon, planetary_nebula_ejecta as pn, remnant_mass


def keys(n):
    return [f"nebula-{i}" for i in range(n)]


# --- the progenitor to ejecta relation -----------------------------------------------------------------------------------

@pytest.mark.parametrize("mass", [1.0, 1.5, 2.0, 3.0, 3.9, 4.0, 5.0, 6.5, 7.0])
@pytest.mark.parametrize("fe_h", [0.0, -0.5, 0.2])
def test_the_ejecta_mass_fractions_add_to_one_and_stay_positive(mass, fe_h):
    e = pn.ejecta_of(pn.H, mass, fe_h)
    assert sum(e.composition.elements.values()) == pytest.approx(1.0, abs=1e-9)
    assert min(e.composition.elements.values()) >= 0.0
    assert e.composition.x + e.composition.y + e.composition.z == pytest.approx(1.0, abs=1e-9)


def test_the_ejected_mass_is_the_initial_mass_less_the_white_dwarf():
    for mass in (1.0, 2.5, 5.0):
        e = pn.ejecta_of(pn.H, mass, 0.0)
        assert e.white_dwarf_mass_sol == pytest.approx(remnant_mass.white_dwarf_mass_sol(mass))
        assert e.ejected_mass_sol == pytest.approx(mass - e.white_dwarf_mass_sol)
        assert 0.5 < e.white_dwarf_mass_sol < 1.35 and e.ejected_mass_sol > 0.3


def test_a_star_too_light_to_dredge_carbon_up_stays_oxygen_rich_and_a_heavier_one_turns_carbon_rich():
    light = pn.ejecta_of(pn.H, 1.0, 0.0)
    assert not light.third_dredge_up and light.pulses == 0 and light.dredged_mass_sol == 0.0
    assert light.c_to_o < 1.0 and not light.carbon_rich
    for mass in (1.6, 2.0, 3.0):
        e = pn.ejecta_of(pn.H, mass, 0.0)
        assert e.third_dredge_up and e.pulses > 10 and e.dredged_mass_sol > 0.01
        assert 1.0 < e.c_to_o < 4.0 and e.carbon_rich and not e.hot_bottom_burning


def test_the_third_dredge_up_mixes_the_intershell_into_the_lost_envelope_by_mass():
    """Each fraction is (1 - share) of the envelope after the first dredge-up plus share of the intershell, share the dredged mass over the
    mass the thermally pulsing branch loses."""
    mass = 2.0
    e = pn.ejecta_of(pn.H, mass, 0.0)
    fe_h = e.birth_fe_h
    envelope = pn.sa.after_first_dredge_up(pn.sa.composition(fe_h, pn.sa.abundance_offsets(fe_h, pn.sa.BRANCH_THIN, (0.0,) * 5)), mass).elements
    pulses = int(agb_carbon.pulses_in_branch(mass))
    assert e.pulses == pulses and e.dredged_mass_sol == pytest.approx(agb_carbon.dredged_core_mass_sol(mass, pulses))
    share = e.dredged_mass_sol / (mass - pn.RGB_MASS_LOSS_SOL - e.white_dwarf_mass_sol)
    assert e.composition.elements["C"] == pytest.approx((1 - share) * envelope["C"] + share * agb_carbon.INTERSHELL_CARBON_FRACTION)
    assert e.composition.elements["He"] == pytest.approx((1 - share) * envelope["He"] + share * pn.INTERSHELL_HE)
    assert e.composition.elements["O"] == pytest.approx((1 - share) * envelope["O"] + share * pn.INTERSHELL_O)
    assert e.composition.elements["N"] == pytest.approx((1 - share) * envelope["N"])
    assert e.composition.elements["H"] == pytest.approx((1 - share) * envelope["H"])


def test_more_dredge_up_raises_carbon_at_lower_metallicity_for_the_same_mass():
    rich, poor = pn.ejecta_of(pn.H, 2.5, 0.0), pn.ejecta_of(pn.H, 2.5, -1.0)
    assert poor.c_to_o > rich.c_to_o                         # the same carbon over less oxygen


def test_hot_bottom_burning_turns_a_heavy_progenitors_carbon_into_nitrogen_and_makes_a_type_i_nebula():
    below, above = pn.ejecta_of(pn.H, 3.9, 0.0), pn.ejecta_of(pn.H, 4.1, 0.0)
    assert not below.hot_bottom_burning and above.hot_bottom_burning
    assert above.n_to_o > 4.0 * below.n_to_o and above.c_to_o < below.c_to_o
    assert above.composition.y - below.composition.y < 0.01              # the second dredge-up has hardly begun at 4.1 solar masses
    assert above.type_i and not below.type_i
    solar_five = pn.ejecta_of(pn.L, 5.0, 0.0)
    assert math.log10(solar_five.n_to_o) > pn.TYPE_I_LOG_N_TO_O and solar_five.carbon_rich is False


@pytest.mark.parametrize("mass, delta_y", [(1.0, 0.0), (4.0, 0.0), (4.5, 0.004), (5.0, 0.018), (5.5, 0.033), (6.0, 0.044), (7.0, 0.059), (8.0, 0.066)])
def test_the_second_dredge_up_adds_the_helium_of_the_table_of_karakas_and_lattanzio(mass, delta_y):
    """[R5098: Table 1], Z = 0.02: the rise in Y from after the first to after the second dredge-up, at the rows of the table."""
    assert pn.second_dredge_up_delta_y(mass) == pytest.approx(delta_y)


def test_the_second_dredge_up_is_linear_between_rows_and_held_above_the_last():
    assert pn.second_dredge_up_delta_y(4.75) == pytest.approx(0.011)
    assert pn.second_dredge_up_delta_y(6.5) == pytest.approx(0.0515)
    assert pn.second_dredge_up_delta_y(12.0) == pn.second_dredge_up_delta_y(8.0) == 0.066
    values = [pn.second_dredge_up_delta_y(4.0 + 0.1 * i) for i in range(50)]
    assert all(b >= a for a, b in zip(values, values[1:]))


def test_the_second_dredge_up_raises_the_helium_of_a_heavy_progenitor_by_its_row(monkeypatch):
    with_it = pn.ejecta_of(pn.L, 7.0, 0.0).composition.y
    monkeypatch.setattr(pn, "SECOND_DREDGE_UP_TABLE", ((4.0, 0.0), (8.0, 0.0)))
    assert with_it - pn.ejecta_of(pn.L, 7.0, 0.0).composition.y == pytest.approx(0.059, abs=0.001)


def test_neon_sulfur_argon_and_iron_follow_the_birth_composition():
    e = pn.ejecta_of(pn.H, 1.0, 0.0)             # no dredge-up reaches the surface of a one solar mass star in this model beyond the first
    for symbol in ("S", "Ar", "Fe", "Mg", "Si"):
        assert e.composition.elements[symbol] == pytest.approx(
            pn.sa.composition(e.birth_fe_h, pn.sa.abundance_offsets(e.birth_fe_h, pn.sa.BRANCH_THIN, (0.0,) * 5)).elements[symbol])


def test_an_older_progenitor_formed_from_poorer_gas():
    assert pn.birth_fe_h(0.0, 1.0) < pn.birth_fe_h(0.0, 3.0) < pn.birth_fe_h(0.0, 6.0) < 0.0
    assert pn.birth_fe_h(0.0, 1.0) == pytest.approx(0.0 + pn.age_metallicity.THIN_SLOPE_DEX_PER_GY * pn.lifespan_gy(1.0, 0.0))


# --- the progenitor draw -----------------------------------------------------------------------------------------------

def test_the_same_nebula_always_has_the_same_progenitor_and_another_nebula_another():
    a = [pn.progenitor_mass_sol(pn.H, 0.0, k) for k in keys(20)]
    assert a == [pn.progenitor_mass_sol(pn.H, 0.0, k) for k in keys(20)]
    assert len(set(a)) > 15
    assert pn.progenitor("K", 0.0, "x") == pn.progenitor("K", 0.0, "x")
    assert pn.progenitor_mass_sol(pn.H, 0.0, "x") != pn.progenitor_mass_sol(pn.H, 0.0, "y")


@pytest.mark.parametrize("nebula_class", pn.CLASSES)
def test_every_progenitor_lies_in_its_classs_window(nebula_class):
    low, high = pn.mass_window_sol(nebula_class, 0.0)
    for key in keys(200):
        mass = pn.progenitor_mass_sol(nebula_class, 0.0, key)
        assert low - 1e-9 <= mass <= high + 1e-9


def test_the_windows():
    assert pn.mass_window_sol(pn.H, 0.0) == pn.mass_window_sol(pn.J, 0.0) == (pn.MIN_MASS_SOL, pn.MAX_MASS_SOL)
    assert pn.mass_window_sol(pn.K, 0.0) == (agb_carbon.tdu_minimum_mass_sol(0.0), agb_carbon.HBB_MIN_MASS_SOL)
    assert pn.mass_window_sol(pn.L, 0.0) == (agb_carbon.HBB_MIN_MASS_SOL, pn.MAX_MASS_SOL)
    assert pn.mass_window_sol(pn.K, -1.0)[0] < pn.mass_window_sol(pn.K, 0.0)[0]        # a poorer star dredges at a lower mass
    with pytest.raises(ValueError):
        pn.mass_window_sol("M", 0.0)


def test_a_class_k_nebula_is_always_carbon_rich_and_a_class_l_always_type_i():
    for key in keys(150):
        k, l = pn.progenitor(pn.K, 0.0, key), pn.progenitor(pn.L, 0.0, key)
        assert k.carbon_rich and k.c_to_o >= 1.0 and not k.hot_bottom_burning
        assert l.type_i and l.hot_bottom_burning and not l.carbon_rich


def test_the_progenitors_follow_the_initial_mass_function():
    """Light stars are far more numerous: the median of a class H or J draw is under 2 solar masses and few are over 4."""
    masses = sorted(pn.progenitor_mass_sol(pn.H, 0.0, k) for k in keys(600))
    assert masses[len(masses) // 2] < 2.0
    assert sum(m > 4.0 for m in masses) / len(masses) < 0.2
    assert masses[0] >= pn.MIN_MASS_SOL and masses[-1] <= pn.MAX_MASS_SOL


def test_a_later_central_star_reads_the_same_mass_from_the_same_key():
    """The progenitor mass is a function of the class, the place's [Fe/H] and the key alone."""
    first = pn.progenitor(pn.J, -0.2, "FE81000A2B-0000005-000")
    again = pn.progenitor_mass_sol(pn.J, -0.2, "FE81000A2B-0000005-000")
    assert again == first.progenitor_mass_sol and first.white_dwarf_mass_sol == pytest.approx(remnant_mass.white_dwarf_mass_sol(again))


def test_the_values_to_store():
    values = pn.progenitor(pn.K, 0.0, "x").values()
    assert {"progenitor_mass_sol", "white_dwarf_mass_sol", "mass_fraction_he", "mass_fraction_c", "mass_fraction_n", "mass_fraction_o", "c_to_o",
            "type_i"} <= set(values)
    assert values["mass_fraction_h"] + values["mass_fraction_he"] + values["mass_fraction_z"] == pytest.approx(1.0)
