# tests/test_nebula_composition.py

"""GEN.1269: the shared composition model of nebulae and remnants (`planetgen/physics/nebula_composition.py`): mass is conserved, a metal-poor place
gives a metal-poor nebula, every class has a row, the dust law reproduces its source's points, and the phases and dust follow the class."""

import math

import pytest

from planetgen import tuning
from planetgen.physics import nebula_composition as nc
from planetgen.physics import planetary_nebula_ejecta as pne

CLASS_LETTERS = sorted(tuning.NEBULA_CLASSES)
COLD = ("F", "M", "N", "P", "Q")


def make(letter, oh=0.0, density=300.0, radius_ly=10.0, **kwargs):
    return nc.makeup(letter, oh, density, radius_ly, **kwargs)


# --- the table ---------------------------------------------------------------------------------------------------------------------------

def test_every_class_of_the_table_has_a_row_a_phase_split_and_a_known_grain_mix():
    assert set(nc.CLASSES) == set(tuning.NEBULA_CLASSES) == set(nc.PHASES)
    for letter, row in nc.CLASSES.items():
        assert row.mix in nc.GRAIN_MIXES, letter
        assert set(row.survival) <= set(nc.GRAINS), letter
        assert all(0.0 <= value <= 1.0 for value in row.survival.values()), letter


def test_every_grain_mix_adds_to_one_and_names_known_grains():
    for name, mix in nc.GRAIN_MIXES.items():
        assert sum(mix.values()) == pytest.approx(1.0), name
        assert set(mix) <= set(nc.GRAINS), name
    for grain, shares in nc.GRAINS.items():
        assert sum(shares.values()) == pytest.approx(1.0), grain


@pytest.mark.parametrize("letter", CLASS_LETTERS)
def test_a_nebula_of_any_class_keeps_its_mass_and_its_shares(letter):
    m = make(letter, density=1e3 if letter in COLD else 100.0)
    assert sum(m.elements.values()) == pytest.approx(1.0, abs=1e-12)
    assert m.x + m.y + m.z == pytest.approx(1.0, abs=1e-12)
    assert sum(m.gas_elements.values()) == pytest.approx(1.0, abs=1e-12)
    assert min(m.elements.values()) >= 0.0 and min(m.gas_elements.values()) >= 0.0
    for symbol, value in m.dust_elements.items():
        assert 0.0 <= value <= m.elements[symbol] + 1e-15, symbol
    assert m.ionized + m.atomic + m.molecular == pytest.approx(1.0, abs=1e-12)
    assert min(m.ionized, m.atomic, m.molecular) >= 0.0
    if m.grains:
        assert sum(m.grains.values()) == pytest.approx(1.0)
    assert set(m.abundances) == set(nc.ABUNDANCE_ELEMENTS)


@pytest.mark.parametrize("letter", CLASS_LETTERS)
def test_no_element_goes_into_dust_past_its_cap(letter):
    m = make(letter, density=1e3 if letter in COLD else 100.0)
    for symbol, limit in nc.MAX_IN_DUST.items():
        assert m.dust_elements.get(symbol, 0.0) <= limit * m.elements[symbol] + 1e-15, symbol
    assert 0.0 <= m.dust_to_gas <= m.dust_to_gas_law + 1e-15


# --- the place ---------------------------------------------------------------------------------------------------------------------------

def test_a_metal_poor_place_gives_a_metal_poor_nebula_in_every_class():
    for letter in "ADEFGMNX":
        rich, poor = make(letter, 0.2), make(letter, -0.7)
        assert poor.z < rich.z and poor.dust_to_gas < rich.dust_to_gas and poor.abundances["O"] < rich.abundances["O"]
        assert poor.y < rich.y + 0.03 and poor.x > rich.x or letter == "X"


def test_the_oxygen_of_a_nebula_follows_the_place():
    for oh in (-0.66, -0.34, 0.0, 0.15):
        m = make("A", oh)
        assert nc.oh_dex(m.elements) == pytest.approx(oh, abs=1e-6)
        assert nc.oxygen_12log(m.elements) == pytest.approx(8.69 + oh, abs=1e-6)


def test_the_gas_gradient_is_its_own_and_zero_at_the_reference_radius():
    from planetgen.galaxy import scales
    reference = scales.reference_radius_pc()
    assert nc.gas_gradient_oh(reference) == pytest.approx(0.0)
    assert nc.gas_gradient_oh(reference + 5000.0) == pytest.approx(-5.0 * nc.GAS_GRADIENT_DEX_PER_KPC)
    assert nc.GAS_GRADIENT_DEX_PER_KPC != tuning.METALLICITY_DISK_FALL_DEX_PER_KPC
    assert 0.03 <= nc.GAS_GRADIENT_DEX_PER_KPC <= 0.04          # Boss: 0.03 to 0.04


def test_the_gas_of_a_place_is_the_same_every_time_and_scatters_around_its_gradient():
    reference = __import__("planetgen.galaxy.scales", fromlist=["x"]).reference_radius_pc()
    assert nc.gas_oh(reference, "n1") == nc.gas_oh(reference, "n1")
    values = [nc.gas_oh(reference, f"nebula-{i}") for i in range(600)]
    mean = sum(values) / len(values)
    sd = math.sqrt(sum((v - mean) ** 2 for v in values) / len(values))
    assert abs(mean) < 0.02 and sd == pytest.approx(nc.GAS_SCATTER_DEX, abs=0.02)


def test_a_satellite_gas_lies_inside_the_published_ranges_of_the_two_clouds():
    lmc = [nc.gas_oh(0.0, f"c{i}", component_fe_h=-0.4) for i in range(200)]
    smc = [nc.gas_oh(0.0, f"c{i}", component_fe_h=-0.9) for i in range(200)]
    assert -0.34 - 1e-9 <= min(lmc) and max(lmc) <= -0.14 + 1e-9
    assert -0.66 - 1e-9 <= min(smc) and max(smc) <= -0.28 + 1e-9
    assert sum(smc) / 200 < sum(lmc) / 200
    assert max(smc) - min(smc) > max(lmc) - min(lmc)            # the wider range of the SMC


def test_a_cloud_is_more_metal_poor_than_the_solar_circle():
    solar = make("D", nc.gas_gradient_oh(__import__("planetgen.galaxy.scales", fromlist=["x"]).reference_radius_pc()))
    smc = make("D", nc.gas_oh(0.0, "smc-nebula", component_fe_h=-0.9))
    assert smc.z < 0.5 * solar.z and smc.dust_to_gas < solar.dust_to_gas


# --- the dust ----------------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("x, gas_to_dust", [(8.69, 162.0), (8.35, 355.0), (8.03, 741.0), (7.70, 5360.0)])
def test_the_dust_law_gives_the_points_of_its_source(x, gas_to_dust):
    assert nc.gas_to_dust_law(x) == pytest.approx(gas_to_dust, rel=0.02)
    assert nc.dust_to_gas_law(x) == pytest.approx(1.0 / gas_to_dust, rel=0.02)


def test_the_dust_law_is_continuous_at_its_break():
    below, above = nc.gas_to_dust_law(nc.DUST_BREAK_12LOG - 1e-9), nc.gas_to_dust_law(nc.DUST_BREAK_12LOG + 1e-9)
    assert above == pytest.approx(below, rel=0.03)


def test_the_dust_of_a_class_follows_what_survives_in_it():
    diffuse, hii, shock = make("A"), make("D"), make("S")
    assert diffuse.dust_to_gas > hii.dust_to_gas > shock.dust_to_gas
    assert make("R").dust_to_gas < hii.dust_to_gas
    assert make("M", density=1e4).dust_to_gas > hii.dust_to_gas
    assert "ice" in make("M", density=1e4).grains and "ice" not in make("A").grains
    assert make("K").grains["carbon"] > make("K").grains["silicate"]
    assert make("H").grains["silicate"] > make("H").grains["carbon"]


def test_iron_stays_in_the_gas_of_a_carbon_rich_nebula_and_leaves_that_of_a_dusty_one():
    assert make("A").abundances["Fe"] < make("K").abundances["Fe"]
    orion = make("D", 0.0)
    assert 6.0 < orion.abundances["Fe"] < 7.2          # Orion's gas holds 6.07 against 7.5 in the Sun: iron mostly on grains
    assert 8.4 < orion.abundances["O"] < 8.75


def test_the_dust_that_survives_the_shock_of_a_remnant_is_the_sourced_share():
    shell = nc.CLASSES["S"].survival
    assert shell["silicate"] == pytest.approx(0.39) and shell["iron"] == pytest.approx(0.60)
    assert nc.CLASSES["D"].survival["iron"] == pytest.approx(0.9)


# --- the phases ----------------------------------------------------------------------------------------------------------------------------

def test_the_molecular_law_reaches_its_limits():
    assert nc.molecular_fraction(1e4, 1.0) > 0.95
    assert nc.molecular_fraction(1.0, 1.0) == 0.0                # a thin solar column is atomic
    assert nc.molecular_fraction(0.0, 1.0) == 0.0 and nc.molecular_fraction(100.0, 0.0) == 0.0
    assert 20.0 < next(s for s in range(1, 200) if nc.molecular_fraction(float(s), 1.0) > 0.5) < 45   # half molecular near 27 to 34 solar masses per pc2


def test_a_metal_poor_cloud_holds_less_molecular_gas_than_a_solar_one_of_the_same_column():
    for sigma in (30.0, 60.0, 100.0):
        solar, lmc, smc = (nc.molecular_fraction(sigma, z) for z in (1.0, 0.46, 0.22))
        assert solar > lmc >= smc
    assert nc.molecular_fraction(30.0, 0.22) == 0.0 < nc.molecular_fraction(30.0, 1.0)


def test_the_column_of_a_cloud_is_four_thirds_of_density_times_radius():
    sigma = nc.column_density_msun_pc2(100.0, 32.616)         # 10 pc
    assert sigma == pytest.approx(46.0, rel=0.05)
    assert nc.column_density_msun_pc2(200.0, 100.0) == pytest.approx(2.0 * nc.column_density_msun_pc2(100.0, 100.0))


@pytest.mark.parametrize("letter", COLD)
def test_a_cold_cloud_is_neutral_and_its_molecular_share_rises_with_its_column(letter):
    thin = make(letter, density=10.0, radius_ly=2.0)
    thick = make(letter, density=1e5, radius_ly=30.0)
    assert thin.ionized == thick.ionized == 0.0
    assert thin.molecular < thick.molecular and thick.molecular > 0.9
    poor = make(letter, -0.7, density=100.0, radius_ly=20.0)
    rich = make(letter, 0.0, density=100.0, radius_ly=20.0)
    assert poor.molecular <= rich.molecular


def test_a_warm_neutral_cloud_has_no_molecular_gas_and_an_h_ii_region_is_ionized():
    assert make("A", density=10.0, radius_ly=100.0).molecular == 0.0
    assert make("D").ionized > 0.9 and make("B").ionized > 0.8
    assert make("G").molecular == pytest.approx(0.20)


# --- enrichment ---------------------------------------------------------------------------------------------------------------------------

def test_a_wolf_rayet_ring_is_richer_in_helium_than_an_h_ii_region_of_the_same_place():
    ring, plain = make("X", wolf_rayet_subtype="WN6"), make("D")
    assert ring.y > plain.y + 0.1 and ring.elements["N"] > plain.elements["N"]
    wc = make("X", wolf_rayet_subtype="WC")
    assert wc.elements["C"] > ring.elements["C"] and wc.x < plain.x


def test_a_planetary_nebula_takes_the_ejecta_of_its_progenitor():
    ejecta = pne.progenitor("K", 0.0, "nebula-k").composition.elements
    carbon_rich = make("K", ejecta=ejecta)
    plain = make("K")
    assert carbon_rich.elements["C"] > plain.elements["C"]
    assert carbon_rich.abundances["C"] > carbon_rich.abundances["O"] - 0.01
    assert sum(carbon_rich.elements.values()) == pytest.approx(1.0)
    assert "ejecta" in plain.open_items[0]
    assert make("D").open_items == ()


def test_a_type_ia_remnant_with_no_hydrogen_says_so():
    ejecta = {"H": 0.0, "He": 0.0, "Fe": 0.5, "Si": 0.2, "S": 0.1, "Ca": 0.05, "O": 0.1, "C": 0.05}
    w = make("W", ejecta=ejecta)
    assert w.x == 0.0 and w.elements["Fe"] == pytest.approx(0.5)
    assert w.dust_to_gas >= 0.0 and math.isnan(w.abundances["O"])        # no hydrogen: an abundance against hydrogen has no meaning


# --- inputs ---------------------------------------------------------------------------------------------------------------------------------

def test_an_unknown_class_is_refused():
    with pytest.raises(ValueError):
        nc.makeup("Z", 0.0, 100.0, 10.0)


def test_a_field_cloud_gets_a_density_of_its_cell_inside_its_class_range_and_the_same_each_time():
    low, high = tuning.NEBULA_CLASSES["N"]["density_range_cm3"]
    values = [nc.field_cloud_density_cm3("N", f"cell-{i}") for i in range(300)]
    assert all(low <= v <= high for v in values) and len(set(values)) > 250
    assert nc.field_cloud_density_cm3("N", "cell-1") == nc.field_cloud_density_cm3("N", "cell-1")


def test_the_values_a_composition_item_stores_are_all_numbers():
    values = make("D").values()
    assert {"mass_fraction_h", "mass_fraction_he", "mass_fraction_z", "mass_fraction_o", "dust_to_gas", "molecular_share", "abundance_fe"} <= set(values)
    assert all(isinstance(v, float) and math.isfinite(v) for v in values.values())


def test_a_galaxy_with_a_very_wide_disk_caps_the_gas_metallicity_and_every_class_still_has_a_makeup(monkeypatch):
    """The gradient is zero at the reference radius, so a reference radius of 95 kpc would put the gas at the center above +3 dex, where the dust took more
    than half the hydrogen and `makeup` raised (the `planetgen galaxy` runs of the fuzz test on the flat test galaxy)."""
    from planetgen.galaxy import scales

    monkeypatch.setattr(scales, "reference_radius_pc", lambda: 96000.0)
    assert nc.gas_gradient_oh(0.0) > 3.0
    for key in range(50):
        assert nc.gas_oh(0.0, f"wide-{key}") <= nc.GAS_OH_MAX_DEX
    for nebula_class in ("M", "N"):
        for density in (0.01, 1.0, 1e4):
            made = nc.makeup(nebula_class, nc.gas_oh(0.0, "wide-0"), density, 10.0)
            assert 0.0 < made.dust_to_gas < 10.0


# --- the dust species (GEN.1376) --------------------------------------------------------------------------------------------------------------

def test_the_species_split_each_kind_without_losing_any_dust():
    assert set(nc.SPECIES) == {"amorphous_olivine", "amorphous_pyroxene", "crystalline_silicate", "iron", "amorphous_carbon", "graphite", "PAH", "SiC", "ice"}
    for letter in CLASS_LETTERS:
        m = make(letter, density=1e3 if letter in COLD else 100.0)
        if not m.grains:
            continue
        assert set(m.species) == set(nc.SPECIES), letter
        assert sum(m.species.values()) == pytest.approx(1.0), letter
        silicate = m.species["amorphous_olivine"] + m.species["amorphous_pyroxene"] + m.species["crystalline_silicate"]
        assert silicate == pytest.approx(m.grains.get("silicate", 0.0)), letter
        assert m.species["amorphous_carbon"] + m.species["graphite"] == pytest.approx(m.grains.get("carbon", 0.0)), letter
        assert m.species["SiC"] == pytest.approx(m.grains.get("SiC", 0.0)), letter


def test_the_amorphous_silicate_is_85_percent_olivine_and_the_diffuse_medium_is_barely_crystalline():
    """The diffuse-medium silicate fit [R5143]: 84.9 percent amorphous olivine and 15.1 percent pyroxene by mass; 0.2 +- 0.2 percent crystalline."""
    m = make("B")
    amorphous = m.species["amorphous_olivine"] + m.species["amorphous_pyroxene"]
    assert m.species["amorphous_olivine"] / amorphous == pytest.approx(0.849)
    assert m.species["crystalline_silicate"] / m.grains["silicate"] == pytest.approx(0.002)


def test_the_stellar_ejecta_silicates_are_11_to_18_percent_crystalline():
    for letter in ("H", "J", "L"):
        m = make(letter)
        assert 0.11 <= m.species["crystalline_silicate"] / m.grains["silicate"] <= 0.18, letter


def test_the_silicon_in_interstellar_sic_is_under_five_percent_of_the_silicon_in_silicates():
    """The dust review [R5144], section 5.3: Si in SiC dust under 5 percent of Si in silicates."""
    for name in ("diffuse", "dense", "reflection"):
        mix = nc.GRAIN_MIXES[name]
        in_sic = mix["SiC"] * nc.GRAINS["SiC"]["Si"]
        in_silicate = mix["silicate"] * nc.GRAINS["silicate"]["Si"]
        assert in_sic <= 0.05 * in_silicate, name
        assert in_sic == pytest.approx(nc.SIC_SI_SHARE_OF_SILICATE_SI * in_silicate, rel=0.05), name


def test_a_carbon_star_wind_is_the_only_mix_with_a_tenth_of_its_dust_in_sic():
    assert nc.GRAIN_MIXES["carbon-rich"]["SiC"] == pytest.approx(0.10)
    assert all(mix.get("SiC", 0.0) < 0.01 for name, mix in nc.GRAIN_MIXES.items() if name != "carbon-rich")
    assert make("K").species["SiC"] > 10 * make("B").species["SiC"]
    assert make("K").species["amorphous_carbon"] > make("K").species["graphite"]


def test_the_ice_mantle_is_boogerts_and_gives_the_ice_its_carbon_and_nitrogen():
    assert nc.ICE_MANTLE["H2O"] == pytest.approx(0.4419, abs=1e-3)
    assert sum(nc.ICE_MANTLE.values()) == pytest.approx(1.0)
    assert nc.ICE_MANTLE["CO2"] > nc.ICE_MANTLE["CO"] > nc.ICE_MANTLE["CH3OH"] > nc.ICE_MANTLE["NH3"] > nc.ICE_MANTLE["CH4"]
    assert {"H", "O", "C", "N"} == set(nc.GRAINS["ice"]) and nc.GRAINS["ice"]["C"] > 0.15
    cloud = make("M", density=1e4)
    assert cloud.species["ice"] > 0.2 and cloud.dust_elements["N"] > 0.0


def test_sic_survives_as_the_silicates_do():
    shocked = make("S")
    assert shocked.grains["SiC"] / shocked.grains["silicate"] == pytest.approx(make("A").grains["SiC"] / make("A").grains["silicate"] * 1.0, rel=0.3)
    assert nc.SURVIVES_AS["SiC"] == "silicate"
