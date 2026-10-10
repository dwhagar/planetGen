# tests/test_habitability.py

"""
The habitability index's reference maths (GEN.84,
`planetgen.physics.habitability`): the worked examples in
docs/design/habitability-index.md section 2 reproduce, and the bands and
derived constants hold.
"""

import pytest

from planetgen.physics import habitability as hab
from planetgen.physics.habitability import World, scores


EARTH = World(101.3, 15, {"N2": 0.78, "O2": 0.21, "CO2": 0.0004}, ph=8.1, relative_humidity=0.7)

EXAMPLES = {
    "Earth": (EARTH, ("Blue", "Blue", "Blue", "Blue"), 1.0, 0.958, 0.957, 1.0, "ideal"),
    "Dense CO2": (World(300, 40, {"CO2": 0.95, "N2": 0.05}, ph=6.5, relative_humidity=0.8, energy_flux_w_m2=40),
                  ("Green", "Yellow", "Red", "Blue"), 0.0, 0.958, 0.0, 0.479, "sealed suit"),
    "Hycean": (World(2000, 60, {"H2": 0.9, "He": 0.1}, ph=7, relative_humidity=0.99, energy_flux_w_m2=20),
               ("Yellow", "Yellow", "Green", "Blue"), 0.7364, 0.958, 0.0, 0.022, "sealed suit"),
    "Photochemical CO": (World(100, 10, {"N2": 0.9, "CO2": 0.05, "CO": 0.05}, ph=7.5),
                         ("Blue", "Blue", "Red", "Blue"), 0.5891, 0.958, 0.0, 1.0, "sealed suit"),
    "Ice-sealed ocean": (World(1e-6, -170, {}, gravity_ms2=1.3, ph=8, energy_flux_w_m2=1e-3, surface_dose_sv_yr=0.05,
                               habitat_dose_sv_yr=0.002, water_source="ice", relative_humidity=0.05),
                         ("Red", "Red", "Green", "Blue"), 0.0, 0.41, 0.0, 0.468, "sealed suit"),
    "Mars-analog brine": (World(0.6, -60, {"CO2": 0.95, "N2": 0.03}, gravity_ms2=3.71, water_activity=0.5,
                                chaotropicity_kj_kg=90, ph=8, inventories={"C": 0.05, "N": 0.01, "H": 0.05},
                                phosphate_umol_l=0.5, energy_flux_w_m2=40, water_source="ice",
                                relative_humidity=0.05),
                          ("Red", "Red", "Red", "Yellow"), 0.0, 0.0, 0.0, 0.458, "sealed suit"),
    "Dune world": (World(90, 45, {"N2": 0.78, "O2": 0.21, "CO2": 0.001}, relative_humidity=0.05, water_activity=0.7,
                         ph=8, water_source="vapour", inventories={"C": 0.5, "N": 0.5, "H": 0.05}),
                   ("Blue", "Green", "Yellow", "Blue"), 0.8289, 0.54, 0.539, 0.795, "ideal"),
    "Europan ocean": (World(1e-9, -160, {}, gravity_ms2=1.31, ph=9, energy_flux_w_m2=1e-4, surface_dose_sv_yr=2000,
                            habitat_dose_sv_yr=0.002, water_source="ice", relative_humidity=0.05),
                      ("Red", "Red", "Green", "Red"), 0.0, 0.274, 0.0, 0.276,
                      "full life support with radiation hardening"),
}


@pytest.mark.parametrize("name", list(EXAMPLES))
def test_the_worked_examples_reproduce(name):
    world, tiers, phi4, bio, cpx, tech, kit = EXAMPLES[name]
    got = scores(world)
    assert tuple(got[d][1] for d in ("Pressure", "Temperature", "Chemistry", "Radiation")) == tiers
    assert got["PHI-4"] == pytest.approx(phi4, abs=1e-3)
    assert (got["PHI_bio"], got["PHI_cpx"], got["Phi_tech"]) == pytest.approx((bio, cpx, tech), abs=2e-3)
    assert got["equipment"] == kit


def test_the_mask_limit_is_derived_from_the_armstrong_limit():
    assert hab.MASK_MIN_KPA == pytest.approx(14.3)
    pure_o2 = World(hab.MASK_MIN_KPA, 20, {"O2": 1.0})
    assert hab.inspired_o2_kpa(pure_o2) == pytest.approx(hab.INSPIRED_O2_MIN_KPA)


def test_band_edges_score_their_tier_boundaries():
    band = hab.TEMPERATURE_BAND
    assert hab.band_score(10, band) == 1.0
    assert hab.band_score(45, band) == pytest.approx(0.75)
    assert hab.band_score(122, band) == pytest.approx(0.4)
    assert hab.band_score(200, band) == 0.0
    assert hab.band_score(-50, band) == pytest.approx(0.4)
    assert [hab.tier(s) for s in (1.0, 0.8, 0.5, 0.2)] == ["Blue", "Green", "Yellow", "Red"]


def test_surface_dose_fits_the_moon_mars_and_earth():
    assert hab.gcr_surface_dose_sv_yr(0) == pytest.approx(0.502)
    mars = hab.gcr_surface_dose_sv_yr(hab.column_mass_g_cm2(0.6, 3.71))
    assert 0.2 < mars < 0.3
    earth = hab.gcr_surface_dose_sv_yr(hab.column_mass_g_cm2(101.3, 9.81))
    assert earth == pytest.approx(0.0024, abs=2e-4)


def test_an_uncompensable_wet_bulb_is_never_better_than_yellow():
    humid = World(40, 38, {"N2": 0.79, "O2": 0.21}, relative_humidity=0.9)
    assert hab.wet_bulb_c(38, 0.9) > hab.WET_BULB_CRIT_C
    assert hab.tier(hab.temperature_score(humid)) == "Yellow"


def test_no_solvent_means_no_microbial_life():
    dry = World(101.3, 15, {"N2": 0.79, "O2": 0.21}, solvent=None)
    assert hab.phi_bio(dry) == 0.0
    assert hab.phi_cpx(dry) == 0.0
