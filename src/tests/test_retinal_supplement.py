# tests/test_retinal_supplement.py

"""GEN.881: retinal pigments are a light supplement to a chemical class that fixes carbon in sunlit water (`life_classes.retinal_supplement_w_m2`); they add energy and
fix no carbon, and a class that eats light already, an acid eater, a class under ice or a world with no star gets none."""

from types import SimpleNamespace

import pytest

from planetgen.physics import life_classes as lc

SUN = {"_teff_k": lc.REFERENCE_STAR_K}
LIT = SimpleNamespace(solvent="water", life_on_surface=True, energy_flux_w_m2=80.0)
SEA = {"hydrosphere": "surface ocean"}


def _supplement(class_id, found=SUN, world=LIT, values=SEA):
    return lc.retinal_supplement_w_m2(lc.BY_ID[class_id], found, world, values)


@pytest.mark.parametrize("class_id", ["3a", "3b", "4a", "4b", "4c"])
def test_a_chemical_class_that_fixes_carbon_in_sunlit_water_gets_its_share_of_the_light(class_id):
    assert lc.BY_ID[class_id].life["metabolism"] in lc.CARBON_FIXERS
    assert _supplement(class_id) == pytest.approx(lc.RETINAL_SHARE * 80.0)


def test_the_supplement_follows_the_light_in_its_band_not_the_total():
    cool = _supplement("4a", {"_teff_k": 3000.0})
    hot = _supplement("4a", {"_teff_k": 9000.0})
    assert 0.0 < cool < _supplement("4a") < hot


@pytest.mark.parametrize("class_id", ["1a", "1b", "1c", "2"])
def test_a_class_that_eats_light_already_gets_none(class_id):
    assert _supplement(class_id) == 0.0


@pytest.mark.parametrize("class_id", ["5a", "5b", "6", "7", "8", "9"])
def test_a_class_that_fixes_no_carbon_or_does_not_live_in_water_gets_none(class_id):
    assert _supplement(class_id) == 0.0


def test_no_light_reaches_water_under_ice_or_underground_or_from_no_star():
    assert _supplement("4a", values={"hydrosphere": "ice-covered ocean"}) == 0.0
    assert _supplement("4a", world=SimpleNamespace(solvent="water", life_on_surface=False, energy_flux_w_m2=80.0)) == 0.0
    assert _supplement("4a", world=SimpleNamespace(solvent="hydrocarbon", life_on_surface=True, energy_flux_w_m2=80.0)) == 0.0
    assert _supplement("4a", found={}) == 0.0


def test_the_supplement_raises_a_chemical_class_energy_term_and_no_other_term():
    row = lc.BY_ID["4a"]
    values = {"hydrosphere": "surface ocean", "internal_heat_flux_w_m2": 1e-4}
    world = SimpleNamespace(solvent="water", life_on_surface=True, energy_flux_w_m2=80.0)
    lit = lc._energy_term(row, SUN, world, values)
    dark = lc._energy_term(row, SUN, SimpleNamespace(solvent="water", life_on_surface=False, energy_flux_w_m2=80.0), values)
    assert lit > dark >= 0.0
