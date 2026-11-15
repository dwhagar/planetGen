# tests/test_rock_assemblage.py

"""GEN.474, part 2: the assemblage of a rock is its norm on the rows of the materials table (`physics/rock_assemblage.py`). The Earth's and the Moon's analyses are test points, never inputs."""

import pytest

from planetgen.physics import materials, mixture_lava, rock_assemblage as ra, rock_types as rt

BULK_SILICATE_EARTH = {oxide: row[0] for oxide, row in mixture_lava.OXIDES.items()}     # McDonough and Sun 1995
UPPER_CRUST = dict(SiO2=66.6, TiO2=0.64, Al2O3=15.4, FeO=5.04, MgO=2.48, CaO=3.59, Na2O=3.27, K2O=2.80)     # [R3842: Table 11]
HIGHLAND_ANORTHOSITE = dict(SiO2=44.0, Al2O3=35.5, FeO=0.2, MgO=0.2, CaO=19.0, Na2O=0.4)     # a ferroan anorthosite of the lunar highlands (recalled)
ALKALI_BASALT = dict(SiO2=44.0, TiO2=2.5, Al2O3=14.0, Fe2O3=3.0, FeO=9.0, MgO=9.0, CaO=10.0, Na2O=3.5, K2O=1.5, P2O5=0.5)     # a silica-undersaturated basalt (illustrative)
PERALKALINE = dict(SiO2=70.0, Al2O3=8.0, FeO=3.0, Na2O=12.0, K2O=5.0)      # more alkali than alumina: metasilicates


@pytest.mark.parametrize("rock", [BULK_SILICATE_EARTH, UPPER_CRUST, HIGHLAND_ANORTHOSITE, ALKALI_BASALT])
def test_the_shares_add_to_one_and_every_one_is_a_row_of_the_materials_table(rock):
    shares, unmapped = ra.assemblage(rock)
    assert sum(shares.values()) == pytest.approx(1.0, rel=1e-12)
    assert unmapped == pytest.approx(0.0, abs=1e-12)
    assert all(materials.has(name) for name in shares)


def test_the_mantle_is_an_olivine_rock_and_the_upper_crust_is_quartz_and_feldspar():
    mantle, _ = ra.assemblage(BULK_SILICATE_EARTH)
    assert mantle["olivine"] > 0.5 and "quartz" not in mantle
    crust, _ = ra.assemblage(UPPER_CRUST)
    assert 0.18 < crust["quartz"] < 0.26
    assert crust["albite"] + crust["anorthite"] + crust["orthoclase"] > 0.5


def test_a_highland_anorthosite_is_nearly_all_plagioclase():
    shares, _ = ra.assemblage(HIGHLAND_ANORTHOSITE)
    assert shares["anorthite"] + shares["albite"] > 0.97


def test_the_orthopyroxene_is_split_by_the_iron_of_the_pool():
    mantle, _ = ra.assemblage(BULK_SILICATE_EARTH)
    crust, _ = ra.assemblage(UPPER_CRUST)
    assert mantle["enstatite"] / (mantle["enstatite"] + mantle["ferrosilite"]) > 0.8      # a magnesian pool
    assert crust["ferrosilite"] > crust["enstatite"]                                        # an iron-rich crust
    assert ra.magnesian_mass_share(dict(SiO2=50.0, MgO=10.0)) == pytest.approx(1.0)
    assert 0.0 < ra.magnesian_mass_share(UPPER_CRUST) < 1.0


def test_the_norm_minerals_without_a_row_are_returned_as_the_unmapped_share():
    shares, unmapped = ra.assemblage(PERALKALINE)
    assert 0.0 < unmapped < 1.0
    assert sum(shares.values()) == pytest.approx(1.0, rel=1e-12)


def test_the_assemblage_of_the_melt_of_a_mantle_is_the_norm_of_the_melt():
    melt, _residue = rt.partial_melt(BULK_SILICATE_EARTH, 0.10)
    shares, _ = ra.assemblage(melt)
    assert rt.ground_type_of(melt) == "basaltic"
    assert shares["anorthite"] > shares["albite"] > 0.0 and shares["diopside"] > 0.0 and shares.get("olivine", 0.0) < 0.01


def test_a_rock_the_norm_refuses_has_no_assemblage():
    with pytest.raises(ValueError):
        ra.assemblage(dict(SiO2=1.0, K2O=10.0, Al2O3=10.0))
