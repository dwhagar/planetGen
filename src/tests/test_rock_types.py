# tests/test_rock_types.py

"""GEN.474: the CIPW norm, the melt of a rock and the ground type of its oxides (`physics/rock_types.py`). The Earth's and the Moon's analyses are test points, never inputs."""

import pytest

from planetgen.physics import mixture_lava, rock_types as rt, surface_composition

BULK_SILICATE_EARTH = {oxide: row[0] for oxide, row in mixture_lava.OXIDES.items()}     # McDonough and Sun 1995
UPPER_CRUST = dict(SiO2=66.6, TiO2=0.64, Al2O3=15.4, FeO=5.04, MgO=2.48, CaO=3.59, Na2O=3.27, K2O=2.80)     # [R3842: Table 11]
MIDDLE_CRUST = dict(SiO2=63.5, TiO2=0.69, Al2O3=15.0, FeO=6.02, MgO=3.59, CaO=5.25, Na2O=3.39, K2O=2.30)
LOWER_CRUST = dict(SiO2=53.4, TiO2=0.82, Al2O3=16.9, FeO=8.57, MgO=7.24, CaO=9.59, Na2O=2.65, K2O=0.61)
TOTAL_CRUST = dict(SiO2=60.6, TiO2=0.72, Al2O3=15.9, FeO=6.71, MgO=4.66, CaO=6.41, Na2O=3.07, K2O=1.81)
HIGHLAND_ANORTHOSITE = dict(SiO2=44.0, Al2O3=35.5, FeO=0.2, MgO=0.2, CaO=19.0, Na2O=0.4)     # a ferroan anorthosite of the lunar highlands (recalled)
ALKALI_BASALT = dict(SiO2=44.0, TiO2=2.5, Al2O3=14.0, Fe2O3=3.0, FeO=9.0, MgO=9.0, CaO=10.0, Na2O=3.5, K2O=1.5, P2O5=0.5)     # a silica-undersaturated basalt (illustrative)


@pytest.mark.parametrize("rock", [BULK_SILICATE_EARTH, UPPER_CRUST, LOWER_CRUST, TOTAL_CRUST, HIGHLAND_ANORTHOSITE, ALKALI_BASALT])
def test_the_norm_adds_up_to_the_oxides(rock):
    assert sum(rt.cipw_norm(rock).values()) == pytest.approx(sum(rock.values()), rel=1e-9)


def test_the_feldspars_of_the_upper_crust_follow_from_its_alkalis():
    """Orthoclase is the K2O with its alumina and silica (K2O 94.2, KAlSi3O8 times two 556.7 per K2O), albite likewise: 2.8 K2O gives 16.6 of orthoclase."""
    norm = rt.cipw_norm(UPPER_CRUST)
    assert norm["orthoclase"] == pytest.approx(2.80 * (94.195 + 101.961 + 6 * 60.083) / 94.195, rel=1e-3)
    assert norm["albite"] == pytest.approx(3.27 * (61.979 + 101.961 + 6 * 60.083) / 61.979, rel=1e-3)
    assert 18.0 < norm["quartz"] < 26.0                      # a granodiorite is quartz-bearing
    assert "olivine" not in norm and "nepheline" not in norm


def test_the_mantle_norm_is_an_olivine_rock_and_an_undersaturated_basalt_has_nepheline():
    mantle = rt.cipw_norm(BULK_SILICATE_EARTH)
    assert mantle["olivine"] > 0.5 * sum(mantle.values()) and "quartz" not in mantle
    basalt = rt.cipw_norm(ALKALI_BASALT)
    assert basalt["nepheline"] > 1.0 and "quartz" not in basalt and basalt["apatite"] > 0.0


def test_a_rock_that_cannot_be_allocated_is_refused_not_trimmed():
    with pytest.raises(ValueError):
        rt.cipw_norm(dict(SiO2=50.0, Al2O3=5.0, Na2O=10.0))             # alkalis over alumina
    with pytest.raises(ValueError):
        rt.cipw_norm(dict(SiO2=1.0, K2O=10.0, Al2O3=10.0))              # no silica to speak of
    with pytest.raises(ValueError):
        rt.cipw_norm(dict(SiO2=50.0, MgO=-1.0))
    with pytest.raises(ValueError):
        rt.cipw_norm(dict(SiO2=50.0, H2O=5.0))


def test_the_earths_crusts_and_the_moons_highlands_get_their_ground_types():
    assert rt.ground_type_of(UPPER_CRUST) == "felsic"
    assert rt.ground_type_of(MIDDLE_CRUST) == "felsic"
    assert rt.ground_type_of(LOWER_CRUST) == "intermediate"
    assert rt.ground_type_of(TOTAL_CRUST) == "intermediate"
    assert rt.ground_type_of(HIGHLAND_ANORTHOSITE) == "anorthositic"
    assert rt.ground_type_of(ALKALI_BASALT) == "basaltic"
    assert rt.ground_type_of(BULK_SILICATE_EARTH) == "ultramafic"


def test_the_melt_of_a_mantle_is_a_basalt_at_the_degrees_of_melting_a_mantle_reaches():
    """A tenth of the bulk silicate Earth melts to a basaltic liquid; the liquid stays at the eutectic composition while the low-melting stage lasts."""
    melt, residue = rt.partial_melt(BULK_SILICATE_EARTH, 0.10)
    assert 45.0 <= melt["SiO2"] < 52.0 and rt.ground_type_of(melt) == "basaltic"
    assert residue["MgO"] > BULK_SILICATE_EARTH["MgO"] * 100.0 / sum(BULK_SILICATE_EARTH.values())
    early, _ = rt.partial_melt(BULK_SILICATE_EARTH, 0.03)
    assert early == pytest.approx(melt)


def test_a_great_degree_of_melting_makes_a_magnesian_liquid_of_komatiite_kind():
    melt, _ = rt.partial_melt(BULK_SILICATE_EARTH, 0.50)
    assert melt["MgO"] > rt.KOMATIITE_MAGNESIA and rt.ground_type_of(melt) == "ultramafic"
    assert rt.partial_melt(BULK_SILICATE_EARTH, 0.05)[0]["MgO"] < rt.partial_melt(BULK_SILICATE_EARTH, 0.50)[0]["MgO"]


def test_the_melt_and_the_residue_add_back_to_the_rock():
    source_total = sum(BULK_SILICATE_EARTH.values())
    for fraction in (0.05, 0.25, 0.6):
        melt, residue = rt.partial_melt(BULK_SILICATE_EARTH, fraction)
        assert sum(melt.values()) == pytest.approx(100.0) and sum(residue.values()) == pytest.approx(100.0)
        for oxide, value in BULK_SILICATE_EARTH.items():
            assert fraction * melt[oxide] + (1.0 - fraction) * residue[oxide] == pytest.approx(100.0 * value / source_total, abs=1e-9)


def test_a_crust_melts_to_a_granite_minimum_melt():
    """The quartz and feldspars of the upper crust melt first and together, so its melt is a felsic liquid at any degree of melting up to their exhaustion."""
    melt, _ = rt.partial_melt(UPPER_CRUST, 0.2)
    assert rt.ground_type_of(melt) == "felsic" and melt["SiO2"] > 68.0


def test_the_melt_fraction_and_the_partition_coefficient_are_bounded():
    for bad in (0.0, 1.5, -0.1):
        with pytest.raises(ValueError):
            rt.partial_melt(BULK_SILICATE_EARTH, bad)
        with pytest.raises(ValueError):
            rt.batch_enrichment(bad)
    with pytest.raises(ValueError):
        rt.batch_enrichment(0.1, partition=1.5)


def test_batch_melting_enriches_what_does_not_enter_the_minerals():
    assert rt.batch_enrichment(0.1) == pytest.approx(10.0)
    assert rt.batch_enrichment(0.1, partition=1.0) == pytest.approx(1.0)
    assert rt.batch_enrichment(1.0) == pytest.approx(1.0)
    assert rt.batch_enrichment(0.1, partition=0.5) == pytest.approx(1.0 / (0.5 + 0.1 * 0.5))


def test_the_heat_producers_of_the_earths_basalt_and_mantle_agree_with_one_degree_of_melting():
    """N-MORB over depleted mantle (K 22.8, Th 31.9, U 25.9 times) is the enrichment of one batch melt of 3 to 4.5 percent, the same F for the three: the elements
    behave as one incompatible group (a test point of the model, never an input)."""
    mantle, basalt = surface_composition.RESERVOIRS["ultramafic"], surface_composition.RESERVOIRS["basaltic"]
    for key in ("k_ppm", "th_ppm", "u_ppm"):
        fraction = 1.0 / (basalt[key] / mantle[key])             # the F whose enrichment 1 / F is the ratio, D being zero
        assert 0.025 < fraction < 0.045
        assert rt.batch_enrichment(fraction) == pytest.approx(basalt[key] / mantle[key])
