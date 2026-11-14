# tests/test_white_dwarf_core.py

"""
GEN.1370 slice 0: the central composition of a carbon-oxygen white dwarf (oxygen by mass, 22Ne by the progenitor's CNO) and the elements and class
of the body a polluted white dwarf accreted. The functions are pure and not yet wired into the stored white dwarf.
"""

import math

import pytest

from planetgen import tuning
from planetgen.physics import condensate_shares, materials
from planetgen.physics import small_body_composition as sb
from planetgen.physics import stellar_abundances as sa
from planetgen.physics import white_dwarf_core as core
from planetgen.util import draw

SUN_W = 3.828e26
WATER_LINE_AU = condensate_shares.ice_line_au(SUN_W)


def _belts(x, count, seed):
    with draw.bound(seed):
        return [sb.draw_belt(x * WATER_LINE_AU, SUN_W, 5.2e-5) for _ in range(count)]


# --- the core ----------------------------------------------------------------------------------------------------

def test_the_tabulated_masses_return_the_tabulated_oxygen():
    for mass, oxygen in tuning.WD_CENTRAL_OXYGEN_BY_MASS_SOL:
        assert core.central_oxygen(mass) == pytest.approx(oxygen)


def test_the_oxygen_is_linear_between_the_points_and_held_outside_them():
    (m0, y0), (m1, y1) = tuning.WD_CENTRAL_OXYGEN_BY_MASS_SOL[:2]
    assert core.central_oxygen((m0 + m1) / 2) == pytest.approx((y0 + y1) / 2)
    assert core.central_oxygen(m0 - 0.01) == pytest.approx(y0)
    assert core.central_oxygen(0.95) == pytest.approx(tuning.WD_CENTRAL_OXYGEN_BY_MASS_SOL[-1][1])


def test_the_oxygen_stays_inside_the_asteroseismic_range():
    low, high = tuning.WD_CENTRAL_OXYGEN_RANGE
    for i in range(0, 41):
        mass = 0.45 + 0.01 * i
        if sa.white_dwarf_core(mass) == "carbon-oxygen":
            assert low <= core.central_oxygen(mass) <= high


def test_a_helium_or_oxygen_neon_core_has_no_carbon_oxygen_split():
    for mass in (0.3, 1.2):
        assert sa.white_dwarf_core(mass) != "carbon-oxygen"
        with pytest.raises(core.WhiteDwarfCoreError):
            core.central_oxygen(mass)


def test_the_neon_of_the_proto_solar_mix_is_the_sum_of_its_cno():
    birth = sa.composition(0.0, {})
    by_hand = sum(birth.elements[s] * 22.0 / core.elements.element(s).atomic_weight for s in ("C", "N", "O"))
    assert core.neon22_fraction(birth) == pytest.approx(by_hand * tuning.WD_NEON22_PER_CNO)
    assert 0.012 < core.neon22_fraction(birth) < 0.018


def test_the_neon_rises_with_the_metallicity_of_the_progenitor():
    poor, solar, rich = (core.neon22_fraction(sa.composition(f, {})) for f in (-1.0, 0.0, 0.3))
    assert poor < solar < rich
    assert rich / solar == pytest.approx(10 ** 0.3, rel=0.1)


def test_a_composition_without_cno_makes_no_neon():
    class Bare:
        elements = {"Fe": 1.0}
    with pytest.raises(core.WhiteDwarfCoreError):
        core.neon22_fraction(Bare())


def test_the_core_fractions_add_to_one_and_keep_the_oxygen_to_carbon_ratio():
    birth = sa.composition(0.0, {})
    for mass in (0.55, 0.6, 0.7, 0.85):
        fractions = core.core_fractions(mass, birth)
        assert fractions.oxygen + fractions.carbon + fractions.neon22 == pytest.approx(1.0)
        oxygen = core.central_oxygen(mass)
        assert fractions.oxygen / (fractions.oxygen + fractions.carbon) == pytest.approx(oxygen)
        assert fractions.neon22 > 0.0


# --- the pollutant -----------------------------------------------------------------------------------------------

def test_every_material_a_belt_draws_states_its_elements():
    seen = set()
    for x in (0.3, 1.0, 2.0, 4.0):
        for belt in _belts(x, 40, 3):
            seen.update(belt.components)
    assert seen
    assert not [name for name in sorted(seen) if core._elements_of(name) is None]


def test_the_element_fractions_of_a_material_add_to_one_and_follow_its_formula():
    olivine = core._elements_of("forsterite")
    assert sum(olivine.values()) == pytest.approx(1.0)
    atoms = core.elements.parse_formula(materials.material("forsterite").formula)
    mass = {s: n * core.elements.element(s).atomic_weight for s, n in atoms.items()}
    assert olivine["Mg"] == pytest.approx(mass["Mg"] / sum(mass.values()))
    for name, shares in tuning.POLLUTANT_ELEMENT_SHARES.items():
        assert sum(shares.values()) == pytest.approx(1.0), name
        assert core._elements_of(name) == shares


def test_the_end_member_formulas_parse_and_hold_what_the_table_names():
    for name, formula in tuning.POLLUTANT_FORMULAS.items():
        atoms = core.elements.parse_formula(formula)
        assert atoms, name
        assert math.isclose(sum(core._elements_of(name).values()), 1.0)
    assert {"Fe", "Ni", "P"} <= set(core._elements_of("schreibersite"))


def test_the_pollutant_fractions_add_to_one():
    for x in (0.3, 1.0, 2.0):
        for belt in _belts(x, 12, 4):
            found = core.pollutant(belt.components)
            assert sum(found.element_fractions.values()) + found.unassigned == pytest.approx(1.0)
            assert found.unassigned < 1e-9
            assert found.kind in core.CLASSES
            assert found.ratios_to_mg["Mg"] == pytest.approx(1.0)


def test_an_iron_nickel_body_is_core_rich_and_has_no_magnesium_to_scale_to():
    found = core.pollutant({"kamacite": 0.9, "taenite": 0.1})
    assert found.kind == core.CORE_RICH
    assert found.ratios_to_mg == {}
    assert found.element_fractions["Fe"] > 0.8 and found.element_fractions["Ni"] == pytest.approx(0.9 * 0.06 + 0.1 * 0.28)


def test_an_icy_body_is_volatile_rich():
    found = core.pollutant({"forsterite": 0.6, "water_ice": 0.4})
    assert found.kind == core.VOLATILE_RICH
    assert found.element_fractions["H"] > 0.04


def test_a_peridotite_body_is_primitive_and_a_calcium_aluminum_rich_one_is_crust_rich():
    assert core.pollutant({"forsterite": 0.7, "enstatite": 0.3}).kind == core.PRIMITIVE
    crust = core.pollutant({"melilite": 0.6, "hibonite": 0.1, "forsterite": 0.3})
    assert crust.kind == core.CRUST_RICH
    assert crust.ratios_to_mg["Ca"] > 1.0 and crust.ratios_to_mg["Al"] > 1.0


def test_a_carbonaceous_belt_shows_carbon_and_a_dry_one_does_not():
    carbon = core.pollutant({"forsterite": 0.6, "complex_organics": 0.2, "tholin": 0.2})
    assert carbon.ratios_to_mg["C"] > 0.1 and "N" not in carbon.ratios_to_mg
    assert "C" not in core.pollutant({"forsterite": 1.0}).ratios_to_mg


def test_a_cold_belt_is_mostly_volatile_rich_and_a_warm_one_is_not():
    cold = [core.pollutant(b.components).kind for b in _belts(2.0, 30, 6)]
    warm = [core.pollutant(b.components).kind for b in _belts(0.5, 30, 6)]
    assert cold.count(core.VOLATILE_RICH) > len(cold) / 2
    assert core.VOLATILE_RICH not in warm


def test_an_empty_body_or_an_unknown_material_is_an_error():
    with pytest.raises(core.WhiteDwarfCoreError):
        core.pollutant({})
    with pytest.raises(materials.UnknownMaterial):
        core.pollutant({"no_such_mineral": 1.0})


# --- the accreted body -------------------------------------------------------------------------------------------

def test_the_accreted_member_is_picked_in_proportion_to_its_mass():
    members = [(1.0, "C"), (3.0, "S"), (0.0, "M")]
    picks = [core.accreted_taxon(members, (i + 0.5) / 1000) for i in range(1000)]
    assert picks.count("C") == 250 and picks.count("S") == 750 and "M" not in picks
    assert core.accreted_taxon(members, 0.999999) == "S" and core.accreted_taxon(members, 0.0) == "C"
    assert core.accreted_taxon([], 0.5) is None and core.accreted_taxon([(0.0, "C")], 0.5) is None


def test_every_taxon_of_a_belt_has_a_pollutant():
    kinds = {taxon: core.pollutant_of_taxon(taxon).kind for taxon in sb.TAXA}
    assert set(kinds.values()) <= set(core.CLASSES)
    assert kinds["M"] == core.CORE_RICH and kinds["S"] == core.PRIMITIVE
    for taxon in sb.TAXA:
        assert core.pollutant_of_taxon(taxon).unassigned < 1e-9, taxon
