"""
The compressed water layer (GEN.461, `planetgen.physics.water_layer`): the depth and the base pressure against the integrated
reference of design part C 2.4, the liquidus against the recomputed melting curves, and determinism.
"""

import math
import pathlib

import pytest

from planetgen.physics import bulk_structure, constants, water_layer as wl

EARTH = constants.EARTH_MASS_TO_KG
R_EARTH_KM = constants.EARTH_RADIUS_KM

# (mass in Earth masses, water mass fraction, radius in Earth radii, depth km, base pressure GPa): C 2.4's ocean-world table,
# CMF 0.325, 300 K isotherm, 1 bar at the top. The x = 0.2 rows were recomputed from this module (DOC.79); the table itself is
# read by the last tests in this file, so the document and these rows cannot drift apart.
REFERENCE = (
    (1, 1e-4, 0.9895, 1.2, 0.012), (1, 1e-3, 0.9907, 11.7, 0.119), (1, 1e-2, 1.0004, 102.0, 1.18),
    (1, 0.1, 1.0583, 746.0, 11.1), (1, 0.2, 1.1281, 1294.0, 19.2), (1, 0.5, 1.2301, 3005.0, 45.4),
    (5, 1e-3, 1.5428, 22.8, 0.51), (5, 1e-2, 1.5571, 173.0, 5.03), (5, 0.1, 1.6419, 1189.0, 48.4), (5, 0.2, 1.7463, 2024.0, 84.7),
    (5, 0.5, 1.8942, 4623.0, 205.0),
)


@pytest.mark.parametrize("mass, water, radius, depth_km, base_gpa", REFERENCE)
def test_the_depth_and_the_base_pressure_match_the_compressed_integration(mass, water, radius, depth_km, base_gpa):
    column = wl.column(mass * EARTH, radius * R_EARTH_KM, water * mass * EARTH, bottom_temperature_k=300.0)
    assert column.depth_km == pytest.approx(depth_km, rel=0.03)
    assert column.base_pressure_pa / 1e9 == pytest.approx(base_gpa, rel=0.05 if water < 1e-3 else 0.03)


def test_a_tenth_of_an_earth_mass_is_far_shallower_than_the_fixed_density_layer():
    mass, radius = EARTH, 1.0583 * R_EARTH_KM
    fixed = 0.1 * mass / (1000.0 * 4 * 3.141592653589793 * (radius * 1000.0) ** 2) / 1000.0
    assert fixed == pytest.approx(1047.0 * (1.0583 ** -2) * (1.0583 ** 2), rel=0.1) or fixed > 1000.0
    assert wl.column(mass, radius, 0.1 * mass, bottom_temperature_k=300.0).depth_km < 0.75 * fixed


@pytest.mark.parametrize("mass, water, radius, depth_km, base_gpa", REFERENCE)
def test_the_integration_has_converged_at_the_default_steps(mass, water, radius, depth_km, base_gpa):
    coarse = wl.column(mass * EARTH, radius * R_EARTH_KM, water * mass * EARTH, bottom_temperature_k=300.0)
    fine = wl.column(mass * EARTH, radius * R_EARTH_KM, water * mass * EARTH, bottom_temperature_k=300.0, steps=400)
    assert coarse.depth_km == pytest.approx(fine.depth_km, rel=5e-3)
    assert coarse.base_pressure_pa == pytest.approx(fine.base_pressure_pa, rel=1e-2)


def test_an_inconsistent_radius_and_water_fraction_stays_finite():
    column = wl.column(EARTH, 1.06 * R_EARTH_KM, 0.9 * EARTH, bottom_temperature_k=300.0)
    assert math.isfinite(column.depth_km) and math.isfinite(column.base_pressure_pa) and column.depth_km > 0.0


def test_the_density_follows_the_pressure_through_the_ice_phases():
    """Liquid at 0.5 GPa, ice VI at 1 GPa and ice VII at 3 GPa (the materials table's 300 K points)."""
    assert bulk_structure.water_density(5.0e8) == pytest.approx(1148.5, rel=0.01)
    assert bulk_structure.water_density(1.0e9) == pytest.approx(1351.1, rel=0.02)
    assert bulk_structure.water_density(3.0e9) == pytest.approx(1537.6, rel=0.01)


def test_the_sealing_water_fraction_is_a_few_thousandths():
    """The liquid stops at the liquidus: a column of about 70 km at Earth gravity (about 6e-3 of the mass) puts ice VI under the
    ocean at 285 K, and 3e-3 does not."""
    shallow = wl.column(EARTH, R_EARTH_KM, 3e-3 * EARTH, bottom_temperature_k=285.0)
    deep = wl.column(EARTH, R_EARTH_KM, 1.2e-2 * EARTH, bottom_temperature_k=285.0)
    assert shallow.hp_ice_km == 0.0 and deep.hp_ice_km > 0.0
    low, high = 3e-3, 1.2e-2
    for _ in range(30):
        mid = 0.5 * (low + high)
        if wl.column(EARTH, R_EARTH_KM, mid * EARTH, bottom_temperature_k=285.0).hp_ice_km > 0.0:
            high = mid
        else:
            low = mid
    assert 4e-3 < high < 9e-3


def test_the_liquidus_is_monotone_and_continues_past_the_last_triple_point():
    values = [wl.liquidus_pressure_pa(t) for t in range(240, 420, 5)]
    assert values == sorted(values)
    assert wl.liquidus_pressure_pa(400.0) > wl.liquidus_pressure_pa(353.8)
    assert wl.liquidus_pressure_pa(200.0) == wl.liquidus_pressure_pa(251.2)


def test_a_lid_heavier_than_the_water_freezes_it_solid_and_a_second_run_is_the_same():
    area = 4 * 3.141592653589793 * (R_EARTH_KM * 1000.0) ** 2
    frozen = wl.column(EARTH, R_EARTH_KM, 1.0e6 * area, shell_km=10.0, bottom_temperature_k=273.15)
    assert (frozen.liquid_km, frozen.hp_ice_km) == (0.0, 0.0)
    assert frozen.ice_shell_km == pytest.approx(1.0e6 / 917.0 / 1000.0)
    again = wl.column(EARTH, R_EARTH_KM, 1.0e6 * area, shell_km=10.0, bottom_temperature_k=273.15)
    assert again == frozen


def test_a_partly_covered_planet_holds_its_water_in_a_deeper_column():
    full = wl.column(EARTH, R_EARTH_KM, 2.34e-4 * EARTH, bottom_temperature_k=285.0)
    covered = wl.column(EARTH, R_EARTH_KM, 2.34e-4 * EARTH, covered=0.5, bottom_temperature_k=285.0)
    assert covered.liquid_km == pytest.approx(2.0 * full.liquid_km, rel=0.03)


DOC = pathlib.Path(__file__).resolve().parents[2] / "docs/design/composition/body-composition-parts/C-volatile-and-giants.md"


def _documented_table():
    rows = []
    for line in DOC.read_text().splitlines():
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if line.startswith("| ") and len(cells) == 7:
            try:
                rows.append(tuple(float(cell) for cell in cells))
            except ValueError:
                continue
    return rows


def test_the_documented_ocean_table_has_the_twelve_rows():
    assert len(_documented_table()) == 12


def test_every_documented_row_is_consistent_with_its_own_radius_and_with_the_code():
    """Depth and base pressure come from the code at the row's own radius; the gravity and the estimate follow from the
    row's mass, water share and radius, so no row can disagree with its neighbours the way the old x = 0.2 rows did."""
    for mass, water, radius, depth_km, base_gpa, gravity, estimate_gpa in _documented_table():
        column = wl.column(mass * EARTH, radius * R_EARTH_KM, water * mass * EARTH, bottom_temperature_k=300.0)
        label = f"M {mass} x {water}"
        assert column.depth_km == pytest.approx(depth_km, rel=0.03), label
        assert column.base_pressure_pa / 1e9 == pytest.approx(base_gpa, rel=0.05 if water < 1e-3 else 0.03), label
        assert gravity == pytest.approx(9.80665 * mass / radius ** 2, rel=0.005), label
        area = 4.0 * math.pi * (radius * R_EARTH_KM * 1000.0) ** 2
        assert estimate_gpa == pytest.approx(gravity * water * mass * EARTH / area / 1e9, rel=0.02, abs=2e-3), label


def test_the_documented_depth_and_pressure_rise_with_the_water_share():
    for mass in (1.0, 5.0):
        rows = [row for row in _documented_table() if row[0] == mass]
        assert [row[3] for row in rows] == sorted(row[3] for row in rows)
        assert [row[4] for row in rows] == sorted(row[4] for row in rows)
