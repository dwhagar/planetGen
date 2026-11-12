"""GEN.1611 slice 0: the snow line a star's disk was built round (`physics/formation_snow_line.py`)."""

import pytest

from planetgen import tuning
from planetgen.physics import formation, formation_snow_line as fsl


def test_a_main_sequence_or_contracting_star_keeps_the_luminosity_it_has():
    for cls in fsl.UNEVOLVED_CLASSES:
        assert fsl.formation_luminosity_sol(1.0, 0.0, 1.1, cls) == 1.1
    assert fsl.formation_snow_line_au(1.0, 0.0, 1.1, "V") == pytest.approx(formation.snow_line_au(1.1 * 3.828e26))


@pytest.mark.parametrize("mass, present, cls, low, high", [
    (4.36, 1.0e-4, "VII", 20.0, 80.0),        # a white dwarf born at 4.36 solar masses: its disk saw a 330 solar-luminosity B star (was 0.02 AU)
    (1.43, 1.0e-4, "VII", 2.0, 8.0),
    (1.04, 60.0, "III", 1.0, 4.0),             # a K giant born at 1.04 solar masses (was 21.8 AU)
    (1.73, 120.0, "III", 4.0, 12.0),
])
def test_an_evolved_star_gets_the_snow_line_of_its_birth_not_of_today(mass, present, cls, low, high):
    line = fsl.formation_snow_line_au(mass, 0.0, present, cls)
    assert low < line < high
    assert line != pytest.approx(formation.snow_line_au(present * 3.828e26), rel=0.5)


def test_it_draws_nothing_and_reads_the_disk_age():
    from planetgen.util import draw
    draw.set_run_seed(7)
    first = draw.random()
    draw.set_run_seed(7)
    value = fsl.formation_luminosity_sol(4.36, 0.0, 1.0e-4, "VII")
    assert draw.random() == first
    assert value == fsl.formation_luminosity_sol(4.36, 0.0, 1.0e-4, "VII")
    assert tuning.FORMATION_DISK_AGE_GY == 0.01
