"""
GEN.1608: a comet's period class agrees with its orbit across the 200-year line that `comet_designation` splits periodic (`P/`) from the rest (`C/`).
"""

import pytest

from planetgen.generation import comet as comet_module
from planetgen.generation.comet import Comet, PERIODIC_COMET_MAX_PERIOD_YEARS, comet_designation, period_class_of
from planetgen.generation.config import SystemConfig
from planetgen.util import draw


@pytest.mark.parametrize("drawn, period, expected", [
    ("long_period", 30.0, "halley_type"),            # the ceiling cut its eccentricity: tens of years is no long-period comet
    ("long_period", 199.9, "halley_type"),
    ("long_period", 200.0, "long_period"),
    ("long_period", 5000.0, "long_period"),
    ("halley_type", 350.0, "long_period"),
    ("jupiter_family", 260.0, "long_period"),
    ("halley_type", 90.0, "halley_type"),
    ("jupiter_family", 12.0, "jupiter_family"),
    ("jupiter_family", 60.0, "jupiter_family"),      # a Jupiter-family comet is placed by its giant, and keeps its family under the line
])
def test_the_class_follows_the_orbit_across_the_periodic_line(drawn, period, expected):
    assert period_class_of(drawn, period) == expected


@pytest.mark.parametrize("ceiling_au, mass_solar", [(None, 1.0), (60.0, 0.3), (20.0, 0.2)])
def test_a_drawn_comet_is_long_period_exactly_when_it_is_not_periodic(ceiling_au, mass_solar):
    config = SystemConfig()
    for index in range(300):
        with draw.bound(draw.Stream(f"comet-class:{ceiling_au}:{index}")):
            comet = Comet(config, mass_solar, orbit_type="elliptical", ceiling_au=ceiling_au)
        periodic = comet.orbital_period_years < PERIODIC_COMET_MAX_PERIOD_YEARS
        assert (comet.period_class == "long_period") != periodic
        assert comet_designation("Star", 1, comet).startswith("P/" if periodic else "C/")
