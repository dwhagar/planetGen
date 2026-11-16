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
    ("jupiter_family", 60.0, "halley_type"),         # GEN.1612: the family is the one whose period range holds the period, whatever it was drawn as
    ("jupiter_family", 2.0, "jupiter_family"),       # shorter than the shortest range: the shortest family
    ("long_period", 5.0e6, "long_period"),           # longer than the longest range: the longest family
    ("halley_type", 20.0, "halley_type"),            # a range takes its lower end
    ("jupiter_family", 20.0, "halley_type"),
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


@pytest.mark.parametrize("ceiling_au, mass_solar", [(None, 1.0), (None, 0.3), (5000.0, 1.0), (5000.0, 0.3)])
def test_a_drawn_comets_period_sits_in_its_familys_range(ceiling_au, mass_solar):
    """GEN.1612: the orbit is drawn within the family's period range (a perihelion that cannot give it is drawn again), so no comet of a star that leaves room has a period outside its family's."""
    from planetgen import tuning
    config = SystemConfig()
    outside = []
    for index in range(400):
        with draw.bound(draw.Stream(f"comet-range:{ceiling_au}:{mass_solar}:{index}")):
            comet = Comet(config, mass_solar, orbit_type="elliptical", ceiling_au=ceiling_au)
        low, high = tuning.COMET_PERIOD_CLASSES[comet.period_class]["period_range_years"]
        if not low <= comet.orbital_period_years < high:
            outside.append((comet.period_class, comet.orbital_period_years))
    assert not outside, outside[:5]


def test_the_families_keep_their_shares_of_the_elliptical_comets():
    """GEN.1612: drawing inside the ranges leaves the family weights as the table gives them (0.5, 0.3, 0.2), where relabelling the strays had left the long-period family at 0.14."""
    from collections import Counter
    from planetgen import tuning
    config = SystemConfig()
    drawn = Counter()
    for index in range(600):
        with draw.bound(draw.Stream(f"comet-share:{index}")):
            drawn[Comet(config, 1.0, orbit_type="elliptical").period_class] += 1
    for name, data in tuning.COMET_PERIOD_CLASSES.items():
        assert drawn[name] / 600 == pytest.approx(data["weight"], abs=0.06), (name, drawn)
