"""
GEN.1765 slice 2: a star's comets read two shapes from its outer disks where it has them (`physics/comet_supply`): the family shares from the disk's scatter step, and the top of the nucleus draw
from its largest body. A star with no scattered outer disk keeps the table's weights and the range's top.
"""

import types

import pytest

from planetgen import tuning
from planetgen.generation.comet import Comet
from planetgen.generation.config import SystemConfig
from planetgen.physics import comet_supply
from planetgen.util import draw


def _belt(origin="outer-disk", fates=None, mass=10.0, largest=300.0):
    return types.SimpleNamespace(origin=origin, scatter_fates=fates, primordial_mass_earth=mass, largest_body_km=largest, body_type="a")


def test_no_outer_disk_or_no_scattered_one_keeps_the_defaults():
    assert comet_supply.family_weights([]) is None
    assert comet_supply.nucleus_top_km([]) is None
    assert comet_supply.family_weights([_belt(origin="inner", fates={"scattered_disk": 0.5})]) is None
    assert comet_supply.family_weights([_belt(fates=None)]) is None
    assert comet_supply.family_weights([_belt(fates={"cold_classical": 1.0})]) is None


def test_the_inward_populations_make_the_jupiter_family_and_the_detached_the_cloud():
    weights = comet_supply.family_weights([_belt(fates={"scattered_disk": 0.3, "centaurs": 0.1, "detached": 0.1, "cold_classical": 0.5})])
    assert weights["jupiter_family"] == pytest.approx(10.0 * 0.4)
    assert weights["halley_type"] + weights["long_period"] == pytest.approx(10.0 * 0.1)
    classes = tuning.COMET_PERIOD_CLASSES
    assert weights["halley_type"] / weights["long_period"] == pytest.approx(classes["halley_type"]["weight"] / classes["long_period"]["weight"])


def test_disks_add_by_the_mass_they_had():
    heavy = _belt(fates={"scattered_disk": 1.0}, mass=9.0)
    light = _belt(fates={"detached": 1.0}, mass=1.0)
    weights = comet_supply.family_weights([heavy, light])
    assert weights["jupiter_family"] == pytest.approx(9.0)
    assert weights["halley_type"] + weights["long_period"] == pytest.approx(1.0)


def test_the_nucleus_top_is_the_largest_body_of_the_outer_disks():
    assert comet_supply.nucleus_top_km([_belt(largest=8.0), _belt(largest=12.0), _belt(origin="inner", largest=500.0)]) == 12.0


def _comets(count, **supply):
    config = SystemConfig()
    out = []
    for index in range(count):
        with draw.bound(draw.Stream(f"comet-supply:{index}")):
            out.append(Comet(config, 1.0, orbit_type="elliptical", **supply))
    return out


def test_a_top_under_the_range_cuts_the_nuclei_and_one_under_its_lower_end_gives_the_lower_end():
    low, high = tuning.BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM
    assert max(c.nucleus_diameter_km for c in _comets(200, nucleus_top_km=6.0)) <= 6.0
    assert {c.nucleus_diameter_km for c in _comets(20, nucleus_top_km=0.5)} == {low}
    assert max(c.nucleus_diameter_km for c in _comets(200, nucleus_top_km=1000.0)) <= high


def test_the_weights_move_the_family_shares():
    comets = _comets(400, family_weights={"jupiter_family": 0.0, "halley_type": 1.0, "long_period": 0.0})
    drawn = [c for c in comets]
    assert sum(1 for c in drawn if c.period_class == "halley_type") > 0.9 * len(drawn)         # a comet keeps the family its orbit gives (GEN.1612), so nearly all
    base = _comets(400)
    assert sum(1 for c in base if c.period_class == "jupiter_family") > sum(1 for c in drawn if c.period_class == "jupiter_family") + 50


def test_a_disk_that_scattered_only_inward_still_draws_a_family_when_the_giant_holds_none():
    comets = _comets(30, family_weights={"jupiter_family": 1.0, "halley_type": 0.0, "long_period": 0.0}, giant_au=5.2)
    assert len(comets) == 30
