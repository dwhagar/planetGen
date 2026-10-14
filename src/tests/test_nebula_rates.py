"""
The nebula rates stay inside what the Milky Way has (TEST.136, GEN.152, GEN.203, GEN.204): the counts and fillings
the tuning produces are held to the bands of docs/design/nebula-density-vs-reality.md, the field keeps its survivors'
positions and object IDs, and an O star hosts an H II region about one time in four.
"""

import collections
import math
import types

import pytest

from planetgen import tuning
from planetgen.galaxy import nebula_field
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.generation import run_sector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw

SHAPE = build_galaxy_shape(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
                           bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4,
                           calibration_radius_pc=None)
DEFAULT_GALAXY_STARS = 2.953e11
"""float: Stars in the default galaxy (the design doc's figure), for the scatter kinds' counts."""
ARM_CREST_GAS_FACTOR = 2.3
"""float: The gas factor at an arm crest on the solar circle."""
SEED = bytes(range(16))


def _field_total():
    """The galaxy's drawn field clouds: density times volume times gas factor, over a coarse grid."""
    rho = tuning.PHENOMENON_DENSITY_PC3["molecular-cloud"] * tuning.PHENOMENON_RATE_SCALE.get("molecular-cloud", 1.0)
    step_xy, step_z, radius, zmax = 1000.0, 100.0, 15000.0, 600.0
    total = 0.0
    x = -radius + step_xy / 2
    while x < radius:
        y = -radius + step_xy / 2
        while y < radius:
            if math.hypot(x, y) < radius:
                z = -zmax + step_z / 2
                while z < zmax:
                    total += rho * step_xy * step_xy * step_z * nebula_field.gas_factor((x, y, z), SHAPE)
                    z += step_z
            y += step_xy
        x += step_xy
    return total


def _class_shares():
    weights = {letter: tuning.NEBULA_CLASSES[letter]["frequency"] for letter in nebula_field.FIELD_CLASSES}
    whole = sum(weights.values())
    return {letter: weight / whole for letter, weight in weights.items()}


def _mean_sphere_volume_pc3(draws=600):
    """Mean bounding-sphere volume of a drawn cloud by class, from a fixed run of draws."""
    draw.set_run_seed(36)
    volumes = collections.defaultdict(list)
    for _ in range(draws):
        nebula = Nebula(SystemConfig(), nebula_type=nebula_field.FIELD_FAMILY)
        volumes[nebula.nebula_class].append(4 / 3 * math.pi * ly_to_pc(nebula.radius_ly) ** 3)
    return {letter: sum(values) / len(values) for letter, values in volumes.items()}


def test_the_nebula_counts_and_fillings_stay_inside_the_observed_bands():
    keep = tuning.NEBULA_FIELD_KEEP
    shares = _class_shares()
    total = _field_total()
    giant = total * shares["M"] * keep["M"]
    assert 5_000 <= giant <= 20_000, f"class M: {giant:,.0f} (8,107 to 9,710 catalogued)"
    small = total * (shares["P"] * keep["P"] + shares["Q"] * keep["Q"])
    assert 5e4 <= small <= 4e5, f"classes P and Q: {small:,.0f}"

    volumes = _mean_sphere_volume_pc3()
    rho = tuning.PHENOMENON_DENSITY_PC3["molecular-cloud"] * tuning.PHENOMENON_RATE_SCALE.get("molecular-cloud", 1.0)
    per_pc3 = rho * ARM_CREST_GAS_FACTOR * sum(shares[letter] * keep[letter] * volumes[letter] for letter in shares)
    crest = 1.0 - math.exp(-per_pc3)
    assert 0.01 <= crest <= 0.04, f"arm-crest filling {crest:.1%} (observed 1 to 2% in the arms)"
    assert 0.01 <= tuning.GMC_ARM_FILLING_FACTOR <= 0.02

    per_star = tuning.phenomenon_rate_per_star
    nebulae = per_star("planetary-nebula") * DEFAULT_GALAXY_STARS
    assert 4_000 <= nebulae <= 46_000, f"planetary nebulae: {nebulae:,.0f}"
    remnants = per_star("supernova-remnant") * DEFAULT_GALAXY_STARS
    assert 400 <= remnants <= 5_000, f"supernova remnants: {remnants:,.0f} (310 catalogued)"
    hosted = [chance for letter, _low, _high, _classes, chance in tuning.NEBULA_HOST_RULES if letter == "O"]
    assert hosted and all(0.1 <= chance <= 0.5 for chance in hosted), hosted


def test_the_field_keeps_the_drawn_share_of_each_class_and_the_survivors_do_not_move():
    drawn = collections.Counter()
    kept = collections.Counter()
    for i in range(20, 50):
        for j in range(-8, 8):
            index = (i, j, 0)
            full = nebula_field.cell_draw(SEED, SHAPE, index)
            survivors = nebula_field.cell_clouds(SEED, SHAPE, index)
            for n, nebula, _center in full:
                drawn[nebula.nebula_class] += 1
            for nebula, _center in survivors:
                kept[nebula.nebula_class] += 1
            by_place = {center: nebula.nebula_class for _n, nebula, center in full}
            assert all(by_place[center] == nebula.nebula_class for nebula, center in survivors)
    assert drawn["M"] > 100
    for letter in "MN":
        share = kept[letter] / drawn[letter]
        assert abs(share - tuning.NEBULA_FIELD_KEEP[letter]) < 0.35 * tuning.NEBULA_FIELD_KEEP[letter] + 0.02, \
            (letter, share)
    for letter in "PQ":
        assert kept[letter] == drawn[letter]


def test_a_surviving_cloud_keeps_its_object_id(monkeypatch):
    edge_pc = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
    survivors = []
    for i in range(38, 46):
        index = (i, 0, 0)
        survivors += [(nebula, center) for nebula, center in nebula_field.cell_clouds(SEED, SHAPE, index)]
        if len(survivors) >= 4:
            break
    survivors = survivors[:4]
    assert len(survivors) == 4
    under_the_cut = [nebula_field.cloud_object_id(SEED, SHAPE, nebula, center, edge_pc) for nebula, center in survivors]
    monkeypatch.setattr(tuning, "NEBULA_FIELD_KEEP", {"M": 1.0, "N": 1.0, "P": 1.0, "Q": 1.0})
    # The same clouds, found again in the uncut draw, hold the same IDs: the cut counts the clouds it drops.
    uncut = {center: nebula for i in range(38, 46) for nebula, center in nebula_field.cell_clouds(SEED, SHAPE, (i, 0, 0))}
    assert all(center in uncut for _nebula, center in survivors)
    without_the_cut = [nebula_field.cloud_object_id(SEED, SHAPE, uncut[center], center, edge_pc)
                       for _nebula, center in survivors]
    assert under_the_cut == without_the_cut


@pytest.mark.parametrize("letter,expected", [("O", 0.25), ("B", 0.1)])
def test_an_o_or_early_b_star_hosts_a_nebula_at_its_rule_chance(monkeypatch, letter, expected):
    class Star:
        type = f"{letter}1V"

    class Entry:
        def __init__(self):
            self.star_system = types.SimpleNamespace(stars=[Star()], name="x")
            self.position = (0.0, 0.0, 0.0)

    class Sector:
        name = "s"

        def __init__(self, count):
            self.entries = [Entry() for _ in range(count)]
            self.added = 0

        def add_phenomenon(self, nebula, kind, position=None):
            self.added += 1
            return nebula

    monkeypatch.setattr(run_sector, "Nebula", lambda config, nebula_class=None: types.SimpleNamespace(name="n"))
    monkeypatch.setattr(run_sector, "choose_weighted_class", lambda classes, label: classes[0])
    draw.set_run_seed(204)
    sector = Sector(4000)
    run_sector.add_star_hosted_nebulae(sector, types.SimpleNamespace(markdown=False))
    assert abs(sector.added / 4000 - expected) < 0.03, sector.added / 4000
