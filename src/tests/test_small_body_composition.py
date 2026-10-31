"""GEN.442: the composition of belts, comets, icy bodies and rings (`planetgen/physics/small_body_composition.py`;
docs/design/small-bodies-composition.md sections 3 to 5 and 9). The Solar System bodies named below are test targets only: the
draw must be able to reach them, and nothing in the module is centred on one."""

import concurrent.futures
import math
import multiprocessing
import statistics

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.physics import compaction, condensate_shares
from planetgen.physics import small_body_composition as sb
from planetgen.util import draw


SUN_W = 3.828e26
WATER_LINE_AU = condensate_shares.ice_line_au(SUN_W)        # where a blackbody round the Sun is as cold as the water line of the derived kernel
AL26 = 5.2e-5


def _belts(x, count, seed, al26=AL26):
    """Belts at `x` times the water line of a Sun-like star, round a star born with `al26`."""
    with draw.bound(seed):
        return [sb.draw_belt(x * WATER_LINE_AU, SUN_W, al26) for _ in range(count)]


def _comets(count, seed, **kwargs):
    with draw.bound(seed):
        return [sb.draw_comet(**kwargs) for _ in range(count)]


# --- the tables --------------------------------------------------------------------------------------------------

def test_the_taxon_table_is_registered_gated_by_kind_and_every_row_is_a_stand_in():
    table = pt.REGISTRY.get("small-body-taxa")
    assert [row.name for row in table.rows] == [f"taxon-{t}" for t in sb.TAXA]
    assert all(row.state == pt.STAND_IN and row.replacement == sb.REPLACEMENT for row in table.rows)
    assert all(set(row.classes) == {sb.BELT, sb.FIELD} for row in table.rows)
    pt.REGISTRY.validate()
    assert sb.REPLACEMENT in pt.REGISTRY.report()
    for kind in (sb.BELT, sb.FIELD):
        total = math.fsum(p for _row, p in table.probabilities(pt.Subject(kind, {"temperature_k": 150.0, "al26_ratio": AL26})))
        assert total == pytest.approx(1.0, abs=1e-12)
    with pytest.raises(pt.PossibilityError):
        table.probabilities(pt.Subject("planet", {"temperature_k": 150.0, "al26_ratio": AL26}))


def test_the_taxa_are_mixes_of_the_shared_table_with_grain_densities_near_their_analogues():
    from planetgen.physics import materials
    for taxon in sb.TAXA:
        data = sb.TAXON_DATA[taxon]
        assert math.fsum(data.fractions.values()) == pytest.approx(1.0, abs=1e-9), taxon
        for name in data.fractions:
            assert materials.has(name), (taxon, name)
    grain = {t: compaction.grain_density_kg_m3(sb.TAXON_DATA[t].fractions) for t in sb.TAXA}
    # Britt and Consolmagno 2003 grain densities (part D 2.2): CI 2.26, CM 2.71, ordinary chondrites 3.48 to 3.64
    assert 2400 <= grain["C"] <= 3100 and 3300 <= grain["S"] <= 3800
    # bulk density at the middle of the family's porosity range against Carry 2012 Table 3 (part D 2.3): C 1.57, S 2.66, M 3.85+,
    # P 1.4 to 1.5, D 1.0 to 1.1 g/cm3 (the Carry values are rho_inf means with large scatter)
    def mid(taxon):
        low, high = sb.POROSITY_RANGES[sb.TAXON_DATA[taxon].family]
        return grain[taxon] * (1.0 - 0.5 * (low + high))
    assert mid("C") == pytest.approx(1570, abs=350) and mid("S") == pytest.approx(2660, abs=450)
    assert mid("M") == pytest.approx(4300, abs=700)
    assert mid("P") == pytest.approx(1400, abs=350) and mid("D") == pytest.approx(1050, abs=350)


# --- a belt ------------------------------------------------------------------------------------------------------

def test_a_belt_has_shares_that_add_up_and_a_density_that_follows_them():
    for belt in _belts(1.0, 40, seed=3):
        shares = belt.detail["taxon_shares"]
        assert set(shares) == set(sb.TAXA) and math.fsum(shares.values()) == pytest.approx(1.0, abs=1e-9)
        assert math.fsum(belt.components.values()) == pytest.approx(1.0, abs=1e-9)
        assert math.fsum(belt.groups().values()) == pytest.approx(1.0, abs=1e-9)
        harmonic = 1.0 / math.fsum(shares[t] / b for t, b in belt.detail["taxon_bulk_density_kg_m3"].items())
        assert belt.bulk_density_kg_m3 == pytest.approx(harmonic, rel=1e-9)
        assert min(belt.detail["taxon_bulk_density_kg_m3"].values()) <= belt.bulk_density_kg_m3 <= max(
            belt.detail["taxon_bulk_density_kg_m3"].values())
        assert 0.0 < belt.porosity < 0.6 and 0.04 < belt.albedo < 0.54
        assert belt.bulk_density_kg_m3 < belt.grain_density_kg_m3


def test_a_belt_moves_with_its_distance_over_the_water_line_in_the_physical_direction():
    water, density = {}, {}
    xs = (0.7, 0.85, 1.0, 1.25, 1.5, 1.95)
    for x in xs:
        belts = _belts(x, 150, seed=11)
        water[x] = statistics.mean(b.water for b in belts)
        density[x] = statistics.mean(b.bulk_density_kg_m3 for b in belts)
    assert water[0.7] < 0.02 and water[1.95] > 0.1                              # dry and igneous inside, ice-rich beyond
    ordered = [water[x] for x in xs]
    assert all(b >= a - 0.005 for a, b in zip(ordered, ordered[1:])), water
    assert density[1.95] < density[0.85] and density[1.0] < density[0.85]       # the ice-rich outer belt is lighter than the dry inner one
    dry = ("S", "A", "K", "L", "M", "V")
    inner = _belts(0.85, 150, seed=5)
    assert statistics.mean(sum(b.detail["taxon_shares"][t] for t in dry) for b in inner) > 0.6
    middle = _belts(1.0, 150, seed=5)
    assert statistics.mean(b.detail["taxon_shares"]["C"] for b in middle) > 0.3
    outer = _belts(1.95, 150, seed=5)
    assert statistics.mean(b.detail["taxon_shares"]["D"] + b.detail["taxon_shares"]["P"] for b in outer) > 0.6


def test_the_aluminium_26_of_the_star_sets_how_igneous_the_dry_belt_is():
    """GEN.473: a star born without aluminium-26 has an undifferentiated dry belt; one born with the Solar System's has a core, mantle and crust share."""
    igneous = ("M", "A", "V")
    cold = statistics.mean(sum(b.detail["taxon_shares"][t] for t in igneous) for b in _belts(0.7, 100, seed=3, al26=1e-7))
    hot = statistics.mean(sum(b.detail["taxon_shares"][t] for t in igneous) for b in _belts(0.7, 100, seed=3, al26=2e-4))
    assert cold < 0.02 and hot > 0.4


def test_the_disk_temperature_and_not_the_distance_sets_the_zone():
    """A brighter star puts the same zone farther out: the belt is the same at the same disk temperature, and the draw says which temperature."""
    near = sb.draw_belt(WATER_LINE_AU, SUN_W, AL26)
    bright = sb.draw_belt(2.0 * WATER_LINE_AU, 4.0 * SUN_W, AL26)
    assert bright.detail["temperature_k"] == pytest.approx(near.detail["temperature_k"])
    assert near.detail["temperature_k"] == pytest.approx(condensate_shares.ICE_CONDENSATION_K["H2O"], rel=1e-3)
    assert near.detail["al26_ratio"] == AL26
    for bad in ((0.0, SUN_W, AL26), (2.0, 0.0, AL26), (2.0, SUN_W, -1e-5)):
        with pytest.raises(ValueError):
            sb.draw_belt(*bad)


# --- a comet -----------------------------------------------------------------------------------------------------

def test_a_comet_has_shares_that_add_up_and_its_ice_is_one_over_one_plus_the_refractory_ratio():
    for comet in _comets(60, seed=4, period_class="jupiter_family"):
        assert math.fsum(comet.components.values()) == pytest.approx(1.0, abs=1e-9)
        assert math.fsum(comet.groups().values()) == pytest.approx(1.0, abs=1e-9)
        delta = comet.detail["refractory_ice_ratio"]
        assert 0.5 <= delta <= 8.0
        assert comet.water + comet.ice == pytest.approx(1.0 / (1.0 + delta), rel=1e-9)
        assert comet.rock == pytest.approx(1.0 - 1.0 / (1.0 + delta), rel=1e-9) and comet.metal == 0.0
        assert 0.5 <= comet.porosity <= 0.8 and 0.0 <= comet.detail["organic_share"] <= 1.0
        assert comet.bulk_density_kg_m3 == pytest.approx(comet.grain_density_kg_m3 * (1.0 - comet.porosity))
        co, co2 = comet.detail["volatiles"]["CO"], comet.detail["volatiles"]["CO2"]
        assert comet.components["co_ice"] / comet.components["water_ice"] == pytest.approx(
            co / 100.0 * 28.0101 / 18.01528, rel=1e-3)
        assert co > 0 and co2 > 0


def test_a_comet_density_reaches_a_measured_nucleus_and_the_mass_follows_the_size():
    comets = _comets(2000, seed=7, period_class="jupiter_family")
    density = sorted(c.bulk_density_kg_m3 for c in comets)
    assert 250 < statistics.median(density) < 900
    # a 67P-like nucleus (537.8 kg/m3, Choukroun 2020; 1.6 km mean radius) is reachable: at least one draw in 14 is within 100
    assert sum(abs(d - 537.8) < 100.0 for d in density) / len(density) > 0.07
    comet = comets[0]
    assert comet.mass_kg(3.4) == pytest.approx(4.0 / 3.0 * math.pi * (1700.0) ** 3 * comet.bulk_density_kg_m3)
    assert comet.mass_kg(6.8) == pytest.approx(8.0 * comet.mass_kg(3.4))


def test_the_ice_of_a_longer_period_comet_carries_more_carbon_monoxide():
    short = _comets(600, seed=9, period_class="jupiter_family")
    long = _comets(600, seed=9, period_class="long_period")
    assert statistics.median(c.detail["volatiles"]["CO"] for c in long) > 2 * statistics.median(c.detail["volatiles"]["CO"] for c in short)
    assert _comets(3, seed=1, period_class=None)[0].detail["volatiles"]
    with draw.bound(5):
        interstellar = sb.draw_comet("long_period", interstellar=True)
    assert interstellar.detail["volatiles"]["CO"] > 0.0


# --- icy bodies and rings ----------------------------------------------------------------------------------------

def test_an_icy_body_draws_its_rock_fraction_about_seventy_percent_and_its_density_from_the_compaction_law():
    with draw.bound(21):
        bodies = [sb.draw_icy_body(sb.KUIPER_BELT_OBJECT) for _ in range(300)]
    rocks = [b.rock for b in bodies]
    assert all(0.45 <= r <= 0.95 for r in rocks) and statistics.median(rocks) == pytest.approx(0.70, abs=0.02)
    assert {b.detail["row"] for b in bodies} == {"icy-intact"}
    for body in bodies[:20]:
        assert body.metal == 0.0 and body.rock + body.water == pytest.approx(1.0) and body.bulk_density_kg_m3 is None
        assert body.density_at(300.0) == compaction.compacted_density_kg_m3(300.0, body.rock)
        assert 700 < body.density_at(300.0) < 2800
    with pytest.raises(ValueError):
        bodies[0].mass_kg(10.0)


def test_a_dwarf_planet_may_be_a_stripped_one_and_a_class_d_body_may_not():
    with draw.bound(22):
        dwarfs = [sb.draw_icy_body(sb.DWARF_PLANET) for _ in range(400)]
        class_d = [sb.draw_icy_body(sb.CLASS_D) for _ in range(100)]
    assert {b.detail["row"] for b in dwarfs} == {"icy-intact", "icy-stripped"}
    stripped = [b.rock for b in dwarfs if b.detail["row"] == "icy-stripped"]
    assert statistics.median(stripped) == pytest.approx(0.88, abs=0.03)
    assert {b.detail["row"] for b in class_d} == {"icy-intact"}
    # design note F3: a class D body whose rock fraction is under 0.9 is lighter than 3000 kg/m3
    assert all(b.density_at(250.0) < 3000.0 for b in class_d if b.rock < 0.9)
    with pytest.raises(pt.PossibilityError):
        sb.draw_icy_body("planet")


def test_a_ring_gains_ice_with_the_distance_past_the_ice_line():
    # the order of the Galilean ramp of decision S8 (Io, Europa, Ganymede, Callisto at r / r_ice 0.295, 0.47, 0.75, 1.3)
    fractions = [sb.disk_ice_fraction(r, 0.35) for r in (0.295, 0.47, 0.75, 1.315)]
    assert fractions[0] == 0.0
    assert fractions[1] == pytest.approx(0.0618, abs=0.002) and fractions[2] == pytest.approx(0.240, abs=0.003)
    assert fractions[3] == pytest.approx(0.35)
    assert fractions == sorted(fractions)
    with draw.bound(31):
        rings = [(r, sb.draw_ring(r)) for r in (0.3, 0.6, 1.0, 2.0, 4.0) for _ in range(30)]
    for r, ring in rings:
        assert ring.rock + ring.water == pytest.approx(1.0) and ring.metal == 0.0
        assert ring.water == pytest.approx(sb.disk_ice_fraction(r, ring.detail["f_inf"]))
    means = [statistics.mean(ring.water for r, ring in rings if r == x) for x in (0.3, 0.6, 1.0, 2.0, 4.0)]
    assert means[0] == 0.0 and means[:3] == sorted(means[:3]) and means[1] > 0.0 and means[2] > 0.25
    assert all(ring.water <= ring.detail["f_inf"] + 1e-12 for _r, ring in rings)       # the ramp saturates at f_inf past 1 r_ice
    with pytest.raises(ValueError):
        sb.disk_ice_fraction(0.0, 0.3)


# --- the groups and determinism ----------------------------------------------------------------------------------

def test_the_groups_split_metal_water_other_ice_and_rock():
    assert [sb.group_of(n) for n in ("kamacite", "water_ice", "co2_ice", "formaldehyde", "forsterite", "complex_organics", "iron",
                                     "troilite", "schreibersite")] == ["metal", "water", "ice", "ice", "rock", "rock", "metal", "rock",
                                                                       "metal"]


def test_a_makeup_converts_to_lane_is_bulk_fractions_when_they_exist():
    pytest.importorskip("planetgen.physics.bulk_composition")
    belt = _belts(1.0, 1, seed=8)[0]
    bulk = belt.bulk_fractions()
    assert bulk.porosity == belt.porosity and bulk.envelope == 0.0 and bulk.water == belt.water


def _worker(seed):
    with draw.bound(seed):
        belt = sb.draw_belt(3.1, SUN_W, AL26)
        comet = sb.draw_comet("halley_type")
        icy = sb.draw_icy_body(sb.DWARF_PLANET)
        ring = sb.draw_ring(1.4)
    return (belt.detail["taxon_shares"], belt.bulk_density_kg_m3, comet.components, comet.bulk_density_kg_m3, icy.rock,
            icy.detail["row"], ring.water)


def test_the_same_seed_gives_the_same_bodies_at_one_and_four_workers():
    seeds = list(range(8))
    serial = [_worker(s) for s in seeds]
    assert serial == [_worker(s) for s in seeds] and serial[0] != serial[1]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4, mp_context=multiprocessing.get_context("spawn")) as pool:
        assert list(pool.map(_worker, seeds)) == serial


def test_a_comets_organics_take_the_shared_tables_density_not_a_separate_one():
    """GEN.712: complex organics are 1.5 g/cm3 in the materials table (the middle of 1.0 to 2.0, Blum and Bischoff 2022) and a comet
    states no other density for them; only the silicate dust, which the table holds none for, carries a stand-in."""
    from planetgen.physics import materials
    row = materials.material("complex_organics")
    prop = row.props["rho_grain_g_cm3"]
    assert prop.value == 1.5
    assert tuple(prop.measured_range) == (1.0, 2.0)
    # GEN.660: the 2 g/cm3 that Choukroun et al. 2020 gives (section 3.3, read) is the top of the table's range, and the row says so.
    assert prop.measured_range[1] >= 2.0 and "Choukroun" in row.note
    comet = sb.draw_comet()
    assert "complex_organics" not in comet.densities
    assert sb.REFRACTORY_DENSITIES == {"silicate_dust": 3500.0}
    organics_only = compaction.grain_density_kg_m3({"complex_organics": 1.0})
    assert organics_only == pytest.approx(1500.0)


def test_a_belt_names_the_chondrite_group_its_mean_grains_are_nearest():
    """GEN.475: the group is written into the draw's detail, and it is a group of the map or `none`; a dry inner belt is never named a carbonaceous group's ice-rich kin."""
    from planetgen.physics import meteorite_groups
    allowed = {*meteorite_groups.CHONDRITE_GROUPS, meteorite_groups.NONE}
    dry = _belts(0.5, 40, 5)
    wet = _belts(3.0, 40, 6)
    assert all(b.detail["group"] in allowed and b.detail["group_distance"] >= 0.0 for b in dry + wet)
    assert len({b.detail["group"] for b in dry}) > 1, "the group follows the draw, it is not one fixed word"
    assert {b.detail["group"] for b in wet} == {meteorite_groups.NONE}, "an icy belt is far from every chondrite group"


def test_the_group_axes_of_a_pure_mix_are_the_volume_shares_of_its_rows():
    pure_metal = sb.group_axes({"kamacite": 0.6, "taenite": 0.4})
    assert pure_metal[0] is None and pure_metal[1] == pytest.approx(100.0) and pure_metal[2:] == (0.0, 0.0)
    olivine = sb.group_axes({"forsterite": 0.8, "fayalite": 0.2, "serpentine": 0.0})
    assert olivine[0] is not None and 5.0 < olivine[0] < 30.0 and olivine[1:] == (0.0, 0.0, 0.0)


def test_the_group_costs_no_draw():
    with draw.bound(11):
        first = sb.draw_belt(2.0, SUN_W, AL26)
        after_first = draw.uniform(0.0, 1.0)
    with draw.bound(11):
        sb.draw_belt(2.0, SUN_W, AL26)
        assert draw.uniform(0.0, 1.0) == after_first
    assert first.detail["group"]
