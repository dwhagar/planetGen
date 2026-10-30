# tests/test_belt_page.py

"""UX.240: the data of a belt's information page (`planetgen.physics.belt_page`)."""

import json
import math

import pytest

from planetgen.physics import belt_page as bp
from planetgen.galaxy import object_uid
from planetgen.physics import belt_ledger, unit_registry
from planetgen.util import draw

LUMINOSITY = 3.828e26
AL26 = 5.0e-5
BELT_ID = object_uid.parse_id("FE81000A2B1234-0000005-003")
INNER, OUTER = 2.1, 3.3


def _record(belt_mass_kg=4.0e21, seed=5, law=belt_ledger.ROCKY_LAW, planets=((5.2, 1.9e27),)):
    ledger = belt_ledger.draw_ledger(law, belt_mass_kg, draw.Stream(seed))
    n = len(ledger.masses_kg)
    placements = tuple(bp.Placement(INNER + (OUTER - INNER) * (i + 0.5) / max(n, 1), (i * 0.37) % 1.0) for i in range(n))
    return bp.BeltRecord(BELT_ID, INNER, OUTER, 1.0, LUMINOSITY, AL26, ledger, placements, planets)


def _quantities(node, path="page"):
    """Every dict with a `value` in the structure, with where it was found."""
    if isinstance(node, dict):
        if "value" in node:
            yield path, node
        else:
            for key, child in node.items():
                yield from _quantities(child, f"{path}.{key}")
    elif isinstance(node, list):
        for i, child in enumerate(node):
            yield from _quantities(child, f"{path}[{i}]")


def test_the_listed_bodies_and_the_bins_add_up_to_the_belts_mass():
    for seed in (1, 5, 9):
        page = bp.build(_record(seed=seed))
        totals = page["totals"]
        belt = totals["belt_mass"]["value"]
        assert totals["listed"]["value"] + totals["bins"]["value"] == pytest.approx(belt, rel=1e-9)
        assert abs(totals["difference"]["value"]) < 1e-9 * belt
        assert sum(b["mass"]["value"] for b in page["bins"]) == pytest.approx(page["bins"][0]["mass"]["value"] + sum(b["mass"]["value"] for b in page["bins"][1:]))
        shares = sum(b["mass_share"]["value"] for b in page["bins"]) + sum(b["mass"]["value"] for b in page["bodies"]) / belt
        assert shares == pytest.approx(1.0, rel=1e-9)


def test_every_number_comes_with_its_unit():
    page = bp.build(_record())
    found = list(_quantities(page))
    assert len(found) > 100
    for where, quantity in found:
        assert set(quantity) == {"value", "unit", "text"}, where
        assert isinstance(quantity["value"], (int, float)) and math.isfinite(quantity["value"]), where
        if quantity["unit"] == "count":
            assert quantity["text"].endswith(" bodies"), where
        elif quantity["unit"] == "frac":
            assert quantity["text"].endswith(" percent"), where
        else:
            assert quantity["text"].endswith(unit_registry.unit(quantity["unit"]).symbol), where
    json.dumps(page)                                               # plain data: the page and the API read the same structure


def test_the_individual_bodies_are_listed_by_designation_heaviest_first():
    record = _record()
    page = bp.build(record)
    bodies = page["bodies"]
    assert len(bodies) == len(record.ledger.masses_kg)
    assert [b["rank"] for b in bodies] == list(range(1, len(bodies) + 1))
    assert [b["mass"]["value"] for b in bodies] == sorted((b["mass"]["value"] for b in bodies), reverse=True)
    for body in bodies[:5] + bodies[-5:]:
        assert body["designation"] == object_uid.designation(BELT_ID, body["rank"])
        assert object_uid.parse_designation(body["designation"]) == (BELT_ID, body["rank"])
    assert bodies[0]["designation"].endswith(".0001")
    assert bodies == bp.build(record)["bodies"]                    # pure: the same record, the same page


def test_a_body_has_a_composition_a_density_a_size_and_a_status_word():
    body = bp.build(_record())["bodies"][0]
    assert body["taxon"] in body["analogue"] or body["analogue"]
    parts = body["composition"]
    assert sum(parts[k]["value"] for k in ("rock", "metal", "water", "ice")) == pytest.approx(1.0)
    assert 500.0 < body["density"]["value"] < 6000.0
    assert body["diameter"]["unit"] == "km" and body["distance"]["unit"] == "au" and INNER <= body["distance"]["value"] <= OUTER
    assert body["status"]["word"] in ("planet", "dwarf planet", "minor body")
    assert body["status"]["roundness"] <= 1.0


def test_a_big_round_body_alone_in_its_zone_is_a_planet_and_a_crowded_one_is_not():
    law = belt_ledger.ROCKY_LAW
    lone = belt_ledger.Ledger((3.0e23,), 3.01e23, 1.0e21, 1, False, law)
    record = bp.BeltRecord(BELT_ID, INNER, OUTER, 1.0, LUMINOSITY, AL26, lone, (bp.Placement(2.7, 0.3),), planets=())
    lone_body = bp.build(record)["bodies"][0]
    assert lone_body["status"]["word"] == "planet" and lone_body["density_held"]       # over the compaction table: held at its top density
    assert not bp.build(_record())["bodies"][0]["density_held"]
    crowded = belt_ledger.Ledger((5.0e21,), 1.0e24, 9.95e23, 1, False, law)
    record = bp.BeltRecord(BELT_ID, INNER, OUTER, 1.0, LUMINOSITY, AL26, crowded, (bp.Placement(2.7, 0.3),), planets=())
    assert bp.build(record)["bodies"][0]["status"]["word"] in ("dwarf planet", "minor body")
    assert bp.build(record)["bodies"][0]["status"]["soter_mu"] < 1.0


def test_the_bins_run_from_the_bottom_to_the_smallest_stored_body():
    record = _record()
    page = bp.build(record)
    bins = page["bins"]
    smallest = page["bodies"][-1]["diameter"]["value"]
    assert bins[0]["from"] is None                                 # the first runs down to the bottom of the table
    assert bins[-1]["to"]["value"] == pytest.approx(smallest, rel=1e-6)
    uppers = [b["to"]["value"] for b in bins]
    assert uppers == sorted(uppers)
    assert all(u in bp.BIN_EDGES_KM or u == pytest.approx(smallest) for u in uppers)
    assert [b["from"]["value"] for b in bins[1:]] == uppers[:-1]
    assert bins[0]["count"]["value"] > bins[-1]["count"]["value"] > 0.0     # many small, few big
    assert all(b["mass"]["value"] > 0.0 for b in bins)


def test_tiers_follow_the_on_demand_floor():
    page = bp.build(_record())
    floor_km = belt_ledger.diameter_km(bp.L2_FLOOR_KG, belt_ledger.ROCKY_LAW.density_kg_m3)
    for b in page["bins"]:
        low = 0.0 if b["from"] is None else b["from"]["value"]
        assert b["tier"] == ("on demand" if low >= floor_km else "total only")
    assert page["bins"][0]["tier"] == "total only" and page["bins"][-1]["tier"] == "on demand"


def test_the_typical_composition_is_the_kernels_over_the_belt():
    page = bp.build(_record())
    composition = page["belt"]["composition"]
    assert sum(c["share"]["value"] for c in composition) > 0.98
    assert [c["share"]["value"] for c in composition] == sorted((c["share"]["value"] for c in composition), reverse=True)
    for b in page["bins"]:
        assert b["composition"] == composition
        assert 800.0 < b["typical_density"]["value"] < 5000.0
    # a belt in the cold outer zone is icier and less dense than one in the warm inner zone
    cold = bp.build(bp.BeltRecord(BELT_ID, 30.0, 40.0, 1.0, LUMINOSITY, AL26, _record().ledger, tuple(bp.Placement(35.0, 0.3) for _ in _record().ledger.masses_kg), ()))
    assert {c["taxon"] for c in cold["belt"]["composition"]} <= {"D", "P", "C", "B", "S", "E", "L", "K", "A", "M", "V"}
    assert cold["belt"]["composition"][0]["taxon"] in ("D", "P")
    assert cold["bins"][0]["typical_density"]["value"] < page["bins"][0]["typical_density"]["value"]


def test_the_density_factor_counts_objects_per_square_au_and_falls_with_size():
    record = _record()
    page = bp.build(record)
    factor = page["density_factor"]
    area = math.pi * (OUTER ** 2 - INNER ** 2)
    assert [f["above"]["value"] for f in factor["by_diameter"]] == list(bp.DENSITY_DIAMETERS_KM)
    objects = [f["objects"]["value"] for f in factor["by_diameter"]]
    assert objects == sorted(objects, reverse=True) and objects[-1] > 0.0
    assert all(f["objects"]["unit"] == "au-2" and f["objects"]["text"].endswith("objects/AU²") for f in factor["by_diameter"])
    assert factor["mass"]["unit"] == "kg_au2"
    assert factor["mass"]["value"] == pytest.approx(record.ledger.belt_mass_kg / area)
    over_100 = sum(1 for b in page["bodies"] if b["diameter"]["value"] >= 100.0)
    assert over_100 > 0 and page["bodies"][-1]["diameter"]["value"] > 100.0
    assert objects[2] == pytest.approx(over_100 / area)           # every body over 100 km is stored here: the count is exact
    assert objects[1] > objects[2]


def test_a_belt_with_no_stored_body_is_all_bins():
    law = belt_ledger.ROCKY_LAW
    ledger = belt_ledger.Ledger((), 1.0e18, 1.0e18, 0, False, law)
    page = bp.build(bp.BeltRecord(BELT_ID, INNER, OUTER, 1.0, LUMINOSITY, AL26, ledger, (), ()))
    assert page["bodies"] == []
    assert page["totals"]["bins"]["value"] == pytest.approx(1.0e18, rel=1e-9)
    assert page["bins"][0]["from"] is None


def test_a_record_that_does_not_hang_together_is_refused():
    record = _record()
    with pytest.raises(bp.BeltPageError):
        bp.build(bp.BeltRecord(BELT_ID, OUTER, INNER, 1.0, LUMINOSITY, AL26, record.ledger, record.placements, ()))
    with pytest.raises(bp.BeltPageError):
        bp.build(bp.BeltRecord(BELT_ID, INNER, OUTER, 1.0, LUMINOSITY, AL26, record.ledger, record.placements[:-1], ()))
    moved = (bp.Placement(50.0, 0.1),) + record.placements[1:]
    with pytest.raises(bp.BeltPageError):
        bp.build(bp.BeltRecord(BELT_ID, INNER, OUTER, 1.0, LUMINOSITY, AL26, record.ledger, moved, ()))
    with pytest.raises(bp.BeltPageError):
        bp.build(bp.BeltRecord(BELT_ID, INNER, OUTER, 0.0, LUMINOSITY, AL26, record.ledger, record.placements, ()))
