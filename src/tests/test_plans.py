# tests/test_plans.py

"""
The change plan engine (ADM.58, `planetgen/admin/plans.py`): every kind of operation, the effect diff, the states and
who made each move, the stale check, the drift an orbital update causes, and the envelope. Pure Python over generated
systems; no database.
"""

import copy
import json

import pytest

from planetgen.admin import edits, plans
from planetgen.admin.plans import (
    Add, Move, Plan, PlanError, PlanStale, PlanStateError, Ref, Regenerate, Remove, SetAspect, World)
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

ADMIN = ("boss", "admin")
EDGE = 100.0


def make_system(label):
    """A single G2V system with at least three planets, a belt and a planet with moons, every body with a row id."""
    for attempt in range(200):
        with draw.bound(draw.Stream(f"plans-{label}-{attempt}")):
            config = SystemConfig()
            config.STAR_TYPE, config.PLANETS, config.MOONS = "G2V", True, True
            config.ASTEROID_BELT, config.BINARY_SYSTEM = True, False
            system = StarSystem(config)
        worlds = [p for p in system.planets if p.body_type != 'a']
        if len(worlds) >= 3 and any(p.body_type == 'a' for p in system.planets) and any(p.moons for p in worlds):
            break
    else:
        pytest.skip("no suitable system generated")
    system.star.db_id = 1
    for index, body in enumerate(system.planets):
        body.db_id = (200 if body.body_type == 'a' else 10) + index
        for moon_index, moon in enumerate(getattr(body, "moons", ())):
            moon.db_id = 100 + index * 10 + moon_index
    return system


def stamp_all(world):
    world.stamps = {key: 1 for system_id, system in world.systems.items() for key in plans._index(system_id, system)}
    return world


_GENERATED = {}


@pytest.fixture
def world():
    """Two generated systems (made once for the file, copied for each test, since a test may move a body)."""
    for label in ("one", "two"):
        if label not in _GENERATED:
            _GENERATED[label] = make_system(label)
    return stamp_all(World({"s1": copy.deepcopy(_GENERATED["one"]), "s2": copy.deepcopy(_GENERATED["two"])},
                           {"s1": (10.0, 10.0, 10.0), "s2": (60.0, 60.0, 60.0)}, EDGE))


def planets_of(world, system_id="s1"):
    return [p for p in world.systems[system_id].planets if p.body_type != 'a']


def belt_of(world, system_id="s1"):
    return next(p for p in world.systems[system_id].planets if p.body_type == 'a')


def ref(world, body, system_id="s1"):
    if body.body_type == 'a':
        return Ref(system_id, plans.BELT, body.db_id)
    return Ref(system_id, plans.MOON if getattr(body, "is_moon", False) else plans.PLANET, body.db_id)


def planned(world, *operations, seed=7):
    plan = Plan(seed=seed, title="test")
    for op in operations:
        plan.add(op)
    plan.propose(world, *ADMIN)
    return plan


def fields_of(plan, key):
    return {c.field: c for c in plan.effect.changes if c.ref.key == key}


# --- references and the envelope ---------------------------------------------------------------------------------

def test_a_reference_names_one_object_and_a_system_has_no_id_of_its_own():
    assert Ref("s1").key == "s1" and Ref("s1", plans.PLANET, 12).key == "s1/planet:12"
    assert Ref.from_dict(Ref("s1", plans.MOON, 4).to_dict()) == Ref("s1", plans.MOON, 4)
    with pytest.raises(PlanError):
        Ref("s1", plans.PLANET)
    with pytest.raises(PlanError):
        Ref("s1", plans.SYSTEM, 3)
    with pytest.raises(PlanError):
        Ref("s1", "galaxy", 3)


def test_every_operation_round_trips_through_its_envelope():
    refs = Ref("s1", plans.PLANET, 12)
    for op in (Move(refs, distance_au=1.5, before={"distance_au": 1.0}), Move(Ref("s1"), position=(1.0, 2.0, 3.0)),
               SetAspect(refs, "name", "Vesper"), Regenerate(refs, keep=("name", "orbit")),
               Add(Ref("s1"), plans.BELT, 4.0, lower_au=3.0, upper_au=5.0, new_id="new-1"), Remove(refs, forced=True)):
        data = json.loads(json.dumps(op.to_dict()))
        assert data["op"] == op.OP and data["schema_version"] == plans.ENVELOPE_VERSION
        assert plans.Operation.from_dict(data) == op
    with pytest.raises(PlanError, match="unknown operation"):
        plans.Operation.from_dict({"op": "teleport"})
    with pytest.raises(PlanError, match="newer"):
        plans.Operation.from_dict({"op": "remove", "target": refs.to_dict(), "schema_version": 999})


# --- move --------------------------------------------------------------------------------------------------------

def test_moving_a_planet_changes_its_distance_directly_and_leaves_the_stored_system_alone(world):
    planet = planets_of(world)[0]
    old = planet.distance
    plan = planned(world, Move(ref(world, planet), distance_au=old * 1.05))
    changes = fields_of(plan, f"s1/planet:{planet.db_id}")
    assert changes["distance"].before == pytest.approx(old) and changes["distance"].after == pytest.approx(old * 1.05)
    assert not changes["distance"].indirect
    assert "Moved " in plan.effect.summaries[0] and plan.effect.direct == [f"s1/planet:{planet.db_id}"]
    assert planet.distance == old, "the plan works on copies"
    assert plan.outcomes["s1"]["planets"] and list(plan.outcomes) == ["s1"]


def test_a_move_that_crowds_a_neighbour_moves_the_neighbour_and_says_so(world):
    inner, next_out = planets_of(world)[0], planets_of(world)[1]
    plan = planned(world, Move(ref(world, inner), distance_au=next_out.distance * 0.999))
    moved = {c.ref.key: c for c in plan.effect.changes if c.field == "distance"}
    assert not moved[f"s1/planet:{inner.db_id}"].indirect
    pushed = moved[f"s1/planet:{next_out.db_id}"]
    assert pushed.indirect and pushed.after > pushed.before, "the spacing pass pushed the neighbour out"
    assert plan.effect.indirect == [pushed.ref.key] or pushed.ref.key in plan.effect.indirect


def test_a_move_just_past_a_neighbour_is_spaced_out_and_the_summary_says_where_it_ended(world):
    inner, next_out = planets_of(world)[0], planets_of(world)[1]
    plan = planned(world, Move(ref(world, inner), distance_au=next_out.distance * 1.001))
    final = fields_of(plan, f"s1/planet:{inner.db_id}")["distance"].after
    assert final > next_out.distance * 1.001 and "to keep the orbits spaced" in plan.effect.summaries[0]


def test_a_planet_cannot_be_moved_past_the_farthest_stable_orbit_or_to_nothing(world):
    planet = ref(world, planets_of(world)[0])
    for distance, message in ((1.0e9, "farthest stable orbit"), (-1.0, "positive"), ("far", "must be a number")):
        plan = Plan()
        plan.add(Move(planet, distance_au=distance))
        with pytest.raises(PlanError, match=message):
            plan.propose(world, *ADMIN)
        assert plan.state == plans.DRAFT and plan.effect is None, "a refused plan stays as it was"


def test_a_belt_moves_with_its_limits(world):
    belt = belt_of(world)
    width = belt.upper_limit - belt.lower_limit
    old = belt.distance
    plan = planned(world, Move(ref(world, belt), distance_au=old * 1.02))
    changes = fields_of(plan, f"s1/belt:{belt.db_id}")
    assert changes["distance"].after == pytest.approx(old * 1.02)
    assert changes["upper_limit"].after - changes["lower_limit"].after == pytest.approx(width)


def test_a_moon_moves_inside_the_orbits_its_planet_can_hold(world):
    planet = next(p for p in planets_of(world) if p.moons)
    moon = planet.moons[0]
    low = plans._km_to_au(plans.planetPhysics.moon_orbit_bounds_km(planet)[0])
    high = plans._km_to_au(plans.planetPhysics.moon_orbit_bounds_km(planet)[1])
    target = Ref("s1", plans.MOON, moon.db_id)
    plan = planned(world, Move(target, distance_au=(low + high) / 2))
    assert fields_of(plan, target.key)["distance"].after == pytest.approx((low + high) / 2)
    outside = Plan()
    outside.add(Move(target, distance_au=high * 3))
    with pytest.raises(PlanError, match="can sit"):
        outside.propose(world, *ADMIN)
    forced = Plan()
    forced.add(Move(target, distance_au=high * 3, forced=True))
    forced.propose(world, *ADMIN)
    assert forced.state == plans.PROPOSED


def test_a_system_moves_inside_its_sector_and_clear_of_its_neighbours(world):
    plan = planned(world, Move(Ref("s1"), position=(30.0, 20.0, 10.0)))
    change = fields_of(plan, "s1")["position"]
    assert change.before == (10.0, 10.0, 10.0) and change.after == (30.0, 20.0, 10.0) and not change.indirect
    assert world.placements["s1"] == (10.0, 10.0, 10.0)
    for position, message in (((EDGE + 1, 0.0, 0.0), "outside the sector"), ((60.0, 60.0, 60.5), "Hill")):
        bad = Plan()
        bad.add(Move(Ref("s1"), position=position))
        with pytest.raises(PlanError, match=message):
            bad.propose(world, *ADMIN)
    crowded = Plan()
    crowded.add(Move(Ref("s1"), position=(60.0, 60.0, 60.5), forced=True))
    crowded.propose(world, *ADMIN)
    with pytest.raises(PlanError, match="three coordinates"):
        planned(world, Move(Ref("s1")))


def test_other_kinds_of_object_cannot_be_moved(world):
    with pytest.raises(PlanError, match="cannot be moved"):
        planned(world, Move(Ref("s1", plans.STAR, 1), distance_au=1.0))


# --- set an aspect -----------------------------------------------------------------------------------------------

def test_renaming_a_planet_renames_its_moons_as_indirect_changes(world):
    planet = next(p for p in planets_of(world) if p.moons)
    plan = planned(world, SetAspect(ref(world, planet), "name", "Vesper"))
    assert fields_of(plan, f"s1/planet:{planet.db_id}")["name"].after == "Vesper"
    moon_changes = [c for c in plan.effect.changes if c.ref.kind == plans.MOON and c.field == "name"]
    assert moon_changes and all(c.indirect and c.after.startswith("Vesper") for c in moon_changes)
    assert plan.effect.indirect


def test_a_moon_the_admin_named_by_hand_keeps_its_name_when_the_planet_is_renamed(world):
    planet = next(p for p in planets_of(world) if p.moons)
    planet.moons[0].name = "Hand-named"
    plan = planned(world, SetAspect(ref(world, planet), "name", "Vesper"))
    assert not [c for c in plan.effect.changes if c.ref.key == f"s1/moon:{planet.moons[0].db_id}"]


def test_a_system_can_be_renamed_and_a_blank_name_is_refused(world):
    plan = planned(world, SetAspect(Ref("s1"), "name", "Hearth"))
    assert fields_of(plan, "s1")["name"].after == "Hearth"
    with pytest.raises(PlanError, match="1 to 120"):
        planned(world, SetAspect(Ref("s1"), "name", "  "))


def test_a_class_is_derived_from_the_body_so_a_plan_cannot_set_it(world):
    """ADM.137: a body's class follows from its conditions (epoch 121), so the class aspect is gone: it is refused as unknown, with
    or without the force."""
    target = ref(world, planets_of(world)[0])
    for forced in (False, True):
        with pytest.raises(PlanError, match="unknown aspect"):
            planned(world, SetAspect(target, "planet_class", "M", forced=forced))


def test_changing_the_star_changes_the_planets_as_indirect_effects(world):
    star = Ref("s1", plans.STAR, 1)
    plan = planned(world, SetAspect(star, "star_type", "K2V"))
    assert fields_of(plan, star.key)["luminosity"].after < fields_of(plan, star.key)["luminosity"].before
    assert any(c.indirect and c.ref.kind == plans.PLANET and c.field == "distance" for c in plan.effect.changes)
    with pytest.raises(PlanError, match="not a spectral type"):
        planned(world, SetAspect(star, "star_type", "banana"))


def test_an_unknown_aspect_or_one_aimed_at_the_wrong_kind_is_refused(world):
    with pytest.raises(PlanError, match="unknown aspect"):
        planned(world, SetAspect(ref(world, planets_of(world)[0]), "mood", "calm"))
    with pytest.raises(PlanError, match="cannot be set"):
        planned(world, SetAspect(ref(world, belt_of(world)), "star_type", "K2V"))


def test_a_registered_aspect_joins_the_list(world, monkeypatch):
    monkeypatch.setitem(plans.ASPECTS, "stub", None)
    plans.register_aspect("stub", (plans.PLANET,), lambda sim, op, system, obj, owner:
                          setattr(obj, "name", "Stubbed") or edits.EditResult("Stubbed.", [], [], [], []))
    planet = planets_of(world)[0]
    plan = planned(world, SetAspect(ref(world, planet), "stub", 1))
    assert fields_of(plan, f"s1/planet:{planet.db_id}")["name"].after == "Stubbed"


# --- regenerate --------------------------------------------------------------------------------------------------

def test_regenerating_a_planet_keeps_its_name_and_its_orbit(world):
    planet = planets_of(world)[0]
    plan = planned(world, Regenerate(ref(world, planet)))
    changes = fields_of(plan, f"s1/planet:{planet.db_id}")
    assert "name" not in changes and changes, "something was redrawn and the name stayed"
    assert "Regenerated" in plan.effect.summaries[0]
    sim, _results = plans.simulate(world, plan.operations, plan.seed)
    redrawn, _owner = edits.find_body(sim.systems["s1"], "planet", planet.db_id)
    assert redrawn.name == planet.name and redrawn.distance == planet.distance


def test_a_regeneration_cannot_keep_the_class_which_is_derived(world):
    """ADM.137: the class is derived from the redrawn body's conditions, so it is not an aspect a regeneration keeps."""
    with pytest.raises(PlanError, match="cannot keep class"):
        planned(world, Regenerate(ref(world, planets_of(world)[0]), keep=("name", "orbit", "class")))


def test_a_belt_and_a_moon_can_be_regenerated(world):
    planet = next(p for p in planets_of(world) if p.moons)
    assert planned(world, Regenerate(ref(world, belt_of(world)))).state == plans.PROPOSED
    assert planned(world, Regenerate(Ref("s1", plans.MOON, planet.moons[0].db_id))).state == plans.PROPOSED


def test_a_regeneration_that_asks_for_what_the_engine_cannot_do_is_refused(world):
    target = ref(world, planets_of(world)[0])
    for keep, message in ((("name",), "keeps the name and the orbit"), (("name", "orbit", "mood"), "cannot keep mood")):
        with pytest.raises(PlanError, match=message):
            planned(world, Regenerate(target, keep=keep))
    with pytest.raises(PlanError, match="cannot be regenerated"):
        planned(world, Regenerate(Ref("s1")))


def test_the_same_plan_on_the_same_stored_state_is_the_same_plan(world):
    planet = planets_of(world)[0]
    operations = [Regenerate(ref(world, planet))]
    first = planned(world, *operations, seed=3)
    again = Plan(seed=3)
    again.add(Regenerate(ref(world, planet)))
    again.propose(world, *ADMIN)
    assert first.effect.to_dict() == again.effect.to_dict() and first.outcomes == again.outcomes
    others = []
    for seed in range(4, 12):
        other = Plan(seed=seed)
        other.add(Regenerate(ref(world, planet)))
        other.propose(world, *ADMIN)
        others.append(other.effect.to_dict() != first.effect.to_dict())
    assert any(others), "another seed draws another planet"


# --- add and remove ----------------------------------------------------------------------------------------------

def free_orbit(world):
    """An AU distance beyond the outermost body and under the farthest stable orbit."""
    system = world.systems["s1"]
    last = max(system.planets, key=plans.validation.edge_au)
    ceiling = plans.validation.orbit_ceiling_au(system.star)
    return min((plans.validation.edge_au(last) + ceiling) / 2, plans.validation.edge_au(last) * 3)


def test_adding_a_planet_gives_it_a_provisional_id_that_later_operations_may_name(world):
    plan = Plan(seed=5)
    add = plan.add(Add(Ref("s1"), plans.PLANET, free_orbit(world), name="Newcomer"))
    assert add.new_id == "new-1"
    plan.add(SetAspect(Ref("s1", plans.PLANET, "new-1"), "name", "Latecomer"))
    plan.propose(world, *ADMIN)
    added = fields_of(plan, "s1/planet:new-1")
    assert added["*"].before is None and added["*"].after["name"] == "Latecomer" and not added["*"].indirect
    assert not any(planet.name in ("Newcomer", "Latecomer") for planet in planets_of(world)), "the plan works on copies"
    assert "s1/planet:new-1" in plan.effect.direct and "s1/planet:new-1" not in plan.base_stamps


def test_adding_a_planet_needs_a_name_and_an_orbit_the_system_can_hold(world):
    with pytest.raises(PlanError, match="needs a name"):
        planned(world, Add(Ref("s1"), plans.PLANET, free_orbit(world)))
    with pytest.raises(PlanError, match="farthest stable orbit"):
        planned(world, Add(Ref("s1"), plans.PLANET, 1.0e9, name="Far"))
    with pytest.raises(PlanError, match="added to a system"):
        planned(world, Add(Ref("s1", plans.PLANET, planets_of(world)[0].db_id), plans.PLANET, 1.0, name="Odd"))
    with pytest.raises(PlanError, match="provisional id"):
        plans.simulate(world, [Add(Ref("s1"), plans.PLANET, 1.0, name="Loose")], 1)


def test_adding_a_moon_to_a_planet_names_it_after_the_planet(world):
    planet = planets_of(world)[0]
    low, high = (plans._km_to_au(v) for v in plans.planetPhysics.moon_orbit_bounds_km(planet))
    plan = planned(world, Add(ref(world, planet), plans.MOON, (low + high) / 2))
    added = [c for c in plan.effect.changes if c.ref.kind == plans.MOON and c.field == "*" and c.after]
    assert added and added[0].after["name"].startswith(planet.name) and not added[0].indirect
    with pytest.raises(PlanError, match="can sit"):
        planned(world, Add(ref(world, planet), plans.MOON, high * 5))
    with pytest.raises(PlanError, match="added to a planet"):
        planned(world, Add(Ref("s1"), plans.MOON, low))


def test_adding_a_belt_needs_its_limits(world):
    distance = free_orbit(world)
    plan = planned(world, Add(Ref("s1"), plans.BELT, distance, lower_au=distance * 0.95, upper_au=distance * 1.05))
    assert fields_of(plan, "s1/belt:new-1")["*"].after["distance"] == pytest.approx(distance)
    with pytest.raises(PlanError, match="lower and upper"):
        planned(world, Add(Ref("s1"), plans.BELT, distance, lower_au=distance, upper_au=distance * 0.9))
    with pytest.raises(PlanError, match="lower limit"):
        planned(world, Add(Ref("s1"), plans.BELT, distance))
    with pytest.raises(PlanError, match="cannot add a comet"):
        planned(world, Add(Ref("s1"), plans.COMET, distance))


def test_removing_a_planet_removes_its_moons_with_it_and_a_moon_alone_leaves_the_planet(world):
    planet = next(p for p in planets_of(world) if p.moons)
    plan = planned(world, Remove(ref(world, planet)))
    removed = {c.ref.key for c in plan.effect.changes if c.field == "*" and c.after is None}
    assert f"s1/planet:{planet.db_id}" in removed and {f"s1/moon:{m.db_id}" for m in planet.moons} <= removed
    assert not fields_of(plan, f"s1/planet:{planet.db_id}")["*"].indirect
    assert all(c.indirect for c in plan.effect.changes if c.ref.kind == plans.MOON and c.field == "*")
    moon = planet.moons[0]
    alone = planned(world, Remove(Ref("s1", plans.MOON, moon.db_id)))
    assert {c.ref.key for c in alone.effect.changes if c.field == "*"} == {f"s1/moon:{moon.db_id}"}
    with pytest.raises(PlanError, match="no such planet"):
        planned(world, Remove(Ref("s1", plans.PLANET, 99999)))
    with pytest.raises(PlanError, match="no such system"):
        planned(world, Remove(Ref("nowhere", plans.PLANET, 1)))
    with pytest.raises(PlanError, match="cannot be removed"):
        planned(world, Remove(Ref("s1")))


def test_operations_run_in_order_over_two_systems(world):
    first, second = planets_of(world)[0], planets_of(world, "s2")[0]
    plan = planned(world, SetAspect(ref(world, first), "name", "One"),
                   Move(Ref("s2"), position=(70.0, 70.0, 70.0)),
                   SetAspect(ref(world, second, "s2"), "name", "Two"))
    assert set(plan.outcomes) == {"s1", "s2"} and len(plan.effect.summaries) == 3
    assert {c.ref.system for c in plan.effect.changes} == {"s1", "s2"}
    with pytest.raises(PlanError, match="operation 2"):
        planned(world, SetAspect(ref(world, first), "name", "One"), Remove(Ref("s1", plans.PLANET, 99999)))


def test_a_plan_with_no_operations_cannot_be_proposed_or_requested(world):
    with pytest.raises(PlanError, match="at least one"):
        Plan().propose(world, *ADMIN)
    with pytest.raises(PlanError, match="at least one"):
        Plan().request("ana", plans.ROLE_USER)


def test_layout_problems_the_plan_leaves_are_warnings_and_the_ones_already_there_are_not(world, monkeypatch):
    planet = planets_of(world)[0]
    Problem = plans.validation.Problem
    monkeypatch.setattr(plans.validation, "check_star_system", lambda system: [Problem("Old", "was already wrong")])
    plan = planned(world, SetAspect(ref(world, planet), "name", "Vesper"))
    assert plan.effect.warnings == []
    calls = iter([[Problem("Old", "was already wrong")], [Problem("Old", "was already wrong"), Problem("New", "is wrong now")]])
    monkeypatch.setattr(plans.validation, "check_star_system", lambda system: next(calls))
    plan = Plan()
    plan.add(SetAspect(ref(world, planet), "name", "Vesper"))
    plan.propose(world, *ADMIN)
    assert plan.effect.warnings == ["New: is wrong now."]


# --- the states ---------------------------------------------------------------------------------------------------

def test_an_admins_plan_goes_from_draft_to_proposed_to_approved_with_no_second_approver(world):
    plan = Plan(title="Rename")
    plan.add(SetAspect(Ref("s1"), "name", "Hearth"))
    assert plan.state == plans.DRAFT
    plan.propose(world, "boss", "admin", now=10.0)
    plan.approve(world, "boss", "admin", now=20.0)
    plan.mark_applied("boss", "admin", now=30.0)
    plan.mark_undone("boss", "admin", now=40.0)
    assert plan.state == plans.UNDONE and plan.author == "boss" and plan.author_role == "admin"
    assert [(t.from_state, t.to_state, t.who, t.role, t.at) for t in plan.history] == [
        ("draft", "proposed", "boss", "admin", 10.0), ("proposed", "approved", "boss", "admin", 20.0),
        ("approved", "applied", "boss", "admin", 30.0), ("applied", "undone", "boss", "admin", 40.0)]


def test_a_request_from_below_the_approving_role_is_proposed_and_approved_by_an_admin(world):
    plan = Plan()
    plan.add(SetAspect(Ref("s1"), "name", "Hearth"))
    plan.request("ana", plans.ROLE_USER, now=1.0)
    assert plan.state == plans.REQUESTED and plan.author == "ana" and plan.author_role == "user"
    with pytest.raises(PlanStateError, match="only admin can propose"):
        plan.propose(world, "ana", plans.ROLE_USER)
    plan.propose(world, "boss", "admin", now=2.0)
    assert plan.author == "ana", "the requester stays the author"
    with pytest.raises(PlanStateError, match="only admin can approve"):
        plan.approve(world, "ana", plans.ROLE_USER)
    plan.approve(world, "boss", "admin", now=3.0)
    assert [(t.who, t.role) for t in plan.history] == [("ana", "user"), ("boss", "admin"), ("boss", "admin")]


def test_a_plan_cannot_take_a_move_the_table_does_not_allow(world):
    plan = planned(world, SetAspect(Ref("s1"), "name", "Hearth"))
    for call in (lambda: plan.mark_applied("boss", "admin"), lambda: plan.mark_undone("boss", "admin"),
                 lambda: plan.request("boss", "admin"), lambda: plan.propose(world, "boss", "admin")):
        with pytest.raises(PlanStateError):
            call()
    plan.reject("boss", "admin", note="no")
    assert plan.state == plans.REJECTED and plan.history[-1].note == "no"
    for call in (lambda: plan.approve(world, "boss", "admin"), lambda: plan.withdraw("boss", "admin"),
                 lambda: plan.add(Remove(Ref("s1", plans.PLANET, 1)))):
        with pytest.raises(PlanStateError):
            call()
    assert all(not plans.TRANSITIONS[state] for state in (plans.UNDONE, plans.REJECTED, plans.WITHDRAWN))
    assert set(plans.TRANSITIONS) == set(plans.STATES)
    with pytest.raises(PlanStateError, match="unknown role"):
        odd = Plan()
        odd.add(Remove(Ref("s1", plans.PLANET, 1)))
        odd.request("x", "emperor")


def test_a_plan_is_withdrawn_by_its_author_or_an_admin(world):
    plan = Plan()
    plan.add(SetAspect(Ref("s1"), "name", "Hearth"))
    plan.request("ana", plans.ROLE_USER)
    with pytest.raises(PlanStateError, match="author or an admin"):
        plan.withdraw("eve", plans.ROLE_USER)
    plan.withdraw("ana", plans.ROLE_USER, note="changed my mind")
    assert plan.state == plans.WITHDRAWN
    other = planned(world, SetAspect(Ref("s1"), "name", "Hearth"))
    other.approve(world, "boss", "admin")
    other.withdraw("boss", "admin")
    assert other.state == plans.WITHDRAWN


# --- stale -------------------------------------------------------------------------------------------------------

def test_a_plan_is_stale_when_an_object_it_touches_was_edited_after_it_was_proposed(world):
    planet = planets_of(world)[0]
    plan = planned(world, SetAspect(ref(world, planet), "name", "Vesper"))
    assert plan.base_stamps[f"s1/planet:{planet.db_id}"] == 1
    world.stamps[f"s1/planet:{planets_of(world)[1].db_id}"] = 2        # another planet was edited: no matter
    assert plan.check_stale(world) == [] and plan.state == plans.PROPOSED
    world.stamps[f"s1/planet:{planet.db_id}"] = 2
    with pytest.raises(PlanStale, match=f"planet:{planet.db_id}"):
        plan.approve(world, "boss", "admin")
    assert plan.state == plans.STALE and plan.history[-1].role == "system"
    assert plan.check_stale(world) == [], "a stale plan is not marked again"
    plan.propose(world, "boss", "admin")
    assert plan.state == plans.PROPOSED and plan.base_stamps[f"s1/planet:{planet.db_id}"] == 2
    plan.approve(world, "boss", "admin")


def test_an_approved_plan_goes_stale_too_and_an_object_added_by_it_has_no_stamp(world):
    plan = planned(world, SetAspect(Ref("s1"), "name", "Hearth"))
    plan.approve(world, "boss", "admin")
    world.stamps["s1"] = 5
    assert plan.check_stale(world) == ["s1"] and plan.state == plans.STALE
    assert plan.stale_keys(world) == ["s1"]


def test_the_orbital_update_moves_what_a_move_starts_from_without_making_the_plan_stale(world):
    planet = planets_of(world)[1]
    start = planet.distance
    plan = planned(world, Move(ref(world, planet), distance_au=start * 1.03))
    assert plan.operations[0].before == {"distance_au": start} and plan.drift(world) == []
    planet.distance = start * 1.01                      # the nightly update moved it; no edit stamp changed
    assert plan.check_stale(world) == [] and plan.state == plans.PROPOSED
    drift = plan.drift(world)
    assert [(d.ref, d.field, d.planned, d.current) for d in drift] == [
        (ref(world, planet), "distance_au", start, start * 1.01)]
    with pytest.raises(plans.PlanError, match="orbital update moved"):
        plan.approve(world, "boss", "admin")
    rebased = plan.rebase(world, "boss", "admin")
    assert rebased == drift and plan.operations[0].before == {"distance_au": start * 1.01}
    assert fields_of(plan, f"s1/planet:{planet.db_id}")["distance"].before == pytest.approx(start * 1.01)
    plan.approve(world, "boss", "admin")


def test_rebasing_an_approved_plan_sends_it_back_to_proposed(world):
    planet = planets_of(world)[1]
    plan = planned(world, Move(ref(world, planet), distance_au=planet.distance * 1.03))
    plan.approve(world, "boss", "admin")
    planet.distance *= 1.01
    plan.rebase(world, "boss", "admin", now=9.0)
    assert plan.state == plans.PROPOSED and plan.history[-1].note == "rebased after the orbital update"
    plan.approve(world, "boss", "admin")
    applied = Plan()
    applied.add(Move(Ref("s1"), position=(20.0, 20.0, 20.0)))
    applied.propose(world, "boss", "admin")
    applied.approve(world, "boss", "admin")
    applied.mark_applied("boss", "admin")
    with pytest.raises(PlanStateError, match="cannot be rebased"):
        applied.rebase(world, "boss", "admin")


def test_a_system_move_drifts_when_the_system_was_placed_somewhere_else(world):
    plan = planned(world, Move(Ref("s1"), position=(30.0, 20.0, 10.0)))
    world.placements["s1"] = (11.0, 10.0, 10.0)
    drift = plan.drift(world)
    assert [(d.field, d.planned, d.current) for d in drift] == [("position", (10.0, 10.0, 10.0), (11.0, 10.0, 10.0))]


# --- the envelope ------------------------------------------------------------------------------------------------

def test_a_plan_survives_a_trip_through_json_with_its_effect_and_outcomes(world):
    planet = next(p for p in planets_of(world) if p.moons)
    plan = Plan(seed=11, title="Round trip")
    plan.add(SetAspect(ref(world, planet), "name", "Vesper"))
    plan.add(Add(Ref("s1"), plans.PLANET, free_orbit(world), name="Newcomer"))
    plan.add(Move(Ref("s2"), position=(70.0, 70.0, 70.0)))
    plan.propose(world, "boss", "admin", now=5.0)
    plan.approve(world, "boss", "admin", now=6.0)
    data = json.loads(json.dumps(plan.to_dict()))
    assert data["plan_version"] == plans.PLAN_VERSION and data["schema_version"] == plans.ENVELOPE_VERSION
    assert data["state"] == "approved" and len(data["operations"]) == 3
    again = Plan.from_dict(data)
    assert json.loads(json.dumps(again.to_dict())) == data
    assert again.state == plans.APPROVED and again.operations == plan.operations
    assert [c.ref.key for c in again.effect.changes] == [c.ref.key for c in plan.effect.changes]
    shown = StarSystem.from_dict(again.outcomes["s1"])
    names = {getattr(p, "name", None) for p in shown.planets}
    assert {"Newcomer", "Vesper"} <= names


def test_a_plan_of_a_newer_version_or_an_unknown_state_is_refused():
    with pytest.raises(PlanError, match="newer"):
        Plan.from_dict({"plan_version": 99, "state": "draft", "id": "x", "operations": []})
    with pytest.raises(PlanError, match="unknown plan state"):
        Plan.from_dict({"state": "dreaming", "id": "x", "operations": []})
