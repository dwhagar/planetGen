# tests/test_galaxy_field_job.py

"""
PERF.225: the galaxy's cluster field is built ahead by a queued job, in small resumable parts, and kept as a file in the
state folder (`galaxy/field_job.py`); the web side only reads it.

Held here: the field made from parts is the field made whole; a part is a file written by a rename and a part that is there is
kept; a run that stopped is carried on; a file that cannot be trusted is not read (and a pickle cannot run anything); the
signature moves with the rows and the code; the store reads the file and a visitor is refused (fast) when a big galaxy has none;
the job is on the queue as the first thing of the warm-up's opening phase and runs end to end on a database.
"""

import dataclasses
import json
import os
import pickle

import pytest

from planetgen.admin import cachebuild
from planetgen.db import store
from planetgen.galaxy import component_field, field_job
from planetgen.queue import commands

from tests.test_component_field import E1, EDGE, _cluster


@pytest.fixture(autouse=True)
def state_folder(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    monkeypatch.setattr(field_job, "PART_AGES", 4)
    return tmp_path


def _many(count=11):
    """`count` clusters, each of an age of its own, so that the parts of 4 ages come to three."""
    return tuple(dataclasses.replace(_cluster(2.0e5 + 1.0e3 * i, 5.0 + 0.1 * i, 5, name=f"C{i}", center=(40.0 * i, 0.0, 0.0)),
                                     age_gy=1.0 + 0.37 * i) for i in range(count))


KEY = ("db", 11, 1, 11, 2.2e6, EDGE, E1)


def _make_parts(sig, found):
    ages = field_job.distinct_ages(c.age_gy for c in found)
    return ages, [field_job.compute_part(sig, i, field_job.part_ages(ages, i)) for i in range(field_job.part_count(ages))]


def _same_field(a, b):
    points = [(0.0, 0.0, 0.0), (40.0, 3.0, 1.0), (200.0, 0.0, 0.0), (5000.0, 1.0, 1.0)]
    return (len(a) == len(b) and a.names == b.names and a.shares == b.shares
            and all(a.densities(p) == b.densities(p) for p in points)
            and a._central.tolist() == b._central.tolist() and a._density_reach.tolist() == b._density_reach.tolist())


# --- The parts and the file --------------------------------------------------------------------------------------------------------

def test_the_field_made_from_parts_is_the_field_made_whole():
    found = _many()
    sig = field_job.signature(KEY)
    ages, done = _make_parts(sig, found)
    assert [d["built"] for d in done] == [1, 1, 1] and sum(d["ages"] for d in done) == len(ages) == 11
    made = field_job.assemble(sig, found, EDGE, E1)
    assert _same_field(made, component_field.ComponentField(found, EDGE, E1))
    assert _same_field(field_job.read(sig), made)


def test_the_parts_cover_every_age_once_in_order():
    ages = list(range(1, 12))
    parts = [field_job.part_ages(ages, i) for i in range(field_job.part_count(ages))]
    assert [x for part in parts for x in part] == ages and all(0 < len(part) <= field_job.PART_AGES for part in parts)
    assert field_job.part_count([]) == 1


def test_a_part_is_a_file_and_one_that_is_there_is_kept(state_folder):
    found = _many()
    sig = field_job.signature(KEY)
    ages = field_job.distinct_ages(c.age_gy for c in found)
    first = field_job.compute_part(sig, 0, field_job.part_ages(ages, 0))
    path = field_job.part_path(sig, 0)
    stamp = os.stat(path).st_mtime_ns
    again = field_job.compute_part(sig, 0, field_job.part_ages(ages, 0))
    assert first["built"] == 1 and again["built"] == 0 and os.stat(path).st_mtime_ns == stamp
    assert not [name for name in os.listdir(os.path.dirname(path)) if name.endswith(".tmp")]
    assert json.load(open(path))["ages"] == field_job.part_ages(ages, 0)


def test_a_run_that_stopped_is_carried_on_not_started_again():
    found = _many()
    sig = field_job.signature(KEY)
    ages = field_job.distinct_ages(c.age_gy for c in found)
    field_job.compute_part(sig, 0, field_job.part_ages(ages, 0))
    with pytest.raises(field_job.PartsMissing):
        field_job.assemble(sig, found, EDGE, E1)
    assert not field_job.is_published(sig)
    rest = [field_job.compute_part(sig, i, field_job.part_ages(ages, i)) for i in (1, 2)]
    assert [d["built"] for d in rest] == [1, 1]
    assert field_job.assemble(sig, found, EDGE, E1) is not None and field_job.is_published(sig)


def test_a_part_for_other_ages_or_cut_short_is_worked_out_again():
    found = _many()
    sig = field_job.signature(KEY)
    ages = field_job.distinct_ages(c.age_gy for c in found)
    field_job.compute_part(sig, 0, field_job.part_ages(ages, 0))
    path = field_job.part_path(sig, 0)
    with open(path, "w") as handle:
        handle.write('{"ages": [9.9], "means": [1.0]')                   # cut short
    assert field_job.compute_part(sig, 0, field_job.part_ages(ages, 0))["built"] == 1
    other = field_job.part_ages(ages, 0)[:-1] + [123.0]
    assert field_job.compute_part(sig, 0, other)["built"] == 1            # the ages changed under it


def test_publishing_removes_the_parts_and_keeps_the_newest_two_fields(state_folder):
    found = _many()
    sigs = []
    for n in range(3):
        key = ("db", 11 + n, 1, 11, 2.2e6, EDGE, E1)
        sigs.append(field_job.signature(key))
        _make_parts(sigs[-1], found)
        field_job.assemble(sigs[-1], found, EDGE, E1)
        os.utime(field_job.field_path(sigs[-1]), (1000 + n, 1000 + n))
    names = sorted(os.listdir(field_job.folder()))
    assert names == sorted(f"{s}.field" for s in sigs[1:])
    assert not any(name.endswith(".parts") for name in names)


# --- A file that cannot be trusted --------------------------------------------------------------------------------------------------

def test_a_missing_foreign_or_cut_short_file_is_not_a_field():
    found = _many()
    sig = field_job.signature(KEY)
    assert field_job.read(sig) is None and not field_job.is_published(sig)
    _make_parts(sig, found)
    field_job.assemble(sig, found, EDGE, E1)
    path = field_job.field_path(sig)
    blob = open(path, "rb").read()
    other = field_job.signature(("db", 12, 1, 12, 1.0, EDGE, E1))
    os.makedirs(field_job.folder(), exist_ok=True)
    open(field_job.field_path(other), "wb").write(blob)                      # a file named for another field
    assert field_job.read(other) is None
    open(path, "wb").write(blob[:len(blob) // 2])                            # cut short
    assert field_job.read(sig) is None
    open(path, "wb").write(b"not a header\n" + blob)
    assert field_job.read(sig) is None


def test_a_pickle_that_names_anything_else_runs_nothing(state_folder):
    marker = state_folder / "ran"

    class Evil:
        def __reduce__(self):
            return (os.system, (f"touch {marker}",))

    sig = field_job.signature(KEY)
    os.makedirs(field_job.folder(), exist_ok=True)
    with open(field_job.field_path(sig), "wb") as handle:
        handle.write(field_job._header(sig, 1) + pickle.dumps(Evil()))
    assert field_job.read(sig) is None and not marker.exists()


def test_the_signature_moves_with_the_rows_and_the_code(monkeypatch):
    assert field_job.signature(KEY) == field_job.signature(tuple(KEY))
    assert field_job.signature(KEY) != field_job.signature(KEY[:-1] + (E1 + 0.1,))
    assert field_job.signature(KEY) != field_job.signature(("db", 12) + KEY[2:])
    before = field_job.signature(KEY)
    monkeypatch.setattr(field_job, "code_version", lambda: "another release")
    assert field_job.signature(KEY) != before and len(before) == 32


def test_the_code_version_reads_the_source_of_every_module_that_decides_the_numbers():
    field_job.code_version.cache_clear()
    assert len(field_job.code_version()) == 64
    assert all(__import__("importlib").import_module(name).__file__ for name in field_job.SOURCES)


# --- The store reads the file --------------------------------------------------------------------------------------------------------

class _Answer:
    def __init__(self, row):
        self.row = row

    def fetchone(self):
        return self.row


def _conn(count):
    row = {"db": "d", "n": count, "first_id": 1, "last_id": count, "mass": 1.0e6}

    class Conn:
        def execute(self, sql, params=()):
            return _Answer({"1": 1} if "information_schema" in sql else dict(row))

    return Conn()


@pytest.fixture
def no_memory(monkeypatch):
    monkeypatch.setattr(store, "_FIELD_BY_SIGNATURE", store.collections.OrderedDict())


def _stored(monkeypatch, found):
    reads = []
    monkeypatch.setattr(store, "get_galaxy_components", lambda conn, kind=None: reads.append(1) or list(found))
    return reads


def test_the_store_reads_the_file_and_not_the_rows(no_memory, monkeypatch):
    found = _many(3)
    key = store._component_key(_conn(3), EDGE, E1)
    sig = field_job.signature(key)
    _make_parts(sig, found)
    field_job.assemble(sig, found, EDGE, E1)
    reads = _stored(monkeypatch, found)
    token = cachebuild.visitor.set(True)
    try:
        field = store._component_field(_conn(3), EDGE, E1)
        assert field is not None and not reads and store._component_field(_conn(3), EDGE, E1) is field
    finally:
        cachebuild.visitor.reset(token)


def test_a_visitor_is_refused_at_once_when_a_big_galaxy_has_no_file(no_memory, monkeypatch):
    reads = _stored(monkeypatch, _many(3))
    token = cachebuild.visitor.set(True)
    try:
        with pytest.raises(cachebuild.CacheRebuilding) as refused:
            store._component_field(_conn(field_job.INLINE_COMPONENTS + 1), EDGE, E1)
    finally:
        cachebuild.visitor.reset(token)
    assert refused.value.state["cause"] == "missing" and not reads


def test_a_small_galaxy_is_built_in_the_request_and_a_worker_may_build_a_big_one(no_memory, monkeypatch):
    reads = _stored(monkeypatch, _many(3))
    token = cachebuild.visitor.set(True)
    try:
        assert store._component_field(_conn(3), EDGE, E1) is not None           # about a second: no file needed
    finally:
        cachebuild.visitor.reset(token)
    assert store._component_field(_conn(field_job.INLINE_COMPONENTS + 1), EDGE, E1) is not None   # no visitor: a worker
    assert len(reads) == 2


# --- The job on the queue ----------------------------------------------------------------------------------------------------------

def test_the_job_is_a_warm_command_and_the_first_step_of_the_opening_phase():
    from planetgen.web import warmup
    spec = commands.REGISTRY["warm-cluster-field"]
    assert spec.unit == ("planetgen.galaxy.field_job", "field_unit") and spec.conclude and spec.unit_seconds
    assert commands.KIND_OF_WARM_FAMILY["cluster-field"] == "warm-cluster-field"
    assert warmup.phase_of("cluster-field") == "opening" and list(warmup.FAMILIES)[0] == "cluster-field"
    args = type("Args", (), {"database": "d", "cap": None, "families": None})()
    steps = commands.build_steps("warm", args)
    assert steps[0].key == "cluster-field" and steps[0].kind == "warm-cluster-field"
    assert "cluster-field" in warmup._WARMERS


def test_the_job_runs_end_to_end_on_a_database_and_a_second_run_finds_the_file(mysql_config, monkeypatch):
    from planetgen.galaxy.density import build_galaxy_shape
    import math
    shape = build_galaxy_shape(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
                               bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=EDGE, outer_shell_index=1500, expected_system_count_at_density_1=E1,
                            config=mysql_config, galaxy_seed=bytes(range(16)))
    found = _many()
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, list(found))
        conn.commit()
        key = store.component_field_key(conn)
    finally:
        conn.close()
    sig = field_job.signature(key)
    params = {"database": mysql_config.database}

    plan = field_job.expand_field(params)
    assert plan.plan["state"] == "parts" and plan.plan["parts"] == 3 and len(plan.units) == 3
    assert all(cost == field_job.PART_SECONDS for _args, cost in plan.units)
    results = [field_job.field_unit(args, None) for args, _cost in plan.units]
    assert [r["built"] for r in results] == [1, 1, 1] and all(r["line"].startswith("cluster-field:") for r in results)
    assert not field_job.is_published(sig)
    summary = field_job.conclude_field(params, plan.plan, results, 1.0)
    assert summary == {"components": 11, "parts": 3} and field_job.is_published(sig)

    conn = store.get_connection(mysql_config)
    try:
        stored = tuple(dataclasses.replace(c, id=None) for c in store.get_galaxy_components(conn))
    finally:
        conn.close()
    assert _same_field(field_job.read(sig), component_field.ComponentField(stored, EDGE, E1))

    again = field_job.expand_field(params)
    assert again.plan["state"] == "ready" and len(again.units) == 1
    assert field_job.field_unit(again.units[0][0], None)["cached"] == 1
    assert field_job.conclude_field(params, again.plan, [], 0.0) == {}

    replaced = list(found[:-1])                                           # the rows change while a run is going
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, replaced)
        conn.commit()
    finally:
        conn.close()
    with pytest.raises(field_job.Changed):
        field_job.field_unit(plan.units[0][0], None)
    assert field_job.expand_field(params).plan["state"] == "parts"        # a new signature: a new field to make


def test_a_galaxy_without_components_has_nothing_to_make(mysql_config):
    plan = field_job.expand_field({"database": mysql_config.database})
    assert plan.plan == {"state": "empty"} and "no field" in field_job.field_unit(plan.units[0][0], None)["line"]
    assert field_job.conclude_field({"database": mysql_config.database}, plan.plan, [], 0.0) == {}


def _plan_a_galaxy(mysql_config):
    from planetgen.galaxy.density import build_galaxy_shape
    import math
    shape = build_galaxy_shape(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
                               bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=EDGE, outer_shell_index=1500, expected_system_count_at_density_1=E1,
                            config=mysql_config, galaxy_seed=bytes(range(16)))
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, list(_many()))
        conn.commit()
        return field_job.signature(store.component_field_key(conn))
    finally:
        conn.close()


def test_ensure_file_makes_the_field_once_for_every_worker_and_leaves_a_file_alone(mysql_config):
    sig = _plan_a_galaxy(mysql_config)
    assert not field_job.is_published(sig)
    assert field_job.ensure_file(mysql_config) == (3, 0) and field_job.is_published(sig)       # made here, by the config's host
    assert field_job.ensure_file(mysql_config) is None                                         # already a file: nothing to do
    assert field_job.ensure_file(mysql_config.database) is None                                 # a database name works too


def test_ensure_file_has_nothing_to_make_for_a_galaxy_without_clusters(mysql_config):
    assert field_job.ensure_file(mysql_config) is None


def test_a_stage_whose_workers_read_the_shape_gets_the_field_made_before_it_is_cut(monkeypatch):
    from planetgen.queue import wrapper
    made = []
    monkeypatch.setattr(field_job, "ensure_file", lambda config: made.append(config.database))
    named = {"host": "h", "port": 3306, "user": "u", "database": "galaxy_db"}
    wrapper._ensure_cluster_field({"database": named})
    assert made == ["galaxy_db"]
    wrapper._ensure_cluster_field({})                                                           # a command with no database
    assert made == ["galaxy_db"]
    monkeypatch.setattr(field_job, "ensure_file", lambda config: (_ for _ in ()).throw(RuntimeError("no room")))
    wrapper._ensure_cluster_field({"database": named})                                          # best effort: the stage goes on


def test_only_the_stages_whose_workers_warm_the_shape_ask_for_the_field():
    warming = {kind for kind, spec in commands.REGISTRY.items() if spec.warm is not None}
    assert {"scatter-mass", "scatter-luminosity", "scatter-phenomena"} <= warming
    assert "plan-skeleton" not in warming
