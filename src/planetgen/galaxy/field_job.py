# planetgen/galaxy/field_job.py

"""
The galaxy's cluster field as a file, built ahead by a queued job in small resumable parts (PERF.225).

`galaxy.component_field.ComponentField` is built from the `galaxy_components` rows, and a galaxy with tens of thousands of
young open clusters (GEN.244) makes that take 20 to 60 seconds, nearly all of it the mean system mass of each distinct age
(`generation.system_mass`, about half a millisecond each). The API's daemon has a 60-second request limit, so a visitor's
request must never build the field. The warm-up (`planetgen warm`, the first phase, `web/warmup.py`) builds it instead, as one
stage on the one waiting line (`queue/runline.py`):

* the **parts**: the distinct ages, sorted, are cut into runs of `PART_AGES`; each unit works out the mean system mass of one
  run and writes it as its own file (`<state folder>/galaxy-field/<signature>.parts/part-NNNNN.json`) by a write and a rename,
  so a unit is a fraction of a second to a few seconds however busy the machine is, and a part that is there is kept: a run
  that stopped is carried on by the next, not started again;
* the **conclusion** (`conclude_field`, once the last part is in): the field is made from the rows with those means (a second or
  so), written as `<signature>.field` by a write and a rename, and the parts and older fields are removed.

The signature names what the field was made from: the rows (the table's signature, `db/store.component_field_key`) and the
code that makes it (`code_version`, a hash of the source files that decide the numbers), so a changed row or a changed
formula is a different file and a stale one is never read. `read` loads the file through an unpickler that allows only the
field's own classes, so a file in the state folder cannot run anything. The web side only reads (`db/store._component_field`);
a miss is `admin.cachebuild.CacheRebuilding` for a visitor.
"""

import functools
import hashlib
import importlib
import io
import json
import math
import os
import pickle
import time

from planetgen.galaxy import component_field
from planetgen.util import settings as settings_module

FORMAT = 1
"""int: The file format; a change makes every earlier file unreadable (a different signature)."""

FOLDER = "galaxy-field"
"""str: The folder of the state folder that holds the fields and their parts."""

PART_AGES = 1000
"""int: How many distinct ages one part works out (about 0.6 s at the test server's speed)."""

PART_SECONDS = 0.6
"""float: About how long one part takes, the weight of its unit in the stage's bar."""

INLINE_COMPONENTS = 1500
"""int: The most stored components for which a visitor's request still builds the field itself, since it takes about a second;
a bigger galaxy waits for the file."""

KEEP_FIELDS = 2
"""int: How many field files are kept (the newest): the one the site reads and the one before it."""

SOURCES = ("planetgen.galaxy.component_field", "planetgen.galaxy.components", "planetgen.galaxy.satellite_history",
           "planetgen.physics.component_profiles", "planetgen.physics.imf", "planetgen.physics.stellar_evolution",
           "planetgen.generation.system_mass", "planetgen.tuning")
"""tuple: The modules whose source decides the field's numbers; the signature hashes them."""

ALLOWED_CLASSES = frozenset({
    ("numpy", "dtype"), ("numpy", "ndarray"),
    ("numpy._core.multiarray", "_reconstruct"), ("numpy.core.multiarray", "_reconstruct"),
    ("numpy._core.multiarray", "scalar"), ("numpy.core.multiarray", "scalar"),
    ("numpy._core.numeric", "_frombuffer"), ("numpy.core.numeric", "_frombuffer"),
    ("planetgen.galaxy.component_field", "ComponentField"),
    ("planetgen.physics.component_profiles", "ProfileTable"), ("planetgen.physics.king", "KingTable"),
})
"""frozenset: The only classes and functions a field file may name."""


class PartsMissing(RuntimeError):
    """The parts a field is made from are not all there (a part was removed or never written)."""


class Changed(RuntimeError):
    """The stored components changed while the job ran: the parts are for rows that are gone."""


def folder():
    """The folder the fields and parts are kept in."""
    return os.path.join(settings_module.state_dir(), FOLDER)


@functools.lru_cache(maxsize=1)
def code_version():
    """A hash of the format and of the source of `SOURCES`: what the code that makes the field is."""
    digest = hashlib.sha256(f"format {FORMAT}".encode())
    for name in SOURCES:
        path = importlib.import_module(name).__file__
        with open(path, "rb") as handle:
            digest.update(name.encode() + b"\0" + handle.read())
    return digest.hexdigest()


def signature(key):
    """The name of the field for the table signature `key` (`db/store.component_field_key`), 32 hex digits."""
    return hashlib.sha256(f"{key!r}|{code_version()}".encode()).hexdigest()[:32]


def field_path(sig):
    return os.path.join(folder(), f"{sig}.field")


def parts_dir(sig):
    return os.path.join(folder(), f"{sig}.parts")


def part_path(sig, index):
    return os.path.join(parts_dir(sig), f"part-{index:05d}.json")


def _write_atomic(path, data):
    """Writes `data` (bytes) to `path` so a reader sees all of it or none: a temporary file in the same folder, then a rename."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    temporary = f"{path}.{os.getpid()}.tmp"
    try:
        with open(temporary, "wb") as handle:
            handle.write(data)
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


class _Unpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) not in ALLOWED_CLASSES:
            raise pickle.UnpicklingError(f"a field file may not name {module}.{name}")
        return super().find_class(module, name)


def _header(sig, components):
    return json.dumps({"format": FORMAT, "signature": sig, "components": components}, sort_keys=True).encode() + b"\n"


def is_published(sig):
    """Whether the field `sig` is a file and its header says it is that field (the body is not read)."""
    try:
        with open(field_path(sig), "rb") as handle:
            head = json.loads(handle.readline())
    except (OSError, ValueError):
        return False
    return head.get("format") == FORMAT and head.get("signature") == sig


def read(sig):
    """The `ComponentField` the file of `sig` holds, or `None` when there is no such file or it cannot be trusted (another
    format or signature, cut short, or naming a class outside `ALLOWED_CLASSES`)."""
    try:
        with open(field_path(sig), "rb") as handle:
            head = json.loads(handle.readline())
            if head.get("format") != FORMAT or head.get("signature") != sig:
                return None
            field = _Unpickler(io.BytesIO(handle.read())).load()
    except (OSError, ValueError, pickle.UnpicklingError, EOFError, AttributeError, ImportError, IndexError):
        return None
    return field if isinstance(field, component_field.ComponentField) else None


def distinct_ages(values):
    """The sorted distinct ages (Gy, rounded as `component_field.age_key` rounds them) of component ages `values`."""
    return sorted({round(float(value), 6) for value in values})


def part_count(ages):
    return max(1, math.ceil(len(ages) / PART_AGES))


def part_ages(ages, index):
    """The ages part `index` works out."""
    return ages[index * PART_AGES:(index + 1) * PART_AGES]


def _read_part(sig, index, ages):
    """The means of part `index` if its file is there and is for exactly these `ages`, else `None`."""
    try:
        with open(part_path(sig, index), "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return None
    if data.get("ages") != list(ages) or len(data.get("means", ())) != len(ages):
        return None
    return data["means"]


def compute_part(sig, index, ages):
    """
    Works out part `index`: the mean system mass at each of `ages`, written as the part's file. A part already there for
    these ages is left as it is.

    Returns:
        dict: `index`, `ages`, `built` (1 when it was worked out, 0 when it was there) and `seconds`.
    """
    started = time.monotonic()
    if _read_part(sig, index, ages) is not None:
        return {"index": index, "ages": len(ages), "built": 0, "seconds": 0.0}
    from planetgen.generation import system_mass
    means = [system_mass.mean_system_mass_at_age_solar(age) for age in ages]
    _write_atomic(part_path(sig, index), json.dumps({"ages": list(ages), "means": means}).encode())
    return {"index": index, "ages": len(ages), "built": 1, "seconds": time.monotonic() - started}


def assemble(sig, found, edge_pc, expected_at_density_1):
    """
    Makes the field from the parts and publishes it: the file `<sig>.field`, then the parts and older fields are removed.

    Args:
        sig (str): The signature (`signature`).
        found (iterable[galaxy.components.Component]): The stored components.

    Returns:
        component_field.ComponentField: The field that was written.

    Raises:
        PartsMissing: A part is not there or is for other ages (the next run makes it).
    """
    held = component_field.gate_components(found)
    ages = distinct_ages(c.age_gy for c in held)
    mean_by_age = {}
    for index in range(part_count(ages)):
        wanted = part_ages(ages, index)
        means = _read_part(sig, index, wanted)
        if means is None:
            raise PartsMissing(f"part {index + 1} of {part_count(ages)} of the cluster field is not there")
        mean_by_age.update(zip(wanted, means))
    field = component_field.ComponentField(found, edge_pc, expected_at_density_1, mean_by_age=mean_by_age)
    _write_atomic(field_path(sig), _header(sig, len(field)) + pickle.dumps(field, protocol=pickle.HIGHEST_PROTOCOL))
    _tidy(keep=sig)
    return field


def _tidy(keep):
    """Removes `keep`'s parts, the field files but the newest `KEEP_FIELDS`, and the parts of other galaxies a day old."""
    import shutil
    try:
        entries = list(os.scandir(folder()))
    except OSError:
        return
    fields = sorted((e for e in entries if e.name.endswith(".field")), key=lambda e: e.stat().st_mtime, reverse=True)
    for entry in fields[KEEP_FIELDS:]:
        _remove(entry.path)
    for entry in entries:
        if entry.name.endswith(".parts") and (entry.name == f"{keep}.parts" or time.time() - entry.stat().st_mtime > 24 * 3600):
            shutil.rmtree(entry.path, ignore_errors=True)


def _remove(path):
    try:
        os.unlink(path)
    except OSError:
        pass


# --- The job: the warm-up's cluster-field stage ----------------------------------------------------------------------

def _connect(database):
    """A connection to `database`: its name (this host's settings), or a `store.MySQLConfig` (a run's named host, port, user)."""
    from planetgen.db import store
    config = database if isinstance(database, store.MySQLConfig) else store.MySQLConfig(database=database)
    return store.get_connection(config, ensure_schema=False)


def _name_of(database):
    """The database's name, whether `database` is a name or a `store.MySQLConfig`."""
    return getattr(database, "database", database)


def _ages_in(conn):
    """The distinct ages of the stored components that count in the field, read as the column only."""
    marks = ", ".join("?" for _ in component_field.KINDS_IN_THE_GATE)
    rows = conn.execute(f"SELECT DISTINCT age_gy FROM galaxy_components WHERE kind IN ({marks})",
                        component_field.KINDS_IN_THE_GATE).fetchall()
    return distinct_ages(row["age_gy"] for row in rows)


def expand_field(params):
    """
    The units of the `warm-cluster-field` stage: one per part of `PART_AGES` ages, or one that says the field is already a
    file (or that the galaxy has no components).
    """
    from planetgen.queue import commands
    database = params["database"]
    conn = _connect(database)
    try:
        from planetgen.db import store
        key = store.component_field_key(conn)
        if key is None:
            return commands.Expansion([({"database": database, "state": "empty"}, 0.1)], {"state": "empty"})
        sig = signature(key)
        if is_published(sig):
            return commands.Expansion([({"database": database, "signature": sig, "state": "ready"}, 0.1)],
                                      {"state": "ready", "signature": sig, "components": key[1]})
        ages = _ages_in(conn)
    finally:
        conn.close()
    units = [({"database": database, "signature": sig, "index": index, "state": "part"}, PART_SECONDS)
             for index in range(part_count(ages))]
    return commands.Expansion(units, {"state": "parts", "signature": sig, "components": key[1], "ages": len(ages),
                                      "parts": len(units)})


def field_unit(args, report):
    """The unit of the `warm-cluster-field` command: one part, or the line that says there is nothing to make."""
    state = args["state"]
    if state == "empty":
        return {"built": 0, "cached": 0, "warmed": 0, "left_out": 0, "seconds": 0.0,
                "line": "cluster-field: the galaxy has no clusters or satellites, so there is no field."}
    if state == "ready":
        return {"built": 0, "cached": 1, "warmed": 1, "left_out": 0, "seconds": 0.0,
                "line": "cluster-field: already a file, nothing to make."}
    conn = _connect(args["database"])
    try:
        from planetgen.db import store
        key = store.component_field_key(conn)
        if key is None or signature(key) != args["signature"]:
            raise Changed("the stored clusters and satellites changed while the cluster field was made; "
                          "the next warm-up starts it again")
        ages = _ages_in(conn)
    finally:
        conn.close()
    done = compute_part(args["signature"], args["index"], part_ages(ages, args["index"]))
    return {"built": done["built"], "cached": 1 - done["built"], "warmed": 1, "left_out": 0, "seconds": done["seconds"],
            "line": f"cluster-field: part {args['index'] + 1} {'worked out' if done['built'] else 'was already there'} "
                    f"({done['ages']} ages, {done['seconds']:.1f} s)"}


def conclude_field(params, plan, results, seconds):
    """Once the last part is in: makes the field from the rows and the parts and publishes the file. A stage with nothing to
    make concludes with nothing."""
    if plan.get("state") != "parts":
        return {}
    from planetgen.admin import activity_log
    from planetgen.db import store
    conn = _connect(params["database"])
    try:
        key = store.component_field_key(conn)
        if key is None or signature(key) != plan["signature"]:
            raise Changed("the stored clusters and satellites changed while the cluster field was made; "
                          "the next warm-up starts it again")
        found = tuple(c for c in store.get_galaxy_components(conn) if c.kind in component_field.KINDS_IN_THE_GATE)
        shape = conn.execute("SELECT edge_pc, expected_system_count_at_density_1 FROM galaxy_shape WHERE id = 1").fetchone()
    finally:
        conn.close()
    import dataclasses
    started = time.monotonic()
    field = assemble(plan["signature"], tuple(dataclasses.replace(c, id=None) for c in found), float(shape["edge_pc"]),
                     float(shape["expected_system_count_at_density_1"]))
    activity_log.event("GEN", "warm.cluster_field", db=_name_of(params["database"]), components=len(field), parts=plan["parts"],
                       seconds=round(seconds + time.monotonic() - started, 1))
    return {"components": len(field), "parts": plan["parts"]}


def warm_inline(database):
    """The same work in this process, part after part (`web/warmup.warm_family`'s inline twin): the number of parts worked
    out and the number that were there."""
    plan = expand_field({"database": database})
    built = cached = 0
    for args, _cost in plan.units:
        result = field_unit(args, None)
        built, cached = built + result["built"], cached + result["cached"]
    conclude_field({"database": database}, plan.plan, [], 0.0)
    return built, cached


def ensure_file(database):
    """
    Makes sure the galaxy's cluster field is a file before the workers of a generation run need it, so that the process
    that cuts the stage builds it once and every worker reads the file, instead of each worker building it for itself (a
    minute and several hundred megabytes in each, ten workers at a time). The warm-up that normally makes the file
    refuses to run while a generation run is live, and a new or re-planned galaxy has no file, so without this the first
    stage of every run met a missing file in every worker.

    Args:
        database: The database's name, or a `store.MySQLConfig` (`wrapper.database_args`).

    Returns:
        tuple | None: `(parts worked out, parts that were there)` when the file was made here, `None` when it was already
            there or the galaxy has no field to make.
    """
    from planetgen.db import store
    conn = _connect(database)
    try:
        key = store.component_field_key(conn)
    finally:
        conn.close()
    if key is None or is_published(signature(key)):
        return None
    return warm_inline(database)
