# planetgen/web/lib/sensor.py

"""
The Sensor registry (UX.275): for every kind the property sheet knows (`describe/sheet.SPECS`), the look builder that turns the stored row into
the record a painter draws, and which painter that is. The Sensor view is the picture of an object as a ship's instruments would show it from
outside, painted from the stored row and the Object ID seed alone (`bodylook.seed_for`), never from a cached sky.

This registry only names and gathers what exists: the look records of `web/maps/bodylook.py` (stars, planets, moons), the parameters of
`web/maps/phenomenonrender.py` (neutron stars, black holes, quasars, rogue planets, interstellar comets) and the nebula view's record. The
painters are in `static/sensorview.js`, which names the same painter ids (`tests/test_sensor.py` checks the two lists agree). A kind no painter
draws yet has an entry that says so, so the guard test fails only for a kind nobody has thought about: a new kind in the property sheet has to
get an entry, with a painter or with the reason it has none.

`GET /api/objects/<ref>/sensor` returns `sensor_for` of the object, so the object page, the map's info panel and the Chart and 3D views read one
record. Nothing here changes a picture.
"""

from dataclasses import dataclass
from typing import Callable, Optional

from planetgen.db import store
from planetgen.describe.sheet import TABLES
from planetgen.web.maps import bodylook, phenomenonrender
from planetgen.web.maps.starmap import NEBULA_SURROUNDINGS_PATH, NEBULA_VOLUME_PATH, nebula_look

SENSOR_VERSION = 1
"""int: The version of the record `sensor_for` returns; a painter refuses a newer one."""

PAINTERS = ("sphere", "phenomenon", "nebula")
"""tuple: The painter ids an entry may name; `static/sensorview.js` has the same list.
- `sphere`: the lit sphere of `static/bodyview.js` (a star, planet or moon look record).
- `phenomenon`: the 3D still of `static/phenomenonrender.js` (the record is `{kind, ...parameters}`).
- `nebula`: the cloud of `static/nebulaview.js` drawn from its points."""

FACILITY_REASON = "a facility has no picture before 9.3, when stations are built"


@dataclass(frozen=True)
class Entry:
    """One kind's place in the registry: its painter and the builder of its record, or the reason it has neither."""
    kind: str
    painter: Optional[str] = None
    build: Optional[Callable] = None
    reason: Optional[str] = None
    inputs: tuple = ()

    def __post_init__(self):
        if (self.painter is None) != (self.build is None):
            raise ValueError(f"{self.kind}: a painter and a builder go together")
        if self.painter is not None and self.painter not in PAINTERS:
            raise ValueError(f"{self.kind}: unknown painter {self.painter!r}")
        if self.painter is None and not self.reason:
            raise ValueError(f"{self.kind}: an entry with no painter says why")


def _star(conn, row):
    look = bodylook.star_look(row)
    return {"look": look, "caption": bodylook.look_note(look), "estimate": _estimates(look)}


def _estimates(look):
    """What the picture guessed: a colour from temperature alone (no stored ground mixture) is an estimate."""
    return ["surface colour"] if look.get("quality") == "estimate" else []


def _with_air(conn, kind, row):
    return store.with_gases(conn, [row], kind == "moon")[0]


def _host_temperature(conn, row):
    star_id = row.get("star_id")
    if star_id is None:
        return None
    found = conn.execute("SELECT temperature_k FROM stars WHERE id = ?", (star_id,)).fetchone()
    return None if found is None else found["temperature_k"]


def _planet(kind):
    def build(conn, row):
        row = _with_air(conn, kind, row)
        row["_star_temperature_k"] = _host_temperature(conn, row)
        if kind == "planet":      # the gaps its moons dig in a ring
            row["moons"] = [dict(moon) for moon in conn.execute("SELECT * FROM moons WHERE planet_id = ?", (row["id"],)).fetchall()]
        look = bodylook.planet_look(row, kind)
        return {"look": look, "caption": bodylook.look_note(look), "estimate": _estimates(look)}
    return build


def _phenomenon(kind):
    def build(conn, row):
        params, caption = phenomenonrender.view_params(kind, row)
        return {"look": {"kind": kind, **params}, "caption": caption, "estimate": []}
    return build


def _nebula(conn, row):
    color, core, edge = nebula_look(row.get("nebula_type"))
    nebula_id = str(row["id"])
    return {"look": {"nebulaId": nebula_id, "volumePath": NEBULA_VOLUME_PATH.replace("{id}", nebula_id),
                     "surroundingsPath": NEBULA_SURROUNDINGS_PATH.replace("{id}", nebula_id),
                     "color": color, "coreOpacity": core, "edgeOpacity": edge},
            "caption": "An irregular cloud drawn from its properties.", "estimate": []}


def _none(kind, reason):
    return Entry(kind, reason=reason)


REGISTRY = {entry.kind: entry for entry in (
    Entry("star", "sphere", _star, inputs=("temperature, luminosity, class, rotation",)),
    Entry("planet", "sphere", _planet("planet"), inputs=("composition, ground, sea, air, albedo, tilt, spin, rings",)),
    Entry("moon", "sphere", _planet("moon"), inputs=("composition, ground, sea, air, tilt, spin",)),
    Entry("rogue_planet", "phenomenon", _phenomenon("rogue_planet"), inputs=("composition, internal heat, temperature",)),
    Entry("neutron_star", "phenomenon", _phenomenon("neutron_star"), inputs=("spin period, field, surface temperature",)),
    Entry("black_hole", "phenomenon", _phenomenon("black_hole"), inputs=("mass class, spin, accretion disk",)),
    Entry("quasar", "phenomenon", _phenomenon("quasar"), inputs=("radio loudness",)),
    Entry("interstellar_comet", "phenomenon", _phenomenon("interstellar_comet"), inputs=("activity",)),
    Entry("nebula", "nebula", _nebula, inputs=("type, size, the galaxy's density",)),
    _none("belt", "no painter draws a belt yet (the grain field is still to build)"),
    _none("comet", "no painter draws a comet of a system yet (a nucleus with its coma and tail is still to build)"),
    _none("asteroid_field", "no painter draws an asteroid field yet"),
    _none("supernova_remnant", "no painter draws a remnant yet; the map draws its shell from the density"),
    _none("globular_cluster", "no painter draws a cluster yet (a star field from its counts is still to build)"),
    _none("superbubble", "no painter draws a superbubble yet"),
    _none("facility", FACILITY_REASON),
)}
"""dict: Kind -> `Entry`, one for every key of `describe.sheet.SPECS`."""


def entry_for(kind):
    """The registry entry of a kind; `ValueError` for a kind with none."""
    if kind not in REGISTRY:
        raise ValueError(f"no Sensor entry for kind {kind!r}")
    return REGISTRY[kind]


def sensor_for(conn, kind, object_id):
    """
    The Sensor record of one stored object, from its row alone.

    Returns:
        dict: `{"v", "kind", "id", "painter", "look", "caption", "estimate", "inputs", "reason"}`: `painter` and `look` are `None` for a
            kind no painter draws (and `reason` says why); `estimate` lists what the picture guessed.

    Raises:
        ValueError: A kind with no entry, or no such row.
    """
    entry = entry_for(kind)
    record = {"v": SENSOR_VERSION, "kind": kind, "id": object_id, "painter": entry.painter, "look": None, "caption": None,
              "estimate": [], "inputs": list(entry.inputs), "reason": entry.reason}
    if entry.build is None:
        return record
    table = TABLES[kind]
    found = conn.execute(f"SELECT * FROM {table} WHERE id = ?", (object_id,)).fetchone()
    if found is None:
        raise ValueError(f"no {table} row with id {object_id}")
    record.update(entry.build(conn, dict(found)))
    return record
