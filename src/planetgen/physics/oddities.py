# planetgen/physics/oddities.py

"""
The oddity table: rare, physically plausible outcomes that enter the generator as a cause at one input (GEN.1036, docs/design/systems/weird-worlds.md)
===================================================================================================================================================

An oddity is never painted on a finished planet. Each row of `physics/data/oddities.json` is a **cause** that changes one input of the chain (a
lever of design section 8.1: the spin axis and rate, the orbit, the bulk composition, the interior heat, the history, the moon and ring set, or an
air state that has two branches); the solvers already built then make the surface. A row holds:

- `id`, `name`, `kind` (planet, moon, giant, rogue, small body or star system), and the `lever` and `input` it changes (the input must belong to
  the lever, so a row can only change what a lever allows);
- `requires`: the gate, conditions on stored quantities in the shape of `planet_classes.json` (`quantity` with `min`, `lt`, `equals` or
  `not_equals`), never empty: a rate is a rate on a gated population, not a global percentage;
- `rate`: the chance given the gate, `value` with its `source`, `label` (READ in the source, 2ND-HAND in a paper that quotes it, RECALLED), GEN.725
  `tier` and `tweak` level and the `registry` id of its stand-in row;
- `apply`: the named `effect` with its `params`, which the module that implements the row defines;
- `limits`: the conditions that must hold on the body after the effect (the plausibility harness checks them), never empty;
- `source` and `text`: where it is read and the one sentence the page shows ("This world turns backwards because ...").

The loader refuses a row without any of these, a quantity the file does not declare, an input outside its lever, a rate outside 0 to 1, a repeated
id, or a tier or tweak that is not in the registry's scale. An empty table loads.

**Draws.** Every roll uses a stream of its own, keyed on the oddity id and the body (`stream`), never the run's stream, so adding, removing or
turning off a row changes no other draw (design 4.4; the isolation test holds this). `probability` multiplies the row's rate by the world's
weirdness dial (`tuning.WEIRDNESS_DIAL`, 0 turns every oddity off) and the row's own scale (`tuning.WEIRDNESS_RATE_SCALE`).

Nothing calls this module yet (GEN.1047 wires the table; each effect has its own module). Use:

    rows = candidates("planet", quantities)          # the rows whose gate holds
    for row in rows:
        if roll(row, body_key=planet_seed):
            ...                                      # the module for row.apply["effect"] changes row.apply["input"]
    broken = limit_failures(row, quantities_after)    # [] when the body is still physical
"""

import json
import os
import re
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import elements, planet_classes
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "oddities.json")
"""str: The data file."""

LABELS = ("READ", "2ND-HAND", "RECALLED")
"""tuple: How far a rate was checked: read in the source, read in a paper that quotes it, or recalled and not opened (a stand-in until confirmed)."""

TIERS = ("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9")
TWEAKS = ("L0", "L1", "L2", "L3", "L4")
CONDITION_KEYS = {"quantity", "min", "lt", "equals", "not_equals"}
_ID = re.compile(r"^[A-Z]{1,2}[0-9]{1,3}$")
_NAME = re.compile(r"^[a-z][a-z0-9_]*$")


class OddityError(DataError):
    """An oddity file that is refused, or a body the table cannot be asked about."""


@dataclass(frozen=True)
class Rate:
    """The chance of an oddity given its gate, and where the number comes from."""
    value: float
    source: str
    label: str
    tier: str
    tweak: str
    registry: str


@dataclass(frozen=True)
class Oddity:
    """One row of the table."""
    id: str
    name: str
    kind: str
    lever: str
    input: str
    effect: str
    params: dict
    requires: tuple
    limits: tuple
    rate: Rate
    source: str
    text: str


def _text(row, field, owner):
    value = row.get(field)
    if not isinstance(value, str) or not value.strip():
        raise OddityError(f"{owner}: '{field}' is missing or empty")
    return value


def _conditions(row, field, owner, quantities):
    found = row.get(field)
    if not isinstance(found, list) or not found:
        raise OddityError(f"{owner}: '{field}' needs at least one condition")
    for condition in found:
        if not isinstance(condition, dict) or "quantity" not in condition or not (set(condition) - {"quantity"}):
            raise OddityError(f"{owner}: {field} condition {condition!r} is not a quantity with min, lt, equals or not_equals")
        if set(condition) - CONDITION_KEYS:
            raise OddityError(f"{owner}: {field} condition {condition!r} uses a key outside {sorted(CONDITION_KEYS)}")
        if condition["quantity"] not in quantities:
            raise OddityError(f"{owner}: {field} reads {condition['quantity']!r}, which the file does not declare")
    return tuple(found)


def _rate(row, owner):
    rate = row.get("rate")
    if not isinstance(rate, dict):
        raise OddityError(f"{owner}: 'rate' is missing")
    value = rate.get("value")
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0.0 <= value <= 1.0:
        raise OddityError(f"{owner}: rate.value {value!r} is not a number from 0 to 1")
    if rate.get("label") not in LABELS:
        raise OddityError(f"{owner}: rate.label {rate.get('label')!r} is not one of {LABELS}")
    if rate.get("tier") not in TIERS:
        raise OddityError(f"{owner}: rate.tier {rate.get('tier')!r} is not one of T1 to T9")
    if rate.get("tweak") not in TWEAKS:
        raise OddityError(f"{owner}: rate.tweak {rate.get('tweak')!r} is not one of L0 to L4")
    return Rate(float(value), _text(rate, "source", f"{owner} rate"), rate["label"], rate["tier"], rate["tweak"],
                _text(rate, "registry", f"{owner} rate"))


def parse(document):
    """
    The table of a loaded document: `(rows, levers, quantities)`.

    Raises:
        OddityError: For any row or vocabulary entry the schema refuses (see the module docstring).
    """
    kinds = tuple(document.get("kinds", ()))
    levers = document.get("levers")
    quantities = document.get("quantities")
    if not kinds or not isinstance(levers, dict) or not levers or not isinstance(quantities, dict):
        raise OddityError("the file needs 'kinds', 'levers' and 'quantities'")
    for key, lever in levers.items():
        if not lever.get("inputs") or not lever.get("label"):
            raise OddityError(f"lever {key}: needs a label and inputs")
    rows = []
    seen = set()
    for row in document.get("rows", ()):
        owner = f"oddity {row.get('id')!r}"
        identifier = row.get("id")
        if not isinstance(identifier, str) or not _ID.match(identifier):
            raise OddityError(f"{owner}: the id is not a letter or two and a number (W10, G1, R4)")
        if identifier in seen:
            raise OddityError(f"{owner}: the id is used by two rows")
        seen.add(identifier)
        if row.get("kind") not in kinds:
            raise OddityError(f"{owner}: kind {row.get('kind')!r} is not one of {kinds}")
        apply = row.get("apply")
        if not isinstance(apply, dict):
            raise OddityError(f"{owner}: 'apply' is missing")
        lever = apply.get("lever")
        if lever not in levers:
            raise OddityError(f"{owner}: apply.lever {lever!r} is not one of {sorted(levers)}")
        if apply.get("input") not in levers[lever]["inputs"]:
            raise OddityError(f"{owner}: apply.input {apply.get('input')!r} is not an input of lever {lever} ({levers[lever]['inputs']})")
        effect = apply.get("effect")
        if not isinstance(effect, str) or not _NAME.match(effect):
            raise OddityError(f"{owner}: apply.effect {effect!r} is not a lower-case name")
        params = apply.get("params", {})
        if not isinstance(params, dict):
            raise OddityError(f"{owner}: apply.params is not an object")
        rows.append(Oddity(
            identifier, _text(row, "name", owner), row["kind"], lever, apply["input"], effect, dict(params),
            _conditions(row, "requires", owner, quantities), _conditions(row, "limits", owner, quantities),
            _rate(row, owner), _text(row, "source", owner), _text(row, "text", owner)))
    return tuple(rows), levers, quantities


def load(path=DATA_FILE):
    """The table in the file at `path`: `(rows, levers, quantities)` (`parse` of its JSON)."""
    with open(path, encoding="utf-8") as handle:
        return parse(json.load(handle))


ROWS, LEVERS, QUANTITIES = load()
"""tuple, dict, dict: The shipped table's `Oddity` rows, the levers (each with its inputs) and the quantities a gate or a limit may read."""

BY_ID = {row.id: row for row in ROWS}


def candidates(kind, found, rows=None):
    """The rows of `kind` whose gate holds for the named quantities in `found` (a quantity the body lacks fails the gate)."""
    return [row for row in (ROWS if rows is None else rows)
            if row.kind == kind and planet_classes.holds(row.requires, found, owner=f"oddity {row.id}", missing=None)]


def limit_failures(row, found):
    """The limits of `row` that do not hold on the body after the effect, as the conditions themselves ([] when it is still physical). A limit
    on a quantity the body lacks counts as failed: an oddity must leave every quantity its limits read."""
    return [limit for limit in row.limits if not planet_classes.holds([limit], found, owner=f"oddity {row.id}", missing=None)]


def probability(row, dial=None, scale=None):
    """The chance that a body that passes the gate gets `row`: its rate times the weirdness dial and the row's own scale, held to 0..1.
    `dial` and `scale` (a dict by id) default to `tuning.WEIRDNESS_DIAL` and `tuning.WEIRDNESS_RATE_SCALE`."""
    dial = tuning.WEIRDNESS_DIAL if dial is None else dial
    scale = tuning.WEIRDNESS_RATE_SCALE if scale is None else scale
    return min(1.0, max(0.0, row.rate.value * dial * scale.get(row.id, 1.0)))


def stream(row_id, body_key):
    """The stream of one oddity on one body: seeded by the id and the body's key alone (a string seed, hashed by SHA-512, so it does not depend on
    `PYTHONHASHSEED`), so no other row and no other draw of the run moves when a row is added, removed or turned off."""
    return draw.Stream(f"oddity:{row_id}:{body_key}")


def roll(row, body_key, dial=None, scale=None):
    """Whether the body keyed `body_key` gets `row`. The draw comes from the row's own stream and is made even at dial 0, so the stream a body
    has for a row never depends on the dial."""
    draw_value = stream(row.id, body_key).random()
    return draw_value < probability(row, dial, scale)


def words(row):
    """The sentence the page shows for `row`: 'Odd: <name>. <text>'."""
    return f"Odd: {row.name}. {row.text}"
