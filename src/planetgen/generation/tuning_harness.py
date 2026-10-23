# planetgen/generation/tuning_harness.py

"""
The shared tuning harness (GEN.503, docs/design/tunable-gaps.md sections 1 and 4)
===================================================================================

A composition value with no sourced weight can be a seeded draw whose knobs are tuned to a published population statistic and
stay tunable. This module is the one build every tuned draw (builders B1 to B16, GEN.509 to GEN.525) uses:

* **The target table** (`physics/data/tuning_targets.json`, `TARGETS`). One `Target` per published statistic: the row it tunes
  (`"<table>:<row>"`, as the possibility-table registry names it), the statistic (`mean`, `sd`, `median` or `proportion`), the
  published value and its 1-sigma uncertainty, the source and the label of how far it was checked (READ, DATA, COMPUTED,
  2ND-HAND, SNIPPET, RECALLED, NOT REACHED), the validity range of the data, the id of the selection model the simulated side
  has to apply, and the fill id and falsifier of the stand-in it replaces (gap-filling-plan.md section 6.3).
* **The comparison** (`compare`, `simulate_and_compare`). The simulated statistic is compared with the target within 2 sigma
  (chi-square of one degree of freedom, at most 4) where sigma combines the published uncertainty and the sampling error of
  the simulation; a proportion is compared by the one-sigma Wilson intervals of both sides, on the side that faces the other.
  Draws are made one per seed, on a stream bound to `tuning:<target id>:<i>`, so a run repeats.
* **The loader check** (`TargetTable.require`, run by `Recipes` at load and by `Registry.validate(targets=...)`). A row marked
  `tunable` with no target is refused with `TargetError`, so a knob cannot be tuned to nothing (extends TEST.214).
* **Scope.** A target's `validity` is the range of each stored quantity the data cover (for example `host_type` WD, [Fe/H]
  -0.5 to +0.5, FGK for the radius valley). `Target.outside(quantities)` names the quantities a subject has outside it, and a
  draw there is a flagged extrapolation (`Comparison.extrapolated`).
* **Solar System values.** A value measured only in the Solar System is not a target (NO FIXED CONTENT): it is a `test-point`
  with the `solar_system` flag stating why, reported and never tuned to (`Comparison.tunes` is false). The loader refuses a
  `target` whose statistic or source names a Solar System body, or any entry with the flag but not the role `test-point`.
* **Selection models.** A statistic of an observed sample needs the same cut on the simulated side. `SELECTION_MODELS` holds the
  five of tunable-gaps.md section 4 by id (built by GEN.504 to GEN.508); a model nobody has registered yet makes
  `simulate_and_compare` refuse (`TargetError`), so a tuned width is never reported as a match without its cut.

A retune changes a table's draws, so it is a table version change and a generator epoch bump (possibility_table.py); adding
a target changes no draw. `TargetTable.digest` hashes the table so a test can pin it.

Use:

    target = TARGETS.get("a4-mean")
    result = simulate_and_compare(target, recipe_sampler("a4", subject), n=2000)
    result.passed, result.z, result.simulated, result.extrapolated
"""

import hashlib
import json
import math
import os
import statistics
from dataclasses import dataclass

from planetgen.physics import elements
from planetgen.physics.elements import DataError
from planetgen.util import draw

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "tuning_targets.json")
"""str: The data file."""

LABELS = ("READ", "DATA", "COMPUTED", "2ND-HAND", "SNIPPET", "RECALLED", "NOT REACHED")
"""tuple: How far a target's number was checked (tunable-gaps.md, "Worker numbers keep their own labels")."""

USABLE_LABELS = ("READ", "DATA", "COMPUTED", "2ND-HAND")
"""tuple: The labels a target may be tuned to (Boss, 3:02 AM EDT 10-17: verified and partly verified are usable). A snippet, a
recalled value and a paper not reached are listed in the table and never tuned to."""

STATISTICS = ("mean", "sd", "median", "proportion")
"""tuple: What a target states about the simulated sample."""

ROLES = ("target", "test-point")

SIGMAS = 2.0
"""float: The comparison passes within this many sigma."""

SELECTION_MODEL_ITEMS = {
    "kepler-recovery": "GEN.504",
    "mass-radius-sample": "GEN.505",
    "jwst-target-selection": "GEN.506",
    "comet-admission": "GEN.507",
    "wd-pollution-detection": "GEN.508",
}
"""dict: The five selection models and the item that builds each (tunable-gaps.md section 4)."""

SOLAR_SYSTEM_WORDS = ("sun", "solar system", "earth", "venus", "mars", "mercury", "moon", "luna", "jupiter", "saturn", "uranus",
                      "neptune", "titan", "pluto", "triton", "io", "europa", "ganymede", "callisto", "enceladus", "ceres",
                      "vesta", "halley", "67p", "hale-bopp", "chondrite")
"""tuple: Names that mark a statistic measured on a Solar System body; read as whole words, lower case."""


class TargetError(ValueError):
    """A target table, a target or a tunable row that the harness refuses."""


@dataclass(frozen=True)
class Target:
    """
    One published statistic a tunable row is tuned to.

    Attributes:
        id (str): The target's name, unique in the table.
        row (str): The tunable row, `"<table>:<row>"` (`"atmosphere-a:a4"`).
        statistic (str): What the statistic is, in words, with the population it is measured on.
        kind (str): One of `STATISTICS`.
        field (str): For a sampler that returns a mapping (a planet with its period, radius and host, which a selection model
            needs whole): the key the statistic reads. Empty when the sampler returns the number itself.
        value (float): The published value (a proportion for `proportion`).
        sigma (float): The 1-sigma uncertainty of the published value; for a `proportion` with `counts` it is derived from them.
        counts (tuple or None): For a `proportion`: `(k, n)`, the published count.
        n (int or None): The number of objects behind a `mean`, `sd` or `median`.
        source (str): Where the number was read.
        label (str): One of `LABELS`.
        validity (dict): `{quantity: (minimum, maximum)}` or `{quantity: value}`; a missing bound is open.
        selection_model (str or None): The id in `SELECTION_MODEL_ITEMS` the simulated side applies.
        role (str): `target` (tuned to) or `test-point` (reported, never tuned to).
        solar_system (str): For a `test-point`: why the value is Solar System only.
        fill_id, falsifier (str): The gap-filling plan's fill id and the check that falsifies the fill (plan section 6.3).
    """

    id: str
    row: str
    statistic: str
    kind: str
    value: float
    sigma: float
    source: str
    label: str
    counts: tuple | None = None
    n: int | None = None
    validity: tuple = ()
    selection_model: str | None = None
    field: str = ""
    role: str = "target"
    solar_system: str = ""
    fill_id: str = ""
    falsifier: str = ""

    @property
    def table(self):
        return self.row.split(":", 1)[0]

    @property
    def row_id(self):
        return self.row.split(":", 1)[1]

    @property
    def tunes(self):
        """True when a knob may be tuned to this target: its role is `target` and its label is a usable one."""
        return self.role == "target" and self.label in USABLE_LABELS

    def outside(self, quantities):
        """The quantities of `quantities` (a dict or a `Subject`) that lie outside the validity range, as `[name]`.
        A quantity the subject lacks counts as outside."""
        get = quantities.get if hasattr(quantities, "get") else quantities.__getitem__
        names = []
        for name, bound in self.validity:
            value = get(name)
            if value is None:
                names.append(name)
            elif isinstance(bound, tuple):
                low, high = bound
                if (low is not None and value < low) or (high is not None and value > high):
                    names.append(name)
            elif value != bound:
                names.append(name)
        return names


def _label_words(text):
    out, word = [], []
    for ch in str(text).lower():
        if ch.isalnum() or ch == "-":
            word.append(ch)
        else:
            out.append("".join(word))
            word = []
    out.append("".join(word))
    return {w for w in out if w}


def _solar_system_hits(text):
    lowered = str(text).lower()
    words = _label_words(lowered)
    hits = [w for w in SOLAR_SYSTEM_WORDS if (" " not in w and w in words) or (" " in w and w in lowered)]
    return sorted(hits)


def wilson_interval(k, n, z=1.0):
    """The Wilson score interval `(low, high)` of a proportion from `k` of `n` at `z` sigma (1 is the 68 percent interval)."""
    if n <= 0 or k < 0 or k > n:
        raise ValueError("a count needs 0 <= k <= n and n > 0")
    p = k / n
    z2 = z * z
    centre = (p + z2 / (2 * n)) / (1 + z2 / n)
    half = z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / (1 + z2 / n)
    return max(0.0, centre - half), min(1.0, centre + half)


class TargetTable:
    """
    The tuning-target table: `Target`s by id, loaded and checked from a document.

    Args:
        document (dict): The data file's content.
    """

    def __init__(self, document):
        self.document = document
        self.version = document.get("version")
        if not isinstance(self.version, int) or self.version < 1:
            raise TargetError("the target table has no version")
        raw = document.get("targets")
        if not raw:
            raise TargetError("the target table has no targets")
        self._targets = {}
        for item in raw:
            target = self._target(item)
            if target.id in self._targets:
                raise TargetError(f"target {target.id} appears twice")
            self._targets[target.id] = target

    @staticmethod
    def _target(raw):
        ident = raw.get("id")
        if not ident or not isinstance(ident, str):
            raise TargetError("a target has no id")
        known = {"id", "row", "statistic", "kind", "value", "sigma", "counts", "n", "source", "label", "validity",
                 "selection_model", "field", "role", "solar_system", "fill_id", "falsifier"}
        extra = sorted(set(raw) - known)
        if extra:
            raise TargetError(f"target {ident}: unknown fields {', '.join(extra)}")
        for key in ("row", "statistic", "kind", "source", "label"):
            if not str(raw.get(key, "")).strip():
                raise TargetError(f"target {ident} has no {key}")
        if ":" not in raw["row"]:
            raise TargetError(f"target {ident}: the row {raw['row']!r} is not '<table>:<row>'")
        if raw["kind"] not in STATISTICS:
            raise TargetError(f"target {ident}: kind {raw['kind']!r} is none of {', '.join(STATISTICS)}")
        if raw["label"] not in LABELS:
            raise TargetError(f"target {ident}: label {raw['label']!r} is none of {', '.join(LABELS)}")
        role = raw.get("role", "target")
        if role not in ROLES:
            raise TargetError(f"target {ident}: role {role!r} is none of {', '.join(ROLES)}")
        counts = raw.get("counts")
        n = raw.get("n")
        if raw["kind"] == "proportion":
            if counts is None:
                raise TargetError(f"target {ident}: a proportion states the published count [k, n]")
            k, total = counts
            if not (isinstance(k, int) and isinstance(total, int) and 0 <= k <= total and total > 0):
                raise TargetError(f"target {ident}: the count {counts!r} is not 0 <= k <= n")
            counts = (k, total)
            value = k / total
            if "value" in raw and abs(raw["value"] - value) > 1e-6:
                raise TargetError(f"target {ident}: the value {raw['value']} is not the count {k}/{total}")
            low, high = wilson_interval(k, total)
            sigma = max(value - low, high - value)
        else:
            if counts is not None:
                raise TargetError(f"target {ident}: only a proportion has a count")
            for key in ("value", "sigma"):
                if not isinstance(raw.get(key), (int, float)) or isinstance(raw.get(key), bool) or not math.isfinite(raw[key]):
                    raise TargetError(f"target {ident}: {key} is not a finite number")
            value, sigma = float(raw["value"]), float(raw["sigma"])
            if not sigma > 0:
                raise TargetError(f"target {ident}: sigma must be positive")
        if n is not None and (not isinstance(n, int) or n < 2):
            raise TargetError(f"target {ident}: n is a whole number of 2 or more")
        validity = []
        for name, bound in (raw.get("validity") or {}).items():
            if isinstance(bound, list):
                if len(bound) != 2 or (bound[0] is not None and bound[1] is not None and bound[0] > bound[1]):
                    raise TargetError(f"target {ident}: the validity range of {name} is not [minimum, maximum]")
                bound = (bound[0], bound[1])
            validity.append((name, bound))
        if not validity:
            raise TargetError(f"target {ident} states no validity range (a knob's scope is the range of the data)")
        model = raw.get("selection_model")
        if model is not None and model not in SELECTION_MODEL_ITEMS:
            raise TargetError(f"target {ident}: selection model {model!r} is none of {', '.join(SELECTION_MODEL_ITEMS)}")
        flag = str(raw.get("solar_system", "")).strip()
        if flag and role != "test-point":
            raise TargetError(f"target {ident}: a Solar System value is a test-point and never a target to tune to")
        if role == "test-point" and not flag:
            raise TargetError(f"target {ident}: a test-point states why it is Solar System only (solar_system)")
        if role == "target":
            hits = _solar_system_hits(raw["statistic"]) + _solar_system_hits(raw["source"])
            if hits:
                raise TargetError(f"target {ident}: names a Solar System body ({', '.join(sorted(set(hits)))}) and has no "
                                  f"solar_system flag; flag it as a test-point or reword the statistic")
        return Target(id=ident, row=raw["row"], statistic=raw["statistic"], kind=raw["kind"], value=value, sigma=sigma,
                      source=raw["source"], label=raw["label"], counts=counts, n=n, validity=tuple(validity),
                      selection_model=model, field=str(raw.get("field", "")), role=role, solar_system=flag, fill_id=raw.get("fill_id", ""),
                      falsifier=raw.get("falsifier", ""))

    def __len__(self):
        return len(self._targets)

    def __iter__(self):
        return iter(self._targets.values())

    def __contains__(self, ident):
        return ident in self._targets

    def get(self, ident):
        try:
            return self._targets[ident]
        except KeyError:
            raise KeyError(f"no tuning target {ident!r}") from None

    def for_row(self, row):
        """The targets of `row` (`"<table>:<row>"`), in file order."""
        return [t for t in self._targets.values() if t.row == row]

    def tunable_targets(self, row):
        """The targets of `row` a knob may be tuned to (`Target.tunes`)."""
        return [t for t in self.for_row(row) if t.tunes]

    def require(self, tables):
        """
        The loader check: every row of `tables` marked `tunable` has at least one target that may be tuned to, and every
        target names a tunable row of a table in `tables` it can see. Raises `TargetError` naming the rows.

        A target for a table not in `tables` is left alone (the tables load one module at a time).
        """
        by_name = {table.name: table for table in tables}
        missing = []
        for table in by_name.values():
            for row in table.rows:
                if row.tunable and not self.tunable_targets(f"{table.name}:{row.name}"):
                    missing.append(f"{table.name}:{row.name}")
        if missing:
            raise TargetError("tunable rows with no tuning target in the target table: " + ", ".join(missing))
        stray = []
        for target in self._targets.values():
            table = by_name.get(target.table)
            if table is None:
                continue
            try:
                row = table.row(target.row_id)
            except KeyError:
                stray.append(f"{target.id} ({target.row}: no such row)")
                continue
            if target.role == "target" and not row.tunable:
                stray.append(f"{target.id} ({target.row}: the row is not marked tunable)")
        if stray:
            raise TargetError("targets that name no tunable row: " + ", ".join(stray))

    def digest(self):
        """A hash of the table's content, so a test can pin it (a retune is a table version change)."""
        return hashlib.sha256(json.dumps(self.document, sort_keys=True).encode("utf-8")).hexdigest()


# --- selection models -------------------------------------------------------------------------------------------------

SELECTION_MODELS = {}
"""dict: id -> function `item -> probability` of the selection models built so far (GEN.504 to GEN.508)."""


def register_selection_model(model_id, function):
    """Registers the selection model `model_id` (one of `SELECTION_MODEL_ITEMS`); `function(item)` returns the probability
    in [0, 1] that a simulated item (the number or the mapping the sampler returns) is in the observed sample."""
    if model_id not in SELECTION_MODEL_ITEMS:
        raise TargetError(f"selection model {model_id!r} is none of {', '.join(SELECTION_MODEL_ITEMS)}")
    SELECTION_MODELS[model_id] = function
    return function


# --- comparison -------------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Comparison:
    """
    The result of comparing a simulated statistic with a target.

    Attributes:
        target (Target): The target.
        simulated (float): The simulated statistic.
        sigma (float): The combined 1-sigma uncertainty the difference is measured in.
        z (float): `(simulated - value) / sigma`.
        passed (bool): `abs(z) <= SIGMAS`.
        draws (int): How many simulated values the statistic used.
        extrapolated (tuple): The validity quantities the subject was outside (empty when none).
        tunes (bool): `Target.tunes`: false means the result is a report and no knob may be tuned to it.
    """

    target: Target
    simulated: float
    sigma: float
    z: float
    passed: bool
    draws: int
    extrapolated: tuple = ()
    tunes: bool = True

    @property
    def chi_square(self):
        return self.z * self.z


def _statistic(kind, values):
    n = len(values)
    if kind == "mean":
        mean = math.fsum(values) / n
        sd = statistics.stdev(values) if n > 1 else 0.0
        return mean, sd / math.sqrt(n)
    if kind == "sd":
        sd = statistics.stdev(values)
        return sd, sd / math.sqrt(2.0 * (n - 1))
    if kind == "median":
        sd = statistics.stdev(values)
        return statistics.median(values), 1.2533141373155003 * sd / math.sqrt(n)
    raise ValueError(f"no sample statistic for {kind!r}")


def values_of(target, items):
    """The numbers the target's statistic reads from simulated `items`: the items themselves, or `item[target.field]` of
    mappings (a proportion's items are truthy or falsy)."""
    out = []
    for item in items:
        if hasattr(item, "keys"):
            if not target.field:
                raise TargetError(f"target {target.id}: the sampler returns mappings and the target names no field")
            item = item[target.field]
        out.append(item)
    return out


def compare(target, values, *, quantities=None, count=None):
    """
    Compares the statistic of the simulated `values` with `target`.

    Args:
        target (Target): What to compare with.
        values (sequence): The simulated values (floats) for a `mean`, `sd` or `median` target; for a `proportion`, the
            booleans of whether each draw had the property (or give `count=(k, n)` instead). Mappings are read at
            `target.field` (`values_of`).
        quantities (dict or Subject, optional): The stored quantities of the subject drawn, to flag an extrapolation.
        count (tuple, optional): `(k, n)` of a proportion already counted.

    Returns:
        Comparison: The simulated statistic and whether it is within `SIGMAS` sigma.
    """
    values = values_of(target, values)
    if target.kind == "proportion":
        if count is None:
            flags = [bool(v) for v in values]
            count = (sum(flags), len(flags))
        k, n = count
        if n < 1:
            raise TargetError(f"target {target.id}: no simulated draws")
        simulated = k / n
        lo_s, hi_s = wilson_interval(k, n)
        k_obs, n_obs = target.counts
        lo_t, hi_t = wilson_interval(k_obs, n_obs)
        if simulated >= target.value:
            sigma = math.hypot(hi_t - target.value, simulated - lo_s)
        else:
            sigma = math.hypot(target.value - lo_t, hi_s - simulated)
        z = (simulated - target.value) / sigma if sigma > 0 else math.inf
        draws = n
    else:
        values = [float(v) for v in values]
        if len(values) < 2:
            raise TargetError(f"target {target.id}: needs at least 2 simulated draws")
        simulated, sampling = _statistic(target.kind, values)
        sigma = math.hypot(target.sigma, sampling)
        z = (simulated - target.value) / sigma
        draws = len(values)
    outside = tuple(target.outside(quantities)) if quantities is not None else ()
    return Comparison(target=target, simulated=simulated, sigma=sigma, z=z, passed=abs(z) <= SIGMAS, draws=draws,
                      extrapolated=outside, tunes=target.tunes)


def simulate(target, sampler, n=2000):
    """
    Draws `n` items with `sampler(i)` (a number, or a mapping of the quantities a selection model reads), one per seed, each on a stream bound to `tuning:<target id>:<i>`, and applies the
    target's selection model: a draw enters the sample with the model's probability (decided on its own stream), and an
    item `None` (the row's gate did not pass) is skipped. Returns the kept items.

    Raises:
        TargetError: The target names a selection model that is not built yet.
    """
    model = None
    if target.selection_model is not None:
        model = SELECTION_MODELS.get(target.selection_model)
        if model is None:
            raise TargetError(f"target {target.id} needs the selection model {target.selection_model!r}, which "
                              f"{SELECTION_MODEL_ITEMS[target.selection_model]} builds; it is not registered")
    values = []
    for i in range(n):
        stream = draw.Stream(f"tuning:{target.id}:{i}")
        with draw.bound(stream):
            value = sampler(i)
            if value is None:
                continue
            if model is not None and not stream.random() < model(value):
                continue
        values.append(value)
    return values


def simulate_and_compare(target, sampler, n=2000, quantities=None):
    """`simulate` then `compare`: the one call a test or a retune makes for a target."""
    return compare(target, simulate(target, sampler, n), quantities=quantities)


def recipe_sampler(row_id, subject, recipes=None):
    """A sampler for `simulate`: the continuous quantity of the `param` row `row_id` of the atmosphere recipes at `subject`."""
    if recipes is None:
        from planetgen.physics.atmosphere_recipes import RECIPES as recipes
    return lambda _i: recipes.draw_value(row_id, subject)


def load(path=DATA_FILE):
    """Reads and checks the target table at `path`."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the tuning-target table {path}: {error}") from error
    try:
        return TargetTable(document)
    except TargetError as error:
        raise DataError(f"{path}: {error}") from error


TARGETS = load()
"""TargetTable: The tuning-target table of the generator."""
