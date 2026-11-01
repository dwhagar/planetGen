# tests/audit/sweep_tables.py

"""
The possibility-table sweep (TEST.214): walks every table of `possibility_table.REGISTRY`, so a table added later is covered without a
change to the test, and checks the sums and gates the composition model promises (docs/design/composition/composition-model.md, 4a and 4b).

For each table and each class the table names, `subject_samples` builds `SAMPLES` subjects: the quantities the gate and the weights read
(found by running the table on a subject that records what is asked of it) get values from a pool that holds the boundaries of every
condition on that quantity, a few values inside each side of them, the names of the rows of the table a `drawn.<letter>` quantity refers to,
and a physical range for the numbers (`RANGES`). The rules:

- `row-states`: every row has a source and a state, an open row has the weight `open`, no other row has, a partly verified row names its
  limit and a stand-in names its replacement item (no licence assertion: the licences are in `DATA-LICENSE.md`).
- `gate-not-empty`: every class a table names is supported by some row for some of the samples.
- `weights-sum-to-one`: where a subject has probabilities, they are above 0, at most 1 and add up to 1 within `SUM_TOLERANCE`.
- `no-failed-subjects`: no sample of a named class fails for a reason but an empty gate (a sourced sum over 1 beside open rows, every
  weight 0, a weight that raises).
- `fractions-sum-to-one`: the fractions each row draws add up to 1 within `SUM_TOLERANCE`, for a sample of draws, and name only
  components of the shared table.
- `draw-follows-weights`: the frequency with which `Table.draw` picks each row is its probability, within five binomial sigma.
- `rows-reachable`: every row is given a probability for every class it names by some subject (the row's own conditions satisfied, then
  400 tries), unless `tests/data/possibility_sweep_notes.json` says why not.
- `open-rows-listed`: the open rows are those of `tests/data/possibility_sweep_notes.json`; a new one needs its entry with a note.
- `all-registered`: every `Table` a freshly started process holds after importing the modules that use the loader is in the registry.
"""

import collections
import gc
import importlib
import json
import math
import pathlib
import pkgutil
import random
import re
import subprocess
import sys
import zlib
from dataclasses import dataclass

SAMPLES = 60
"""int: Subjects built per table and class."""

DRAWS = 1500
"""int: Draws of `Table.draw` per subject for the frequency check."""

FRACTION_DRAWS = 200
"""int: Draws of the fractions of each row that has them."""

NOTES_FILE = pathlib.Path(__file__).resolve().parents[1] / "data" / "possibility_sweep_notes.json"
"""Path: The open rows with a note each, and the tables whose rows are independent flags (an all-zero answer is allowed there)."""

RANGES = {
    "age_gy": (0.01, 13.0), "mass_earth": (0.01, 4000.0), "radius_km": (100.0, 80000.0), "teq_k": (20.0, 3500.0), "ts_k": (20.0, 2500.0),
    "tp_k": (20.0, 3500.0), "tday_k": (50.0, 4000.0), "teff_k": (300.0, 40000.0), "t_layer_k": (20.0, 3000.0),
    "profile_t_min_k": (20.0, 500.0), "profile_t_max_k": (500.0, 4000.0), "insolation": (1e-4, 1e5), "luminosity_solar": (1e-4, 1e5),
    "host_mass_solar": (0.08, 60.0), "pressure_bar": (1e-9, 1e5), "photosphere_bar": (1e-3, 100.0), "v_esc_kms": (0.05, 300.0),
    "xuv_index": (0.01, 100.0), "snow_ratio": (0.01, 100.0), "fe_h": (-2.0, 0.8), "redox_diw": (-4.0, 4.0), "c_over_o": (0.05, 3.0),
    "diameter_m": (1.0, 3e6), "x": (0.2, 6.0), "prevalence": (0.01, 1.0), "perihelion_au": (0.01, 1000.0), "o2_pal": (0.0, 2.0),
    "roche_ratio": (0.2, 20.0), "runaway_ratio": (0.01, 30.0), "sulfur_flux_earth": (0.0, 100.0), "u_star": (0.0, 5.0),
    "u_star_threshold": (0.1, 2.0), "temperature_k": (10.0, 2500.0),
}
"""dict: `{quantity: (low, high)}` for the numbers a table reads, drawn log-uniformly (uniformly where the low end is 0 or below)."""

CLASS_RANGES = {("brown-dwarf", "teff_k"): (250.0, 2400.0), ("comet", "perihelion_au"): (0.05, 60.0)}
"""dict: `{(class, quantity): (low, high)}` where a class has a narrower physical range than `RANGES`."""

STRINGS = {
    "host_type": ("M", "K", "G", "F", "A", "WD", "PSR"),
    "zone": ("h", "e", "r", "c"),
    "period_class": ("halley_type", "jupiter_family", "long_period"),
    "family": ("carbonaceous", "stony", "metallic", "icy", "basaltic", "mixed", "dust", "collisional"),
}
"""dict: `{quantity: values}` for the quantities whose values are names."""


def seed_for(*parts):
    return zlib.crc32(":".join(str(part) for part in parts).encode()) or 1


# --- Loading the registry ------------------------------------------------------------------------------------------------------


def load_tables():
    """Imports every module that uses the loader, so each registers its tables; returns the registry."""
    import planetgen
    from planetgen.generation import possibility_table

    root = pathlib.Path(planetgen.__file__).parent
    skip = {"web"}
    for info in pkgutil.walk_packages(planetgen.__path__, "planetgen."):
        name = info.name
        if any(part in skip for part in name.split(".")):
            continue
        relative = root.joinpath(*name.split(".")[1:])
        source = relative.with_suffix(".py") if relative.with_suffix(".py").exists() else relative / "__init__.py"
        try:
            text = source.read_text(encoding="utf-8")
        except OSError:
            continue
        if "possibility_table" in text:
            importlib.import_module(name)
    return possibility_table.REGISTRY


# --- Subjects --------------------------------------------------------------------------------------------------------------------


def _recording_subject_class():
    from planetgen.generation import possibility_table as pt

    class Recording(pt.Subject):
        __slots__ = ("asked",)

        def __init__(self, body_class, quantities=None):
            super().__init__(body_class, quantities)
            self.asked = []

        def __getitem__(self, name):
            if name not in self.asked:
                self.asked.append(name)
            return super().__getitem__(name)

        def get(self, name, default=None):
            if name not in self.asked:
                self.asked.append(name)
            return super().get(name, default)

    return Recording


def _condition_values(table, quantity):
    values = set()
    for row in table.rows:
        for condition in row.conditions:
            if condition.quantity != quantity:
                continue
            for value in (condition.minimum, condition.maximum, condition.equals, *(condition.one_of or ())):
                if value is not None:
                    values.add(value)
    return values


def _pool(table, quantity, registry):
    """The values a sample may give `quantity`: strings and row names as they are, numbers as the boundaries of the conditions and
    their neighbours, with the physical range of `RANGES` filled in by the caller."""
    if quantity.startswith("drawn."):
        other = f"atmosphere-{quantity.split('.', 1)[1]}"
        names = [row.name for row in registry.get(other).rows] if other in registry else []
        return list(names) or ["unset"]
    values = _condition_values(table, quantity)
    if quantity in STRINGS:
        return sorted(set(STRINGS[quantity]) | {value for value in values if isinstance(value, str)})
    if any(isinstance(value, str) for value in values):
        return sorted(value for value in values if isinstance(value, str))
    pool = set()
    for value in values:
        if isinstance(value, bool):
            pool |= {True, False}
        elif isinstance(value, (int, float)):
            pool |= {value, value * 0.999, value * 1.001, value * 0.5, value * 2.0} if value else {0, 0.001}
    return sorted(pool, key=lambda item: (str(type(item)), item))


def _pick(table, quantity, rng, registry, body_class=None):
    pool = _pool(table, quantity, registry)
    window = CLASS_RANGES.get((body_class, quantity)) or RANGES.get(quantity)
    if window and (not pool or rng.random() < 0.6):
        low, high = window
        return rng.uniform(low, high) if low <= 0.0 else math.exp(rng.uniform(math.log(low), math.log(high)))
    if pool:
        return rng.choice(pool)
    return rng.choice((0.0, 0.01, 0.3, 1.0, 3.0, 10.0, 100.0, 1000.0))


def _satisfying(row, quantity, rng, body_class):
    """A value of `quantity` that passes every condition of `row` on it, or `None` if the row has none on it."""
    conditions = [condition for condition in row.conditions if condition.quantity == quantity]
    if not conditions:
        return None
    low, high, fixed, members = None, None, None, None
    for condition in conditions:
        low = condition.minimum if low is None else (low if condition.minimum is None else max(low, condition.minimum))
        high = condition.maximum if high is None else (high if condition.maximum is None else min(high, condition.maximum))
        fixed = condition.equals if condition.equals is not None else fixed
        members = condition.one_of if condition.one_of is not None else members
    if fixed is not None:
        return fixed
    if members:
        return rng.choice(sorted(members, key=str))
    window = CLASS_RANGES.get((body_class, quantity)) or RANGES.get(quantity) or (0.0, 1000.0)
    low = window[0] if low is None else low
    high = max(window[1], low * 2.0 if low > 0 else 1.0) if high is None else high
    return rng.uniform(low, high)


def _couple(quantities, forced, rng):
    """Makes the numbers a sample did not fix agree as a body's do: the escape speed follows the mass (a rocky mass-radius law), the
    equilibrium temperature the insolation, the surface and layer temperatures the equilibrium one."""
    if "mass_earth" in quantities and "v_esc_kms" in quantities and "v_esc_kms" not in forced:
        mass = max(float(quantities["mass_earth"]), 1e-4)
        radius = max(0.3, min(mass ** 0.28, 12.0))
        quantities["v_esc_kms"] = 11.186 * math.sqrt(mass / radius)
    if "insolation" in quantities and "teq_k" in quantities and "teq_k" not in forced:
        quantities["teq_k"] = 278.6 * max(float(quantities["insolation"]), 1e-6) ** 0.25
    if "teq_k" in quantities:
        for name in ("ts_k", "tp_k", "t_layer_k", "tday_k"):
            if name in quantities and name not in forced:
                quantities[name] = float(quantities["teq_k"]) * rng.uniform(0.8, 1.6)
    return quantities


def _asked(table, body_class, registry):
    """The quantities `table` reads for `body_class`: run it on a recording subject, adding every quantity it reports missing, until it
    runs (or fails some other way); what was asked on the way is the list. A quantity read only on another path is added when a
    sample asks for it (`run_samples`)."""
    Recording = _recording_subject_class()
    rng = random.Random(seed_for("asked", table.name, body_class))
    values = {}
    for _ in range(120):
        subject = Recording(body_class, dict(values))
        try:
            table.probabilities(subject)
        except KeyError as error:
            match = re.search(r"no quantity '([^']*)'", str(error))
            if match:
                values[match.group(1)] = _pick(table, match.group(1), rng, registry)
                continue
        except Exception:    # an empty gate, a sum over 1 or a weight that raises: what it asked is still the list
            pass
        return sorted(name for name in set(values) | set(subject.asked) if not name.startswith("const:"))
    return sorted(name for name in values if not name.startswith("const:"))


@dataclass
class Outcome:
    """One sample: the subject's class and quantities, and `kind`: `ok` (with `probabilities`), `gate` (no row supports it),
    `weights` (every weight is zero, or the sourced weights go over 1) or `error` (a weight raised)."""

    body_class: str
    quantities: dict
    kind: str
    detail: str = ""
    probabilities: tuple = ()


def run_samples(table, body_class, registry, count=SAMPLES, targets=None):
    """`count` seeded subjects of `body_class` run through `table.probabilities`, as `Outcome`s, after two for each row that supports the
    class that satisfy the row's own conditions (so a row with a narrow gate is reached)."""
    from planetgen.generation import possibility_table as pt

    rng = random.Random(seed_for("samples", table.name, body_class))
    names = _asked(table, body_class, registry)
    rows = [row for row in table.rows if pt.ANY_CLASS in row.classes or body_class in row.classes]
    plans = list(targets) if targets is not None else [row for row in rows for _ in range(2)] + [None] * count
    outcomes = []
    for target in plans:
        quantities, forced = {}, set()
        for name in names:
            value = _satisfying(target, name, rng, body_class) if target is not None else None
            if value is None:
                value = _pick(table, name, rng, registry, body_class)
            else:
                forced.add(name)
            quantities[name] = value
        for condition in (target.conditions if target is not None else ()):
            if condition.quantity not in quantities:
                quantities[condition.quantity] = _satisfying(target, condition.quantity, rng, body_class)
                forced.add(condition.quantity)
        _couple(quantities, forced, rng)
        for _attempt in range(30):
            try:
                found = table.probabilities(pt.Subject(body_class, quantities))
                outcomes.append(Outcome(body_class, dict(quantities), "ok", probabilities=tuple((row.name, p) for row, p in found)))
                break
            except KeyError as error:
                match = re.search(r"no quantity '([^']*)'", str(error))
                if not match:
                    outcomes.append(Outcome(body_class, quantities, "error", f"KeyError {error}"))
                    break
                quantities[match.group(1)] = _pick(table, match.group(1), rng, registry, body_class)
                _couple(quantities, forced, rng)
            except pt.PossibilityError as error:
                reason = str(error)
                kind = "gate" if "no row supports the class" in reason else "weights"
                outcomes.append(Outcome(body_class, dict(quantities), kind, reason))
                break
            except Exception as error:    # a weight that raises is a finding, whatever it raised
                outcomes.append(Outcome(body_class, dict(quantities), "error", f"{type(error).__name__}: {error}"))
                break
        else:
            outcomes.append(Outcome(body_class, quantities, "error", "the quantities it needs never ended"))
    return outcomes


def table_classes(table):
    """The classes `table` names; a table whose rows all take any class is run for one class that no row names."""
    classes = table.classes()
    return classes or ["any"]


# --- The checks --------------------------------------------------------------------------------------------------------------------


def check_row_states(table):
    from planetgen.generation import possibility_table as pt

    found = []
    for row in table.rows:
        if not str(row.source).strip():
            found.append(f"{table.name}/{row.name}: no source")
        if row.state not in pt.STATES:
            found.append(f"{table.name}/{row.name}: state {row.state!r}")
        if (row.state == pt.OPEN_STATE) != (row.weight == pt.OPEN):
            found.append(f"{table.name}/{row.name}: state {row.state!r} with weight {row.weight!r}")
        if row.state == pt.PARTLY_VERIFIED and not row.limit.strip():
            found.append(f"{table.name}/{row.name}: partly verified with no limit")
        if row.state == pt.STAND_IN and not row.replacement.strip():
            found.append(f"{table.name}/{row.name}: a stand-in with no replacement item")
        if not row.classes:
            found.append(f"{table.name}/{row.name}: an empty class gate")
    return found


def zero_allowed(notes, table_name, body_class):
    entry = notes.get("all_zero_allowed", {}).get(table_name)
    return bool(entry) and ("*" in entry["classes"] or body_class in entry["classes"]) and bool(str(entry.get("note", "")).strip())


def check_samples(table, outcomes_by_class, notes=None):
    """The three subject rules over the outcomes of one table: gate not empty, weights sum to one, no failed subject. A subject whose
    weights are all zero is not a failure for a class `notes["all_zero_allowed"]` lists for the table (with a note saying why)."""
    notes = notes or {}
    from planetgen.generation import possibility_table as pt

    gate, sums, failed = [], [], []
    for body_class, outcomes in outcomes_by_class.items():
        if body_class != "any" and not any(outcome.kind == "ok" for outcome in outcomes):
            gate.append(f"{table.name}/{body_class}: none of {len(outcomes)} subjects reached a row ({collections.Counter(o.kind for o in outcomes)})")
        for outcome in outcomes:
            if outcome.kind == "ok":
                values = [p for _name, p in outcome.probabilities]
                if not values or any(not 0.0 < p <= 1.0 + pt.SUM_TOLERANCE for p in values) or abs(math.fsum(values) - 1.0) > pt.SUM_TOLERANCE:
                    sums.append(f"{table.name}/{body_class}: probabilities {values} add up to {math.fsum(values)!r}")
            elif outcome.kind in ("weights", "error"):
                if outcome.kind == "weights" and "every weight is zero" in outcome.detail and zero_allowed(notes, table.name, body_class):
                    continue
                failed.append(f"{table.name}/{body_class}: {outcome.detail} at {outcome.quantities}")
    return gate, sums, failed


EXTRA_SAMPLES = 400
"""int: More subjects tried for a row that no sample of a class reached, before it is reported unreached."""


def unreached(table, outcomes_by_class, registry):
    """`["table/row/class"]` for each row and class it names that no subject gave a probability, after `EXTRA_SAMPLES` more tries
    that satisfy the row's own conditions. A row that names no class is checked on the class `any`."""
    from planetgen.generation import possibility_table as pt

    reached = collections.defaultdict(set)
    for body_class, outcomes in outcomes_by_class.items():
        for outcome in outcomes:
            for name, _p in outcome.probabilities:
                reached[name].add(body_class)
    found = []
    for row in table.rows:
        classes = [c for c in row.classes if c != pt.ANY_CLASS] or ["any"]
        for body_class in classes:
            if body_class in reached[row.name] or (body_class == "any" and reached[row.name]):
                continue
            more = run_samples(table, body_class, registry, count=0, targets=(row,) * EXTRA_SAMPLES)
            if not any(name == row.name for outcome in more for name, _p in outcome.probabilities):
                found.append(f"{table.name}/{row.name}/{body_class}")
    return found


def check_fractions(table, registry):
    from planetgen.generation import possibility_table as pt
    from planetgen.util import draw

    found = []
    shared = {row.name for row in registry.get(pt.SHARED_COMPONENTS).rows} if table.components and pt.SHARED_COMPONENTS in registry else None
    for row in table.rows:
        if not row.fractions:
            continue
        draw.set_run_seed(seed_for("fractions", table.name, row.name))
        for _ in range(FRACTION_DRAWS):
            drawn = row.draw_fractions()
            total = math.fsum(drawn.values())
            if abs(total - 1.0) > pt.SUM_TOLERANCE or any(value < 0.0 for value in drawn.values()):
                found.append(f"{table.name}/{row.name}: a drawn set of fractions adds up to {total!r}")
                break
        if shared is not None and not set(row.fractions) <= shared:
            found.append(f"{table.name}/{row.name}: components not in the shared table: {sorted(set(row.fractions) - shared)}")
    return found


def check_draw_follows_weights(table, outcomes_by_class, draws=DRAWS, subjects_per_table=2):
    """For up to `subjects_per_table` subjects with two or more rows, `Table.draw` picks each row as often as its probability says."""
    from planetgen.generation import possibility_table as pt
    from planetgen.util import draw

    found, tried = [], 0
    for body_class, outcomes in outcomes_by_class.items():
        for outcome in outcomes:
            if outcome.kind != "ok" or len(outcome.probabilities) < 2 or tried >= subjects_per_table:
                continue
            tried += 1
            draw.set_run_seed(seed_for("draws", table.name, body_class, tried))
            counts = collections.Counter(table.draw(body_class, outcome.quantities).row.name for _ in range(draws))
            for name, probability in outcome.probabilities:
                sigma = math.sqrt(draws * probability * (1.0 - probability)) or 1.0
                if abs(counts.get(name, 0) - draws * probability) > 5.0 * sigma + 1.0:
                    found.append(f"{table.name}/{body_class}/{name}: drawn {counts.get(name, 0)} of {draws}, probability {probability:.4f}")
            unknown = set(counts) - {name for name, _p in outcome.probabilities}
            if unknown:
                found.append(f"{table.name}/{body_class}: drew rows outside the gate: {sorted(unknown)}")
    return found


def open_rows(registry):
    """`{table: [open row names]}` over the registry."""
    from planetgen.generation import possibility_table as pt

    return {table.name: [row.name for row in table.rows if row.state == pt.OPEN_STATE]
            for table in registry.tables() if any(row.state == pt.OPEN_STATE for row in table.rows)}


def load_notes(path=NOTES_FILE):
    return json.loads(pathlib.Path(path).read_text(encoding="utf-8"))


def check_open_rows(registry, notes):
    """Returns `(new, stale)`: open rows with no entry in `notes["open_rows"]` (`{table: {row: note}}`) or with an empty note, and entries
    that no longer match an open row (report only)."""
    found = open_rows(registry)
    listed = notes.get("open_rows", {})
    new = [f"{table}/{row}" for table, names in found.items() for row in names if row not in listed.get(table, {})]
    missing_note = [f"{table}/{row}" for table, rows in listed.items() for row, note in rows.items() if not str(note).strip()]
    stale = [f"{table}/{row}" for table, rows in listed.items() for row in rows if row not in found.get(table, [])]
    return new + [f"{item} (no note)" for item in missing_note], stale


def check_unreached(found, notes):
    """Returns `(new, stale)`: unreached rows (`unreached`) with no note in `notes["unreached_allowed"]`, and notes for rows that are
    reached now (report only)."""
    listed = notes.get("unreached_allowed", {})
    new = [item for item in found if not str(listed.get(item, "")).strip()]
    stale = [item for item in listed if item not in found]
    return new, stale


def unregistered(live_tables, registry):
    """The names of the tables in `live_tables` that are not the registered table of their name."""
    return sorted({table.name for table in live_tables if table.name not in registry or registry.get(table.name) is not table})


def live_tables():
    from planetgen.generation import possibility_table as pt

    return [obj for obj in gc.get_objects() if isinstance(obj, pt.Table)]


def unregistered_in_fresh_process():
    """Starts a process, imports the modules that use the loader, and returns the names of the tables it holds that are not registered."""
    code = ("import json; from tests.audit import sweep_tables as s; r = s.load_tables(); "
            "print(json.dumps(s.unregistered(s.live_tables(), r)))")
    run = subprocess.run([sys.executable, "-I", "-c", f"import sys; sys.path.insert(0, {str(pathlib.Path(__file__).resolve().parents[2])!r}); {code}"],
                         capture_output=True, text=True, timeout=600, check=False)
    if run.returncode:
        raise RuntimeError(f"the fresh process failed: {run.stderr[-800:]}")
    return json.loads(run.stdout.strip().splitlines()[-1])
