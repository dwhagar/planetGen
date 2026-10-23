# planetgen/generation/possibility_table.py

"""
The gated possibility table and its loader (GEN.472)
=====================================================

Every composition category (the gases of an atmosphere, the materials of a body, the abundance pattern of a star, the
origin of a moon) is a table of possibilities. A row names one possibility, the classes it supports, a weight, where the
weight comes from. A body's composition is drawn by taking the rows whose gate passes the
body, evaluating their weights, normalising to 100 percent and picking one. The format and the draw are written up in
docs/design/composition-model.md section 4b; this module is the one loader every table uses.

The draw, in order:

1. `Table.passing`: the rows whose gate passes the subject: its class is in the row's `classes` (or the row says
   `ANY_CLASS`) and every `Condition` on a stored quantity holds.
2. Each weight is evaluated at the subject: a number, or a function of the subject (its stored quantities: [Fe/H], mass,
   age, the formation distance, the snow line). Rows whose weight is zero are dropped.
3. `Table.probabilities` normalises what is left to add up to 1 (within 1e-12) and `Table.draw` picks one with the one
   random wrapper (`util.draw`, GEN.56), in the table's fixed row order.
4. The drawn row supplies the proportions of its components: `Row.draw_fractions` (below). The fractions are
   normalised separately from the weights, so two sums are checked: the weights of the possibilities and the fractions of
   the drawn one.

Open weights. A weight with no source is `OPEN`. The sourced weights of the rows that pass are probabilities; the open rows
of that gate share what is left over equally, so the table still adds up to 100 percent. If the sourced weights already add up
to 1 or more and an open row passes, the open row gets nothing and is named in the report; if the sourced weights add up to more than 1
while an open row passes, the table is refused for that subject (an error), because the open rows' share would be negative.

Fractions inside a row. A row with a sourced kernel (`fractions` with no `concentration`) gives its components' proportions
directly. A row with only a range (`fractions` and a `concentration`) draws them from a Dirichlet distribution whose mean is
the fractions and whose concentration is given (`concentration_from_spread` derives it from a known standard deviation);
both ways the result adds up to 1 within 1e-12.

Rules that follow, each an error and never a silent fallback:

- No row passes the gate, or every weight is zero: `PossibilityError`, naming the class and the table.
- A row without a source, a partly verified row without its limit, a stand-in without its replacement
  item, a row whose state and weight disagree (`OPEN` is the open state and no other), a duplicate row name, or a
  component that is not a row of the one shared components table: `TableError`, at load.
- A row's source state (`VERIFIED`, `PARTLY_VERIFIED`, `OPEN`, `STAND_IN`) is counted per table and listed in
  `Registry.report`.

A table has a `version`; adding, removing or reweighting a row changes the draws of every seed in its classes, so the version
goes up with a generator epoch bump (`src/tests/golden/epochs.json`). `Table.digest` is a hash of everything that decides a draw,
and `Registry.versions` lists `(version, digest)` per table, so a test can pin them.

Tables register themselves (`REGISTRY.register(table)`), so a sweep (TEST.214) walks every table, including ones added later.
This module imports nothing from the generator beyond `util.draw`, so any table module can import it.
"""

import collections
import hashlib
import json
import math
from dataclasses import dataclass, field

from planetgen.util import draw

ANY_CLASS = "*"
"""str: In a row's `classes`: the row supports every class."""

OPEN = "open"
"""str: A weight with no source."""

VERIFIED = "verified"
PARTLY_VERIFIED = "partly verified"
OPEN_STATE = "open"
STAND_IN = "stand-in"
STATES = (VERIFIED, PARTLY_VERIFIED, OPEN_STATE, STAND_IN)
"""tuple: The source states of a row (composition-model.md section 4a item 10, section 4b)."""

SHARED_COMPONENTS = "materials"
"""str: The name of the one materials-and-components table every kind's fractions refer to (section 4a item 8)."""

SUM_TOLERANCE = 1e-12
"""float: How far a normalised sum may be from 1."""


class TableError(ValueError):
    """A table or a row that the loader refuses at load."""


class PossibilityError(LookupError):
    """
    A draw that cannot be made.

    Attributes:
        table (str): The table's name.
        body_class (str): The class of the subject.
    """

    def __init__(self, table, body_class, reason):
        self.table = table
        self.body_class = body_class
        super().__init__(f"table {table!r}, class {body_class!r}: {reason}")


class Subject:
    """
    What a row's gate and weight are evaluated on: a class and the stored quantities of the body or star.

    Attributes:
        body_class (str): The planet, moon, star or object class.
        quantities (dict): Named stored quantities (`fe_h`, `mass_solar`, `age_gy`, `formation_distance_au`,
            `snow_line_au`, `condensation_temperature_k`, ...).
    """

    __slots__ = ("body_class", "quantities")

    def __init__(self, body_class, quantities=None):
        self.body_class = body_class
        self.quantities = dict(quantities or {})

    def __getitem__(self, name):
        try:
            return self.quantities[name]
        except KeyError:
            raise KeyError(f"the subject of class {self.body_class!r} has no quantity {name!r}") from None

    def get(self, name, default=None):
        return self.quantities.get(name, default)


@dataclass(frozen=True)
class Condition:
    """
    A predicate of a row's gate on one stored quantity: it holds when `minimum <= value <= maximum` (a missing bound is
    open), when `equals` is given the value is equal to it, and when `one_of` is given the value is one of its members.
    A quantity the subject lacks is an error and never a pass.
    """

    quantity: str
    minimum: float | None = None
    maximum: float | None = None
    equals: object = None
    one_of: tuple | None = None

    def holds(self, subject):
        value = subject[self.quantity]
        if self.equals is not None and value != self.equals:
            return False
        if self.one_of is not None and value not in self.one_of:
            return False
        if self.minimum is not None and value < self.minimum:
            return False
        if self.maximum is not None and value > self.maximum:
            return False
        return True


@dataclass(frozen=True)
class Row:
    """
    One possibility.

    Attributes:
        name (str): What the row is, named for what it is and never for a Solar System object.
        classes (frozenset or tuple): The class gate: the classes the row supports, or `(ANY_CLASS,)`.
        weight: A number, `OPEN`, or a function of the `Subject`.
        source (str): Where the weight comes from; `"open"` for an open weight.
        state (str): One of `STATES`.
        conditions (tuple): `Condition`s that must also hold.
        limit (str): For a partly verified row, the limit the check named (validity window or caveat).
        replacement (str): For a stand-in row, the research item that replaces it (`"GEN.473"`).
        fractions (dict): The row's components and their proportions, by name in the shared components table.
        concentration (float or None): `None` takes `fractions` directly (a sourced kernel); a number draws them from a
            Dirichlet distribution around them.
        notes (str): The ranges the row's drawn fractions take, the validity range, the check.
        tunable (bool): The row's weight or draw has knobs that are tuned to a published statistic (GEN.503): the
            tuning-target table must hold at least one target for it (`Registry.validate(targets=...)`). It does not decide
            a draw, so it is not part of `Table.digest`.
    """

    name: str
    classes: tuple
    weight: object
    source: str
    state: str = VERIFIED
    conditions: tuple = ()
    limit: str = ""
    replacement: str = ""
    fractions: dict = field(default_factory=dict)
    concentration: float | None = None
    notes: str = ""
    tunable: bool = False

    def __post_init__(self):
        object.__setattr__(self, "classes", tuple(sorted(set(self.classes))))
        object.__setattr__(self, "conditions", tuple(self.conditions))
        object.__setattr__(self, "fractions", dict(self.fractions))
        name = self.name
        if not isinstance(name, str) or not name.strip():
            raise TableError("a row has no name")
        if not self.classes:
            raise TableError(f"row {name!r} lists no supported class (an empty gate)")
        if self.state not in STATES:
            raise TableError(f"row {name!r}: state {self.state!r} is none of {', '.join(STATES)}")
        if not str(self.source).strip():
            raise TableError(f"row {name!r} has no source")
        if self.state == OPEN_STATE:
            if self.weight != OPEN or self.source != OPEN:
                raise TableError(f"row {name!r}: an open row has the weight and the source 'open' and no other")
        else:
            if self.weight == OPEN or self.source == OPEN:
                raise TableError(f"row {name!r}: a weight with no source is written 'open' with the state 'open'")
            if not callable(self.weight):
                if isinstance(self.weight, bool) or not isinstance(self.weight, (int, float)) or \
                        not math.isfinite(self.weight) or self.weight < 0:
                    raise TableError(f"row {name!r}: the weight {self.weight!r} is not a finite number of 0 or more")
        if self.state == PARTLY_VERIFIED and not self.limit.strip():
            raise TableError(f"row {name!r} is partly verified and names no limit")
        if self.state == STAND_IN and not self.replacement.strip():
            raise TableError(f"row {name!r} is a stand-in and names no replacement item")
        if self.fractions:
            bad = {k: v for k, v in self.fractions.items()
                   if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0}
            if bad or not sum(self.fractions.values()) > 0:
                raise TableError(f"row {name!r}: the fractions must be finite, not negative and not all zero")
        if self.concentration is not None:
            if not self.fractions:
                raise TableError(f"row {name!r}: a concentration needs fractions to centre on")
            if not (math.isfinite(self.concentration) and self.concentration > 0):
                raise TableError(f"row {name!r}: the concentration must be a positive number")

    def supports(self, subject):
        """The gate: the class passes and every condition holds."""
        if ANY_CLASS not in self.classes and subject.body_class not in self.classes:
            return False
        return all(condition.holds(subject) for condition in self.conditions)

    def evaluate(self, subject):
        """The weight at `subject`: `OPEN`, or a number of 0 or more."""
        if self.weight == OPEN:
            return OPEN
        value = self.weight(subject) if callable(self.weight) else self.weight
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise TableError(f"row {self.name!r}: the weight evaluated to {value!r}, not a finite number of 0 or more")
        return float(value)

    def normalised_fractions(self):
        """The row's fractions as given, scaled to add up to 1."""
        total = math.fsum(self.fractions.values())
        return {name: value / total for name, value in self.fractions.items()}

    def draw_fractions(self):
        """
        The proportions of this possibility's components for one body, adding up to 1 within `SUM_TOLERANCE`: the
        fractions directly for a sourced kernel, or a Dirichlet draw around them (concentration times the fractions) where
        the row has only a range. The components keep the row's order.
        """
        if not self.fractions:
            raise TableError(f"row {self.name!r} has no fractions to draw")
        mean = self.normalised_fractions()
        if self.concentration is None:
            return mean
        names = list(mean)
        drawn = draw.dirichlet([self.concentration * mean[n] for n in names])
        return dict(zip(names, drawn))


def concentration_from_spread(mean, sd):
    """
    The Dirichlet concentration for one component whose mean fraction is `mean` and whose standard deviation is `sd`:
    variance = mean (1 - mean) / (concentration + 1). Call it for each component with a known spread and take their median
    for the row.
    """
    if not 0.0 < mean < 1.0:
        raise ValueError("the mean must lie strictly between 0 and 1")
    if not sd > 0.0:
        raise ValueError("the spread must be positive")
    value = mean * (1.0 - mean) / (sd * sd) - 1.0
    if value <= 0.0:
        raise ValueError(f"a spread of {sd} is more than a fraction of {mean} can have")
    return value


@dataclass(frozen=True)
class Possibility:
    """A drawn row: the row, its probability among the rows that passed, and all the probabilities."""

    row: Row
    probability: float
    probabilities: tuple

    @property
    def name(self):
        return self.row.name


class Table:
    """
    One category's table of possibilities.

    Args:
        name (str): The table's registered name (`"atmosphere-recipes"`).
        version (int): Goes up with every change to a row (and a generator epoch bump).
        rows (iterable): The `Row`s, in their fixed order.
        components (str or None): For a table whose rows have `fractions`: the shared table their component names refer to
            (`SHARED_COMPONENTS`); `None` for a table with none.
        description (str): What the table decides.
    """

    def __init__(self, name, version, rows, components=None, description=""):
        if not str(name).strip():
            raise TableError("a table has no name")
        if not isinstance(version, int) or isinstance(version, bool) or version < 1:
            raise TableError(f"table {name!r}: the version is a whole number from 1")
        self.name = name
        self.version = version
        self.rows = tuple(rows)
        self.components = components
        self.description = description
        if not self.rows:
            raise TableError(f"table {name!r} has no rows")
        seen = set()
        for row in self.rows:
            if not isinstance(row, Row):
                raise TableError(f"table {name!r}: {row!r} is not a Row")
            if row.name in seen:
                raise TableError(f"table {name!r}: the row {row.name!r} appears twice")
            seen.add(row.name)
            if row.fractions and components is None:
                raise TableError(f"table {name!r}: row {row.name!r} has fractions but the table names no shared "
                                 f"components table; a component list owned by one kind is refused")

    def __len__(self):
        return len(self.rows)

    def row(self, name):
        for row in self.rows:
            if row.name == name:
                return row
        raise KeyError(f"table {self.name!r} has no row {name!r}")

    def classes(self):
        """The classes the rows name, sorted (`ANY_CLASS` left out)."""
        return sorted({c for row in self.rows for c in row.classes if c != ANY_CLASS})

    def passing(self, subject):
        """The rows whose gate passes `subject`, in row order."""
        return [row for row in self.rows if row.supports(subject)]

    def probabilities(self, subject):
        """
        `[(row, probability)]` for the rows that pass `subject`, in row order, with a weight above zero; the
        probabilities add up to 1 within `SUM_TOLERANCE`.

        Raises:
            PossibilityError: No row passes the gate, every weight is zero, or the sourced weights leave a negative share
                for the open rows.
        """
        passed = self.passing(subject)
        if not passed:
            raise PossibilityError(self.name, subject.body_class, "no row supports the class")
        values = [(row, row.evaluate(subject)) for row in passed]
        sourced = [(row, value) for row, value in values if value != OPEN]
        open_rows = [row for row, value in values if value == OPEN]
        total = math.fsum(value for _row, value in sourced)
        if open_rows:
            if total > 1.0 + SUM_TOLERANCE:
                raise PossibilityError(self.name, subject.body_class,
                                       f"the sourced weights add up to {total}, over 1, and open rows "
                                       f"({', '.join(r.name for r in open_rows)}) share the rest")
            share = max(0.0, 1.0 - total) / len(open_rows)
            weights = [(row, share if value == OPEN else value) for row, value in values]
        else:
            weights = sourced
        weights = [(row, w) for row, w in weights if w > 0.0]
        norm = math.fsum(w for _row, w in weights)
        if not weights or not norm > 0.0:
            raise PossibilityError(self.name, subject.body_class, "every weight is zero")
        return [(row, w / norm) for row, w in weights]

    def draw(self, body_class, quantities=None):
        """
        Picks one possibility for a body of `body_class` with the stored `quantities`, with the one random wrapper
        (`util.draw`) in the table's row order. Returns a `Possibility`.
        """
        subject = body_class if isinstance(body_class, Subject) else Subject(body_class, quantities)
        table = self.probabilities(subject)
        rows = [row for row, _p in table]
        chosen = draw.choices(rows, weights=[p for _row, p in table], k=1)[0]
        probability = dict((row.name, p) for row, p in table)[chosen.name]
        return Possibility(chosen, probability, tuple((row.name, p) for row, p in table))

    def state_counts(self):
        """`{state: rows}` over the table, every state present."""
        counts = collections.Counter(row.state for row in self.rows)
        return {state: counts.get(state, 0) for state in STATES}

    def digest(self):
        """A hash of what decides a draw: the version, every row's name, gate, weight (a function by its qualified
        name), state, fractions and concentration. A change to any of them changes it."""
        parts = {"name": self.name, "version": self.version, "components": self.components, "rows": []}
        for row in self.rows:
            weight = row.weight
            if callable(weight):
                weight = f"fn:{getattr(weight, '__module__', '')}.{getattr(weight, '__qualname__', repr(weight))}"
            parts["rows"].append({
                "name": row.name, "classes": list(row.classes), "weight": weight, "state": row.state,
                "conditions": [[c.quantity, c.minimum, c.maximum, c.equals] + ([list(c.one_of)] if c.one_of is not None else [])
                               for c in row.conditions],
                "fractions": sorted(row.fractions.items()), "concentration": row.concentration})
        text = json.dumps(parts, sort_keys=True, default=str)
        return hashlib.sha256(text.encode("utf-8")).hexdigest()


class Registry:
    """The tables, by name, so a sweep can walk every one (TEST.214) and a table added later is covered."""

    def __init__(self):
        self._tables = {}

    def register(self, table):
        """Registers `table`; a second table of the same name is an error. Returns the table."""
        if table.name in self._tables:
            raise TableError(f"a table named {table.name!r} is already registered")
        self._tables[table.name] = table
        return table

    def unregister(self, name):
        del self._tables[name]

    def get(self, name):
        try:
            return self._tables[name]
        except KeyError:
            raise KeyError(f"no table named {name!r}") from None

    def tables(self):
        """Every registered table, in registration order."""
        return list(self._tables.values())

    def __contains__(self, name):
        return name in self._tables

    def __len__(self):
        return len(self._tables)

    def validate(self, shared=SHARED_COMPONENTS, targets=None):
        """
        Checks what needs the whole registry: every component named in a row's fractions is a row of the shared components
        table, and every table with fractions names that table. Raises `TableError` listing the unknown components.

        With `targets` (a `generation.tuning_harness.TargetTable`) it also refuses a tunable row that has no target in it.
        """
        if targets is not None:
            targets.require(self._tables.values())
        for table in self._tables.values():
            if table.components is None:
                continue
            if table.components != shared:
                raise TableError(f"table {table.name!r} takes its components from {table.components!r}, "
                                 f"not from the one shared table {shared!r}")
            if shared not in self._tables:
                raise TableError(f"table {table.name!r} needs the shared components table {shared!r}, "
                                 f"which is not registered")
            known = {row.name for row in self._tables[shared].rows}
            for row in table.rows:
                missing = sorted(set(row.fractions) - known)
                if missing:
                    raise TableError(f"table {table.name!r}, row {row.name!r}: components not in the shared table "
                                     f"{shared!r}: {', '.join(missing)}")

    def versions(self):
        """`{name: (version, digest)}` for every table."""
        return {table.name: (table.version, table.digest()) for table in self._tables.values()}

    def report(self):
        """
        One text: per table the count of rows in each source state, then every open weight and every stand-in with
        the classes it applies to (and its replacement item), so each is found and replaced.
        """
        lines = ["Possibility tables", ""]
        open_lines, stand_in_lines, limited = [], [], []
        for table in self._tables.values():
            counts = table.state_counts()
            lines.append(f"{table.name} (version {table.version}): " + ", ".join(f"{n} {s}" for s, n in counts.items()))
            for row in table.rows:
                where = ", ".join(row.classes)
                if row.state == OPEN_STATE:
                    open_lines.append(f"  {table.name}: {row.name} [{where}]")
                elif row.state == STAND_IN:
                    stand_in_lines.append(f"  {table.name}: {row.name} [{where}] replaced by {row.replacement}")
                elif row.state == PARTLY_VERIFIED:
                    limited.append(f"  {table.name}: {row.name} [{where}] limit: {row.limit}")
        for title, found in (("Open weights", open_lines), ("Stand-ins", stand_in_lines),
                             ("Partly verified, with their limit", limited)):
            lines += ["", f"{title}: {len(found)}"] + found
        return "\n".join(lines) + "\n"


REGISTRY = Registry()
"""Registry: The tables of the generator; each table module registers its own at import."""
