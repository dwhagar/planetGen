# planetgen/galaxy/dials.py

"""
The lever table (GEN.601; docs/design/galaxy-size-levers.md): the one definition of the generator's galaxy shape inputs.

Every input the user can give the shape of a galaxy is one `Lever` row here: its name, command-line flag, label, unit,
range, step, what Auto means, the law it moves and the source that sets its range. The `planetgen plan` flags
(`add_arguments`), the Generate page's fields (`web/generate_page.py`, `templates/generate.html`), the out-of-range
message (`Lever.problem`), the help text and the table in `docs/cli.md` (`docs_table`) are all built from `LEVERS`, and
`tests/test_galaxy_dials.py` fails when a flag, a field, a help line or a docs row is missing for a row, or exists
without one.

**Auto.** Every row's default is Auto: the flag is left out (`None`) and the galaxy's seed draws the value through the
law the row moves (`galaxy.shape_draw`). A given value replaces its own draw and what follows from it is recomputed. A
row whose Auto is a fixed number rather than a draw says so in `auto_value`.

**Range.** A given value outside `minimum` to `maximum` is refused with the allowed range and its source. Every range
carries a `status` saying how well its source is known (the labels of docs/design/gap-filling-plan.md): `verified`
(read in a paper or file in the repository), `derived` (calculated from a verified source), `bracket` (a range around
values from sources that disagree or were read second hand), `theory` or `placeholder` (no source yet; GEN.609 replaces
these as sources arrive).
"""

import argparse
import dataclasses
import math
import re

DIALS = (
    ("shape", "Shape"),
    ("scale", "Scale"),
    ("structure", "Structure"),
    ("edge", "Edge"),
    ("halo", "Halo"),
    ("census", "Census"),
    ("pattern", "Pattern"),
)
"""tuple: `(key, label)` of each group a row belongs to, in the order the Generate page shows them. Shape is the
selector; Scale, Structure, Edge, Halo and Census are the five dials; Pattern is the spiral arms."""

STATUSES = ("verified", "derived", "bracket", "theory", "placeholder")
"""tuple: The labels a range's source can carry."""

KINDS = ("float", "int", "enum", "length")


@dataclasses.dataclass(frozen=True)
class Unit:
    """One unit a `length` row accepts. An absolute unit has a `factor` to the row's base unit (`Lever.unit`) and takes
    the row's range; a relative one (`factor` is `None`, such as the scale length `L`) is kept as it was given and
    carries its own `minimum` and `maximum`."""

    name: str
    factor: float = None
    minimum: float = None
    maximum: float = None


class RelativeLength(tuple):
    """`(value, unit name)` of a length given in a relative unit, such as `12L`; it resolves against the galaxy being
    planned."""

    __slots__ = ()

    def __new__(cls, value, unit):
        return super().__new__(cls, (value, unit))

    value = property(lambda self: self[0])
    unit = property(lambda self: self[1])


_LENGTH_PATTERN = re.compile(r"^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*([A-Za-z]+)\s*$")


@dataclasses.dataclass(frozen=True)
class Lever:
    """One row of the table. `minimum`, `maximum` and `step` are `None` for an enum, whose `choices` are its range."""

    name: str
    flag: str
    dial: str
    label: str
    kind: str
    unit: str
    minimum: float
    maximum: float
    step: float
    law: str
    auto_text: str
    help: str
    source: str
    status: str
    choices: tuple = ()
    units: tuple = ()
    auto_value: float = None
    nonzero: bool = False
    advanced: bool = False

    @property
    def dest(self):
        """str: The name the parsed value has on `argparse.Namespace`."""
        return self.flag.lstrip("-").replace("-", "_")

    def range_text(self):
        """The allowed values in words: `1 to 100000 pc` or `spiral`."""
        if self.kind == "enum":
            return ", ".join(self.choices)
        text = f"{self.minimum:g} to {self.maximum:g}"
        text += f" {self.unit}" if self.unit else ""
        for unit in self.units:
            if unit.factor is None:
                text += f", or {unit.minimum:g} to {unit.maximum:g} {unit.name}"
        return text

    def problem(self, value):
        """Why `value` is refused, naming the allowed range and its source, or `None` when it is allowed."""
        if self.kind == "enum":
            if value in self.choices:
                return None
            return f"{value!r} is not one of {self.range_text()} ({self.source}; {self.status})"
        if isinstance(value, RelativeLength):
            unit = next((u for u in self.units if u.name == value.unit and u.factor is None), None)
            if unit is None:
                return f"{value.unit!r} is not one of the units {self.unit_names()}"
            if not unit.minimum <= value.value <= unit.maximum:
                return (f"{value.value:g} {unit.name} is outside {unit.minimum:g} to {unit.maximum:g} {unit.name} "
                        f"({self.source}; {self.status})")
            return None
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            return f"{value!r} is not a finite number"
        if not self.minimum <= value <= self.maximum or (self.nonzero and value == 0):
            return (f"{value:g} is outside {self.range_text()}{', not 0' if self.nonzero else ''} "
                    f"({self.source}; {self.status})")
        return None

    def unit_names(self):
        """The units a `length` row takes, as `L, pc, kpc`."""
        return ", ".join(unit.name for unit in self.units)

    def convert(self, raw):
        """`raw` (text) as this row's kind: a float, an int, the choice, or for a `length` row a float in the base
        unit (or a `RelativeLength`) from a number with its unit, `2600pc`.

        Raises:
            ValueError: If `raw` is not that kind of value; a length without a unit is refused naming the units.
        """
        if self.kind == "enum":
            return raw.strip().lower()
        if self.kind == "length":
            match = _LENGTH_PATTERN.match(raw)
            if not match:
                raise ValueError(f"a length needs its unit: {self.unit_names()}")
            unit = next((u for u in self.units if u.name == match.group(2)), None)
            if unit is None:
                raise ValueError(f"{match.group(2)!r} is not one of the units {self.unit_names()}")
            number = float(match.group(1))
            if not math.isfinite(number):
                raise ValueError("not finite")
            return RelativeLength(number, unit.name) if unit.factor is None else number * unit.factor
        if self.kind == "int":
            return int(raw)
        value = float(raw)
        if not math.isfinite(value):
            raise ValueError("not finite")
        return value

    def parse(self, raw):
        """The `argparse` type of the row's flag: converts and range-checks `raw`."""
        try:
            value = self.convert(raw)
        except ValueError as exc:
            kind = {"int": "a whole number", "float": "a number", "enum": "text", "length": "a length"}[self.kind]
            raise argparse.ArgumentTypeError(
                f"{raw!r} is not {kind} ({exc}); allowed: {self.range_text()}") from None
        message = self.problem(value)
        if message:
            raise argparse.ArgumentTypeError(message)
        return value

    def help_text(self):
        """The row's `--help` line: what it does, the range with its source, and what Auto is."""
        advanced = " Advanced." if self.advanced else ""
        return (f"{self.help} Allowed: {self.range_text()}{', not 0' if self.nonzero else ''} "
                f"({self.source}; {self.status}). Default Auto: {self.auto_text}.{advanced}")


LEVERS = (
    Lever(
        name="shape", flag="--shape", dial="shape", label="Shape", kind="enum", unit="", minimum=None, maximum=None,
        step=None, choices=("spiral",),
        law="the galaxy family: which laws draw the rest of the shape (galaxy.shape_draw)",
        auto_text="the family the seed draws, which is spiral, the only family built",
        help="The galaxy family. A selector, not a dial: it never sets a size or a count.",
        source="the galaxy families; only the spiral exists, later families fill the list", status="theory"),
    Lever(
        name="disk_scale_length", flag="--disk-scale-length", dial="scale", label="Disk scale length",
        kind="length", unit="pc", minimum=0.9, maximum=8000.0, step=0.1, advanced=True,
        units=(Unit("pc", 1.0), Unit("kpc", 1000.0)),
        law="the thin disk's scale length L; every other length of the galaxy is a ratio of it",
        auto_text="drawn from the size-mass relation of the galaxy's stellar mass",
        help="Thin disk exponential scale length with its unit, such as 2600pc or 2.6kpc (the thick disk's is 0.77x).",
        source="Boss, 10-21 (7:04 and 7:08 PM EDT): floor half of Ursa Major III / UNIONS 1 (recalled half-light "
               "radius 3 pc, so 1.5 pc, over 1.678 for a scale length), ceiling 1.5 x Andromeda (recalled disc "
               "scale length 5.3 kpc); neither is read in a file of the repository", status="placeholder"),
    Lever(
        name="flatness", flag="--flatness", dial="structure", label="Flatness (scale length over scale height)",
        kind="float", unit="", minimum=3.0, maximum=20.0, step=0.1,
        law="the thin disk's thickness: its scale height is the scale length over this (the thick disk's is 3x); a "
            "ratio keeps its meaning when the scale length moves",
        auto_text="drawn from the seed, median 7.3 with a 0.13 dex scatter (Kregel et al. 2002)",
        help="How flat the thin disk is: its scale length over its scale height. Larger is thinner.",
        source="median 7.3 and 0.13 dex read in Kregel, van der Kruit and de Grijs 2002 (shape_draw); the range 3 to "
               "20 is the proposal's suggestion (galaxy size levers, 10-21), not read in a source", status="bracket"),
    Lever(
        name="bulge_share", flag="--bulge-share", dial="structure", label="Bulge share (bulge-to-total mass)",
        kind="float", unit="", minimum=0.0, maximum=0.9, step=0.01,
        law="the share of the stellar mass in the bulge and bar; the bulge's central density is solved from it, and "
            "0 leaves a disc with no bulge",
        auto_text="drawn from the seed by Hubble stage (about 0.25 early, 0.06 late; Graham and Worley 2008)",
        help="The bulge-to-total mass ratio of the bulge, thin disk and thick disk.",
        source="stage medians from Graham and Worley 2008 Table 4 (shape_draw); the range 0 to 0.9 is the "
               "proposal's suggestion (galaxy size levers, 10-21), not read in a source", status="bracket"),
    Lever(
        name="bulge_size", flag="--bulge-size", dial="structure", label="Bulge size (scale radius over scale length)",
        kind="float", unit="", minimum=0.02, maximum=1.0, step=0.01,
        law="the bar bulge's scale length along the bar over the disk's scale length (the bar ends at three of them)",
        auto_text="drawn from the seed, a third of a bar of 0.2 to 3 scale lengths (Erwin 2005)",
        help="The bar bulge's scale radius as a fraction of the thin disk's scale length (0.39x across the bar, "
             "0.27x vertically).",
        source="bar length from Erwin 2005 Table 5 and Gadotti 2011 (shape_draw), which draws 0.067 to 1; the "
               "floor of 0.02 is the proposal's suggestion (galaxy size levers, 10-21)", status="bracket"),
    Lever(
        name="nuclear_amplitude", flag="--nuclear-amplitude", dial="structure", label="Nuclear structures",
        kind="float", unit="", minimum=0.0, maximum=10.0, step=0.1, auto_value=1.0, advanced=True,
        law="the nuclear star cluster and nuclear stellar disc, as a multiple of the ones the galaxy's own laws draw "
            "(galaxy.structure_draw)",
        auto_text="1, the galaxy's drawn pair",
        help="How much of the galaxy's own nuclear star cluster and nuclear stellar disc the centre holds, scaled to "
             "the bulge: 1 is the pair its stellar mass draws, 0 is none.",
        source="hand option of the first plan, no source", status="placeholder"),
    Lever(
        name="arm_count", flag="--arm-count", dial="pattern", label="Spiral arms", kind="int", unit="",
        minimum=1, maximum=10, step=1,
        law="the number of spiral arms",
        auto_text="drawn from the arm-count shares (two arms in about two thirds of spirals)",
        help="Number of spiral arms.",
        source="Hart et al. 2017 counts 1 to 5 or more arms; the ceiling of 10 is a placeholder",
        status="placeholder"),
    Lever(
        name="pitch_angle_deg", flag="--pitch-angle-deg", dial="pattern", label="Arm pitch angle (degrees)",
        kind="float", unit="degrees", minimum=-90.0, maximum=90.0, step=0.5, nonzero=True,
        law="the arms' pitch angle (a negative one winds the other way)",
        auto_text="drawn from the seed, tighter for heavier galaxies",
        help="Spiral arm pitch angle, degrees.",
        source="Hart 2017 draws 10 to 40; the full quarter turn is a placeholder", status="placeholder"),
    Lever(
        name="arm_density", flag="--arm-density", dial="pattern", label="Arm density", kind="float", unit="",
        minimum=0.01, maximum=10.0, step=0.1, auto_value=1.4,
        law="the thin disk's density on an arm crest; retired for arm contrast",
        auto_text="1.4",
        help="The thin disk's density on a spiral arm crest, in the disk's own units. At least --interarm-density.",
        source="hand option of the first plan, no source", status="placeholder"),
    Lever(
        name="interarm_density", flag="--interarm-density", dial="pattern", label="Inter-arm density",
        kind="float", unit="", minimum=0.0, maximum=10.0, step=0.1, auto_value=0.6,
        law="the thin disk's density midway between arms; retired for arm contrast",
        auto_text="0.6",
        help="The thin disk's density midway between arms, in the same units.",
        source="hand option of the first plan, no source", status="placeholder"),
)
"""tuple: Every lever, in the order the page and the docs show them."""


def lever(name):
    """The row named `name` (or whose flag is `name`).

    Raises:
        KeyError: If there is none.
    """
    for row in LEVERS:
        if name in (row.name, row.flag):
            return row
    raise KeyError(name)


def add_arguments(group):
    """Adds the flag of every lever to `group` (an argument group or parser); each takes its row's range check and
    defaults to `None`, which is Auto."""
    for row in LEVERS:
        group.add_argument(row.flag, type=row.parse, default=None, metavar=_metavar(row), help=row.help_text())


def _metavar(row):
    if row.kind == "enum":
        return "{" + ",".join(row.choices) + "}"
    return {"int": "N", "float": "X", "length": "LENGTH"}[row.kind]


def resolved(args, name):
    """The value of lever `name` in parsed `args`: the given value, else the row's fixed Auto value (`None` for a row
    the seed draws)."""
    row = lever(name)
    value = getattr(args, row.dest, None)
    return row.auto_value if value is None else value


def docs_table():
    """The markdown table of every lever that `docs/cli.md` carries between its dials markers."""
    lines = ["| Dial | Flag | What it sets | Allowed | Default (Auto) | Source of the range |",
             "|---|---|---|---|---|---|"]
    labels = dict(DIALS)
    for row in LEVERS:
        flag = f"`{row.flag}`" + (" (advanced)" if row.advanced else "")
        allowed = row.range_text() + (", not 0" if row.nonzero else "")
        lines.append(f"| {labels[row.dial]} | {flag} | {row.law} | {allowed} | {row.auto_text} | "
                     f"{row.source} ({row.status}) |")
    return "\n".join(lines)


DOCS_BEGIN = "<!-- dials:begin (generated by planetgen.galaxy.dials; python3 -m planetgen.galaxy.dials --write) -->"
DOCS_END = "<!-- dials:end -->"


def docs_block():
    """The text between (and with) the markers in `docs/cli.md`."""
    return f"{DOCS_BEGIN}\n{docs_table()}\n{DOCS_END}"


def _main(argv=None):
    import pathlib
    import re
    parser = argparse.ArgumentParser(description="Print the lever table, or rewrite its block in docs/cli.md.")
    parser.add_argument("--write", action="store_true", help="Rewrite the block in docs/cli.md.")
    args = parser.parse_args(argv)
    if not args.write:
        print(docs_block())
        return 0
    path = pathlib.Path(__file__).resolve().parents[3] / "docs" / "cli.md"
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(re.escape(DOCS_BEGIN.split(" (")[0]) + r".*?" + re.escape(DOCS_END), re.S)
    if not pattern.search(text):
        raise SystemExit(f"{path} has no dials block")
    path.write_text(pattern.sub(lambda _match: docs_block(), text), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
