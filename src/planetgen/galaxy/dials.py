# planetgen/galaxy/dials.py

"""
The lever table (GEN.601; docs/design/galaxy-size-levers.md): the one definition of the generator's galaxy shape inputs.

Every input the user can give the shape of a galaxy is one `Lever` row here: its name, command-line flag, label, unit,
range, step, what Auto means, the law it moves and the source that sets its range. The `planetgen plan` flags
(`add_arguments`), the Generate page's fields (`web/generate_page.py`, `templates/generate.html`), the out-of-range
message (`Lever.problem`), the help text and the table in `docs/reference/cli.md` (`docs_table`) are all built from `LEVERS`, and
`tests/test_galaxy_dials.py` fails when a flag, a field, a help line or a docs row is missing for a row, or exists
without one.

**Auto.** Every row's default is Auto: the flag is left out (`None`) and the galaxy's seed draws the value through the
law the row moves (`galaxy.shape_draw`). A given value replaces its own draw and what follows from it is recomputed. A
row whose Auto is a fixed number rather than a draw says so in `auto_value`.

**Range.** A given value outside `minimum` to `maximum` is refused with the allowed range and its source. Every range
carries a `status` saying how well its source is known (the labels of docs/design/research/gap-filling-plan.md): `verified`
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

TIER_OF_STATUS = {"verified": "T1", "derived": "T2", "bracket": "T3", "theory": "T5", "placeholder": "T4"}
"""dict: The source tier (Boss's nine, docs/design/stand-in-values.md) each label gives a range; GEN.725. A range whose label is
`placeholder` has no source at all (T4, unverified); the tier is the owner's to change when the range is read."""

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
    min_abs: float = None
    """float: The closest to 0 a value may be (pitch angle: a quarter turn each way, but not a flat line), or `None`."""
    advanced: bool = False
    log: bool = False
    """bool: Whether a slider for the row runs on a logarithmic scale (a mass, a density)."""
    tweak: str = "L0"
    """str: How far the range was moved to fit outcomes (GEN.725): L0 as first envisioned from its source, up to L4 fitted to an
    outcome with no independent source."""

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
        if not self.minimum <= value <= self.maximum or (self.min_abs is not None and abs(value) < self.min_abs):
            return (f"{value:g} is outside {self.range_text()}{self._abs_text()} "
                    f"({self.source}; {self.status})")
        return None

    def _abs_text(self):
        return f", at least {self.min_abs:g} away from 0" if self.min_abs is not None else ""

    def units_text(self):
        """How to type a `length` row, for the text under its box: the number with its unit straight after it, with
        an example and what each unit is."""
        examples = {"L": "12L", "pc": "160pc", "kpc": "40kpc", "Mpc": "0.1Mpc"}
        names = {"pc": "parsecs", "kpc": "kiloparsecs", "Mpc": "megaparsecs"}
        shown = [examples.get(unit.name, f"1{unit.name}") for unit in self.units]
        meanings = []
        for unit in self.units:
            if unit.factor is None:
                meanings.append(f"{unit.name} is scale lengths of the thin disc ({unit.minimum:g} to {unit.maximum:g})")
            else:
                meanings.append(f"{unit.name} is {names.get(unit.name, unit.name)}")
        return (f"Type the number with its unit straight after it, such as {', '.join(shown[:-1])} or {shown[-1]}. "
                f"{'; '.join(meanings)}.")

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

    def format(self, value):
        """`value` as it is typed: `12L` for a length in scale lengths, `160pc` for one in the base unit, else `g`."""
        if isinstance(value, RelativeLength):
            return f"{value.value:g}{value.unit}"
        if self.kind == "length":
            return f"{value:g}{self.unit}"
        return f"{value:g}" if isinstance(value, (int, float)) else str(value)

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
        return (f"{self.help} Allowed: {self.range_text()}{self._abs_text()} "
                f"({self.source}; {self.status}). Default Auto: {self.auto_text}.{advanced}")


HALO_MIN_HALO_MASS_SOLAR = 1.0e10
"""float: The lowest halo mass, solar masses (M200), our halo functions are used at (Boss, 7:33 PM EDT 10-21: a galaxy
below the lowest size they work for reliably runs none of them). Each law in the chain was fitted from here up:
the [R1200] concentration is quoted for M200 of 1e10 and over ([R4609, R4610, R4611, R4616] say the same of formulas derived
there, [R4611]), [R0319]'s stellar-to-halo fit holds from 10^10.5 and [R3689] agree with
abundance matching down to 5e9 ([R3689], read by the galaxy size research thread). The stellar mass it
corresponds to is `galaxy.mass_chain.HALO_MIN_STELLAR_MASS_SOLAR`, about 4.6e6 solar masses. Status: bracket, since
the bound is read from papers that agree to within half a dex and is not itself printed in one."""

FULL_STRUCTURE_MIN_STELLAR_MASS_SOLAR = 3.0e8
"""float: The stellar mass, solar masses, under which a galaxy is one disc plus a bulge (Boss, 7:58 PM EDT 10-21): no
thick disc, no bar, no stellar halo and no nuclear star cluster or disc. The laws that draw them were fitted to
galaxies of this mass and up: [R4824]'s sample starts at a circular speed of about 50 km/s, which
the baryonic Tully-Fisher law of `structure_draw` puts at 2.9e8 (density-limits report 10-21, section 2.3), and below
it the thick disc outweighs the thin one (906 times at 281 solar masses) and the nuclear disc, 61 times the cluster,
holds 80 percent of the stars in a space smaller than a cell. Status: derived."""

HALO_MIN_STATUS = "bracket"
"""str: The range label of `HALO_MIN_HALO_MASS_SOLAR`."""

LEVERS = (
    Lever(
        name="shape", flag="--shape", dial="shape", label="Shape", kind="enum", unit="", minimum=None, maximum=None,
        step=None, choices=("spiral",),
        law="the galaxy family: which laws draw the rest of the shape (galaxy.shape_draw)",
        auto_text="the family the seed draws, which is spiral, the only family built",
        help="The galaxy family. A selector, not a dial: it never sets a size or a count.",
        source="the galaxy families; only the spiral exists, later families fill the list", status="theory"),
    Lever(
        name="mass", flag="--mass", dial="scale", label="Stellar mass", kind="float",
        unit="solar masses", minimum=281.0, maximum=1.55e11, step=0.01, log=True,
        law="the total mass of the galaxy's stars (not its dark halo, which the Halo dial sets); the size, the halo, "
            "the thick disc, the bulge and the black hole all follow it through their laws; under about 4.6e6 the "
            "galaxy has no dark halo (no halo draw, no halo law runs), and under 3e8 it is one disc plus a bulge",
        auto_text="drawn from the galaxy stellar mass function [R0238] between 3e9 and 6e9 solar masses, a small spiral (Boss, 10-21)",
        help="The total mass of the galaxy's stars, solar masses. Size and density follow it unless pinned; pin any "
             "two of mass, size and density and the third is shown.",
        source="Boss, 10-21 (7:04 to 7:31 PM EDT): floor half the stellar mass of Segue 1, 562 solar masses (Local "
               "Volume Database v1.1.1; [R4098]); ceiling 1.5 x Andromeda's 1.03e11 "
               "[R4084], read by the galaxy size research thread", status="derived"),
    Lever(
        name="compactness", flag="--compactness", dial="scale", label="Compactness", kind="float",
        unit="sigma", minimum=-3.0, maximum=3.0, step=0.1,
        law="the size at the same mass, as standard deviations from the size-mass law: 0 is the median galaxy, a "
            "positive value a larger disc (the scale length is the law's median times e^(sigma x z))",
        auto_text="drawn from the seed (a normal draw of the scatter, 0.20 dex falling to 0.15), drawn again when the shape would raise a grid warning",
        help="How large the disc is for its mass, in standard deviations of the size-mass scatter (negative is more "
             "compact). Sets the same size as --disk-scale-length; give only one of them.",
        source="scatter of the late-type size-mass law, [R4048] (read in the repository); three sigma is "
               "the draw's own cap", status="derived"),
    Lever(
        name="density", flag="--density", dial="scale", label="Central surface density",
        kind="float", unit="solar masses per pc^2", minimum=0.04, maximum=1.0e5, step=0.01, log=True,
        law="the thin disc's central surface density, mass over 2 pi L^2: the same thing as the size at a fixed mass",
        auto_text="a read-out of the mass and size",
        help="The thin disc's central surface density, solar masses per square parsec. Pin it with the mass to set "
             "the size, or with the size to set the mass.",
        source="density-limits report (10-21): the smallest density with a planned sector at a 4 pc scale height is "
               "0.043; the largest found in dense stellar systems is 1e5 [R2001]; the "
               "floor also grows with the scale height (check_outline)", status="derived"),
    Lever(
        name="disk_scale_length", flag="--disk-scale-length", dial="scale", label="Disk scale length",
        kind="length", unit="pc", minimum=28.0, maximum=8250.0, step=0.1, advanced=True,
        units=(Unit("pc", 1.0), Unit("kpc", 1000.0)),
        law="the thin disk's scale length L; every other length of the galaxy is a ratio of it",
        auto_text="drawn from the size-mass relation of the galaxy's stellar mass, drawn again (with the flatness) when the shape would raise a grid warning",
        help="Thin disk exponential scale length (the thick disk's is 0.77x).",
        source="Boss, 10-21 (7:04 PM to 8:08 PM EDT): floor 28 pc (Segue 1's half, 7.1 pc, is under the 4 pc cell grid's "
               "reach: a variable sector edge, a later item, would bring it back); ceiling 1.5 x Andromeda's 5.5 kpc "
               "[R0909, R1145], read by the galaxy size "
               "research thread", status="derived"),
    Lever(
        name="flatness", flag="--flatness", dial="structure", label="Flatness (scale length over scale height)",
        kind="float", unit="", minimum=3.0, maximum=20.0, step=0.1,
        law="the thin disk's thickness: its scale height is the scale length over this (the thick disk's is 3x); a "
            "ratio keeps its meaning when the scale length moves",
        auto_text="drawn from the seed, median 7.3 with a 0.13 dex scatter [R2484], drawn again with the size when the shape would raise a grid warning",
        help="How flat the thin disk is: its scale length over its scale height. Larger is thinner.",
        source="median 7.3 and 0.13 dex read in Kregel, van [R2484] (shape_draw); the range 3 to "
               "20 is the proposal's suggestion (galaxy size levers, 10-21), not read in a source", status="bracket"),
    Lever(
        name="bulge_share", flag="--bulge-share", dial="structure", label="Bulge share (bulge-to-total mass)",
        kind="float", unit="", minimum=0.0, maximum=0.9, step=0.01,
        law="the share of the stellar mass in the bulge and bar; the bulge's central density is solved from it, and "
            "0 leaves a disc with no bulge",
        auto_text="drawn from the seed by Hubble stage (about 0.25 early, 0.06 late; [R1656])",
        help="The bulge-to-total mass ratio of the bulge, thin disk and thick disk.",
        source="stage medians from [R1656: Table 4] (shape_draw); the range 0 to 0.9 is the "
               "proposal's suggestion (galaxy size levers, 10-21), not read in a source", status="bracket"),
    Lever(
        name="bulge_size", flag="--bulge-size", dial="structure", label="Bulge size (scale radius over scale length)",
        kind="float", unit="", minimum=0.02, maximum=1.0, step=0.01,
        law="the bar bulge's scale length along the bar over the disk's scale length (the bar ends at three of them)",
        auto_text="drawn from the seed, a third of a bar of 0.2 to 3 scale lengths [R1289]",
        help="The bar bulge's scale radius as a fraction of the thin disk's scale length (0.39x across the bar, "
             "0.27x vertically).",
        source="bar length from [R1289: Table 5] and [R1486] (shape_draw), which draws 0.067 to 1; the "
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
        name="edge", flag="--edge", dial="edge", label="Edge", kind="length", unit="pc", minimum=4.0,
        maximum=262140.0, step=1.0,
        units=(Unit("L", None, 3.0, 20.0), Unit("pc", 1.0), Unit("kpc", 1000.0), Unit("Mpc", 1.0e6)),
        law="how far out the plan walks the outline of the galaxy: the outline stops here even where sectors still "
            "clear the one-star threshold (a relative length is in scale lengths of the thin disc)",
        auto_text="the outline's own edge, where no sector clears the one-star-per-sector threshold, up to the grid's "
                  "reach",
        help="How far from the centre the plan walks the outline. Replaces --max-shell.",
        source="relative: Boss, 10-21 (7:01 PM EDT), 3 to 20 scale lengths; absolute: one 4 pc shell up to the grid's "
               "reach of 65,535 shells, 262 kpc (spherical.MAX_ADDRESSABLE_SHELL)", status="derived"),
    Lever(
        name="halo_ratio", flag="--halo-ratio", dial="halo", label="Dark halo mass over stellar mass", kind="float",
        unit="", minimum=10.0, maximum=3000.0, step=1.0, log=True,
        law="the dark halo's mass (M200) as a multiple of the galaxy's stellar mass; the virial radius, the tidal "
            "field and the disc's gas inflow follow it. Refused for a galaxy under the halo threshold",
        auto_text="drawn from the stellar-to-halo relation of the stellar mass [R3181], with its scatter",
        help="The dark halo's mass over the stellar mass. A galaxy under about 4.6e6 solar masses of stars has no halo "
             "and refuses it.",
        source="the relation's own range: a halo holds 26 times its stars at the peak of star formation efficiency "
               "([R3181], read in the repository) and 2200 times at the halo threshold of 1e10 solar masses; "
               "the ends are those, rounded", status="derived"),
    Lever(
        name="stellar_halo_share", flag="--stellar-halo-share", dial="halo", label="Stellar halo share", kind="float",
        unit="", minimum=0.001, maximum=0.3, step=0.001, log=True,
        law="the share of the galaxy's stars in its diffuse stellar halo (the rest are the disc, bulge and bar); "
            "refused for a galaxy under 3e8 solar masses of stars, which has none",
        auto_text="drawn from a power law of the stellar mass (a placeholder, [R1816]) with 1.01 dex "
                  "scatter",
        help="The stellar halo's share of all the galaxy's stars, 0.001 to 0.3 (0.1 to 30 percent). A galaxy under "
             "3e8 solar masses of stars has none and refuses it.",
        source="the draw's own cut: 0.1 percent (Dragonfly galaxies with no detected halo) to 30 percent (the "
               "disc-dominated galaxies of [R0329]), read in structure_draw", status="bracket"),
    Lever(
        name="richness", flag="--richness", dial="census", label="Star richness", kind="float",
        unit="systems per solar mass", minimum=1.0, maximum=3.0, step=0.01, advanced=True,
        law="how many systems each solar mass of the galaxy's stars makes, so the galaxy's number of systems and with it "
            "the density in systems and the cost of a full generation; the mean system mass is its inverse",
        auto_text="the mean system mass the generator derives for the galaxy's populations (about 2 systems per solar "
                  "mass; galaxy.system_mass)",
        help="Systems made per solar mass of stars (the inverse of the mean mass of a system). Auto is the value the "
             "generator derives from its own mass function, remnants and multiplicity.",
        source="the derived mean system mass runs 0.47 (old stars) to 0.79 solar masses (stars just formed), 1.3 to 2.1 "
               "systems per solar mass (the generator's own laws); the band 1 to 3 widens it by about a factor of 1.4", status="bracket"),
    Lever(
        name="arm_count", flag="--arm-count", dial="pattern", label="Spiral arms", kind="int", unit="",
        minimum=1, maximum=10, step=1,
        law="the number of spiral arms",
        auto_text="drawn from the arm-count shares (two arms in about two thirds of spirals)",
        help="Number of spiral arms.",
        source="[R1836] counts 1 to 5 or more arms; the ceiling of 10 is a placeholder",
        status="placeholder"),
    Lever(
        name="pitch_angle_deg", flag="--pitch-angle-deg", dial="pattern", label="Arm pitch angle",
        kind="float", unit="degrees", minimum=-90.0, maximum=90.0, step=0.5, min_abs=1.0,
        law="the arms' pitch angle (a negative one winds the other way)",
        auto_text="drawn from the seed, tighter for heavier galaxies",
        help="How tightly the spiral arms wind; a negative angle winds them the other way.",
        source="[R1836] draws 10 to 40; the full quarter turn is a placeholder; under 1 degree from 0 is refused "
               "because a flat pitch divides by zero (density-limits report, 10-21)", status="placeholder"),
    Lever(
        name="arm_contrast", flag="--arm-contrast", dial="pattern", label="Arm contrast", kind="float", unit="",
        minimum=0.0, maximum=1.0, step=0.01, auto_value=0.4,
        law="how much denser the thin disk is on an arm crest than between arms: the arm factor is 1 + contrast x cos, "
            "so the crest is 1 + contrast and the trough 1 - contrast of the disk's own density (0 is no arms, 1 an "
            "empty trough)",
        auto_text="0.4 (a crest 1.4 and a trough 0.6 of the disk's mean)",
        help="Arm contrast: the thin disk's density on an arm crest over its mean is 1 + this, between arms 1 - this.",
        source="0 to 1 because the density cannot go negative (theory); the 0.4 default is the first plan's hand value, "
               "no source", status="placeholder"),
)
"""tuple: Every lever, in the order the page and the docs show them."""


DERIVED_LIMITS = {
    "disk_scale_height_pc": (
        4.0, 1155.0, "pc",
        "Boss, 10-21 (7:13 to 8:08 PM EDT): floor one 4 pc cell (the density-limits report: a thinner disc is "
        "overcounted by the grid); ceiling 1.5 x Andromeda's 770 pc [R0979], read by "
        "the galaxy size research thread", "derived"),
    "bulge_scale_radius_pc": (
        2.0, 1.0e9, "pc",
        "density-limits report (10-21): half a 4 pc cell; a bulge smaller than that puts five times too many systems "
        "in its cell", "derived"),
}
"""dict: `{quantity: (low, high, unit, source, status)}` for a quantity no row sets but the Structure dial gives (the
scale height is the scale length over the flatness): a value outside it, pinned or solved, is refused with `check`."""


WARN_HEIGHT_PC, WARN_BULGE_PC, WARN_CREST = 8.0, 4.0, 20.0
"""float: The scale height and bulge radius (pc) and the crest density (in units of the one-star density) under which the
4 pc grid represents the galaxy badly; a warning, not a refusal (the scale length floor is 28 pc, where the planned
sectors run 4 to 15 percent high)."""


def shape_warnings(crest, share, height_pc, bulge_radius_pc):
    """The sentences for a shape the 4 pc grid represents badly (the planned sectors expect more systems than the
    model holds), empty for a shape it holds well. `crest` and `share` are the plan's figures for the galaxy
    (`galaxy.plan_figures`)."""
    notes = []
    if height_pc < WARN_HEIGHT_PC:
        notes.append(f"a scale height of {height_pc:.3g} pc is under {WARN_HEIGHT_PC:g}: the disc is under two cells deep")
    if bulge_radius_pc < WARN_BULGE_PC:
        notes.append(f"a bulge radius of {bulge_radius_pc:.3g} pc is under {WARN_BULGE_PC:g}: its cell expects too many "
                     f"systems")
    if crest < WARN_CREST:
        notes.append(f"the disc's crest is {crest:.1f} times the one-star density (under {WARN_CREST:g}): about "
                     f"{100 * share:.0f} percent of its stars sit in planned sectors (approximate)")
    return notes


def check_derived(name, value):
    """Why `value` of the derived quantity `name` is outside its limits (naming the range and source), or `None`."""
    low, high, unit, source, status = DERIVED_LIMITS[name]
    if low <= value <= high:
        return None
    return f"{name.replace('_', ' ')} of {value:.4g} {unit} is outside {low:g} to {high:g} {unit} ({source}; {status})"


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


def restored(value):
    """`value` of a lever as it comes back from a job's JSON parameters: a length in scale lengths is a list there."""
    if isinstance(value, (list, tuple)) and not isinstance(value, RelativeLength):
        return RelativeLength(*value)
    return value


def resolved(args, name):
    """The value of lever `name` in parsed `args`: the given value, else the row's fixed Auto value (`None` for a row
    the seed draws)."""
    row = lever(name)
    value = getattr(args, row.dest, None)
    return row.auto_value if value is None else value


def describe_pins(pins):
    """The sentence the galaxy page shows for the levers a plan was told (DB.81), `pins` being `{lever name: value}` as
    `shape_draw.pins_of` stores it: what was pinned, in the table's order and units, and that the rest was drawn."""
    shown = [f"{row.label} {pins[row.name]:.4g}{' ' + row.unit if row.unit else ''}"
             for row in LEVERS if row.name in (pins or {})]
    if not shown:
        return "Shape drawn from the seed"
    return f"Shape pinned: {', '.join(shown)}; the rest drawn from the seed"


def docs_table():
    """The markdown table of every lever that `docs/reference/cli.md` carries between its dials markers."""
    lines = ["| Dial | Flag | What it sets | Allowed | Default (Auto) | Source of the range |",
             "|---|---|---|---|---|---|"]
    labels = dict(DIALS)
    for row in LEVERS:
        flag = f"`{row.flag}`" + (" (advanced)" if row.advanced else "")
        allowed = row.range_text() + row._abs_text()
        lines.append(f"| {labels[row.dial]} | {flag} | {row.law} | {allowed} | {row.auto_text} | "
                     f"{row.source} ({row.status}) |")
    return "\n".join(lines)


DOCS_BEGIN = "<!-- dials:begin (generated by planetgen.galaxy.dials; python3 -m planetgen.galaxy.dials --write) -->"
DOCS_END = "<!-- dials:end -->"


def halo_note():
    """The sentence under the docs table that says where the halo laws stop."""
    return (f"A galaxy whose halo would weigh under {HALO_MIN_HALO_MASS_SOLAR:.0e} solar masses (a stellar mass under "
            f"about 4.6e6) has no dark halo: none is drawn and no law that reads one runs ({HALO_MIN_STATUS}; see "
            f"`dials.HALO_MIN_HALO_MASS_SOLAR`). A galaxy under {FULL_STRUCTURE_MIN_STELLAR_MASS_SOLAR:.0e} solar masses of "
            f"stars is one disc plus a bulge: no thick disc, bar, stellar halo or nuclear structures "
            f"(`dials.FULL_STRUCTURE_MIN_STELLAR_MASS_SOLAR`).")


def docs_block():
    """The text between (and with) the markers in `docs/reference/cli.md`."""
    return f"{DOCS_BEGIN}\n{docs_table()}\n\n{halo_note()}\n{DOCS_END}"


REGISTRY_BEGIN = "<!-- dial-ranges:begin (generated by planetgen.galaxy.dials; python3 -m planetgen.galaxy.dials --write) -->"
REGISTRY_END = "<!-- dial-ranges:end -->"


def registry_block():
    """The registry rows for the shape dials' ranges in `docs/design/research/gap-filling-plan.md` (GEN.609): one per lever."""
    lines = [REGISTRY_BEGIN, "| Lever | Range | Label | Tier | Tweak | Source |", "|---|---|---|---|---|---|"]
    for row in LEVERS:
        lines.append(f"| `{row.flag}` | {row.range_text()}{row._abs_text()} | {row.status} | {TIER_OF_STATUS[row.status]} | {row.tweak} | {row.source} |")
    lines.append(REGISTRY_END)
    return "\n".join(lines)


def _rewrite(path, begin, end, block):
    import re
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(re.escape(begin.split(" (")[0]) + r".*?" + re.escape(end), re.S)
    if not pattern.search(text):
        raise SystemExit(f"{path} has no {end} block")
    path.write_text(pattern.sub(lambda _match: block, text), encoding="utf-8")


def _main(argv=None):
    import pathlib
    parser = argparse.ArgumentParser(description="Print the lever table, or rewrite its blocks in the docs.")
    parser.add_argument("--write", action="store_true", help="Rewrite the blocks in docs/reference/cli.md and the gap-filling plan.")
    args = parser.parse_args(argv)
    if not args.write:
        print(docs_block())
        return 0
    docs = pathlib.Path(__file__).resolve().parents[3] / "docs"
    _rewrite(docs / "reference" / "cli.md", DOCS_BEGIN, DOCS_END, docs_block())
    _rewrite(docs / "design" / "research" / "gap-filling-plan.md", REGISTRY_BEGIN, REGISTRY_END, registry_block())
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
