# tests/audit/close_binary.py

"""
The close-binary check (TEST.196; docs/design/generation/placement-by-population.md, section 3.1, "the check"):
the share of primaries with a companion inside 10 AU, by primary-mass bin and [Fe/H] bin, drawn the way the generator
draws it, against the published ranges, with every out-of-range bin flagged.

Boss (2026-10-16, 6:15 PM EDT, on the GEN.103 cards): the metallicity effect on close binaries applies to every star
mass, "but let's build into that a check to see if the numbers act as expected". The ranges are reference values for the
check and are never read by the draw.

A draw is the generator's own: `StarSystem._companion` (the kind, the mass ratio and the [Fe/H] odds factor) and, for a
wide pair, `wide_binary.sample_wide_binary_separation_au` at the limit the place keeps. A primary has a companion inside
10 AU when its pair is close, or wide with an axis under `companion_classes.INNER_LIMIT_AU`. The default place is the
Sun's neighbourhood (the system a command line draws).

What fails (the rows of the design note's table):

- **G dwarfs (1 solar mass), the shape**: the share at [Fe/H] -3.0, -1.0, -0.2 and +0.5 over the share at 0 against
  Moe, Kratter and Badenes 2019's 53 +-12, 40 +-6, 24 +-4 and 10 +-3 percent over their line at 0 (0.20), widened by
  three standard errors of the sampled ratio. The level at 0 is printed against Table 13's 0.15 +-0.03 and never fails
  (the ratios are compared, not the absolutes: Boss, 2026-10-21, 7:26 AM EDT, "Ratios, report levels").
- **7 solar masses (mid B)** at [Fe/H] -0.2 and 0: at least 23 - 7 = 16 percent (Moe and Di Stefano 2017).
- **28 solar masses (O)** at [Fe/H] -0.2 and 0: between 35 - 12 = 23 and 69 + 9 = 78 percent (Moe and Di Stefano 2017;
  Sana et al. 2012), widened by three standard errors.
- **O and early B, flatness**: the share at [Fe/H] -0.7 over the share at +0.1 inside `FLATNESS_EDGE`, widened by three
  standard errors. The edge is a stand-in (the sources' own range once they are read, Boss, 2026-10-17, 9:54 PM EDT: the
  edge widens to it, not narrowed to O only); the first run gave 1.32 for O and 1.43 for B.
- **The outer share (10 AU and out)**: the share with a wide companion beyond 10 AU does not move with [Fe/H] for the
  M, G and A bins (the O and B bins are skipped where the companion chance saturates): two [Fe/H] bins differ by no more
  than three standard errors of their difference.

The bins of the other masses are printed and never fail.
"""

import math
from dataclasses import dataclass
from types import SimpleNamespace

from planetgen.generation import companion_classes, wide_binary
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

SIGMA = 3.0
"""float: How many standard errors a bin may stray from its range."""

DRAWS = 4000
"""int: Primaries drawn in each bin."""

FE_H_BINS = (-3.0, -1.0, -0.7, -0.2, 0.0, 0.1, 0.5)
"""tuple: The [Fe/H] of each bin, dex."""

MASS_BINS = (
    ("M dwarf", 0.4, "M"), ("G dwarf", 1.0, "G"), ("A star", 3.0, "A"), ("mid B", 7.0, "B"), ("early B", 12.0, "B"),
    ("O star", 28.0, "O"),
)
"""tuple: `(name, born mass in solar masses, spectral letter)` of each primary-mass bin."""

G_LINE = ((-3.0, 0.53, 0.12), (-1.0, 0.40, 0.06), (-0.2, 0.24, 0.04), (0.5, 0.10, 0.03))
"""tuple: `(fe_h, share, error)` of Moe, Kratter and Badenes 2019 for solar-type primaries."""

G_LINE_AT_ZERO = 0.20
"""float: The pooled line at [Fe/H] 0 (between its -0.2 and +0.5 points), the denominator of the ratios."""

G_LEVEL = (0.15, 0.03)
"""tuple: Moe and Di Stefano 2017's Table 13, companions inside log P 3.7 per solar-type primary, and its error."""

MID_B_FLOOR = 0.23 - 0.07
"""float: Moe and Di Stefano 2017, mid-B primaries (7 +-2 solar masses), at least 23 +-7 percent."""

O_RANGE = (0.35 - 0.12, 0.69 + 0.09)
"""tuple: O primaries: Moe and Di Stefano 2017's 35 +-12 percent (periods 100 to 1500 d) to Sana et al. 2012's 69 +-9."""

FLATNESS_EDGE = 1.5
"""float: How far the share at [Fe/H] -0.7 may sit from the share at +0.1 for O and early B stars, as a ratio either way.
A stand-in for the range the sources give (the draw gave 1.32 and 1.43)."""

OUTER_FLAT = ("M dwarf", "G dwarf", "A star")
"""tuple: The bins whose outer share must not follow [Fe/H]."""

FLATNESS = ("early B", "O star")
"""tuple: The bins whose share must stay within `FLATNESS_EDGE` across [Fe/H] -0.7 to +0.1."""


@dataclass(frozen=True)
class Bin:
    """The primaries drawn at one mass and [Fe/H]: how many, and how many had a companion inside 10 AU, close, or beyond."""
    name: str
    mass_solar: float
    fe_h: float
    draws: int
    inside: int
    close: int
    outer: int

    @property
    def share(self):
        """The share with a companion inside 10 AU."""
        return self.inside / self.draws

    @property
    def outer_share(self):
        """The share with a wide companion beyond 10 AU."""
        return self.outer / self.draws

    def error(self, successes=None):
        """The binomial standard error of the share (of `successes`, the inside count when omitted)."""
        successes = self.inside if successes is None else successes
        p = successes / self.draws
        return math.sqrt(max(p * (1.0 - p), 1.0 / self.draws) / self.draws)


@dataclass(frozen=True)
class Flag:
    """A bin or ratio that left its range."""
    row: str
    where: str
    value: float
    low: float
    high: float

    def label(self):
        return f"{self.row} at {self.where}: {self.value:.3f} outside {self.low:.3f} to {self.high:.3f}"


def _primary(letter, mass, fe_h):
    system = object.__new__(StarSystem)
    config = SystemConfig()
    config.PLACE_TERMS = None
    system.system_config = config
    system.primary_star = SimpleNamespace(type=letter + "2V", metallicity_fe_h=fe_h, age=5.0, initial_mass_sol=mass,
                                          mass=mass * 1.98841e30)
    return system


def draw_bin(name, mass, letter, fe_h, draws=DRAWS):
    """The bin of `draws` primaries of `mass` and `fe_h`, each decided by the generator's own companion draw."""
    system = _primary(letter, mass, fe_h)
    inside = close = outer = 0
    with draw.bound(f"close-binary-check:{name}:{fe_h}"):
        for _ in range(draws):
            kind, ratio, _reason = system._companion()
            if kind == "close":
                inside += 1
                close += 1
            elif kind == "wide":
                secondary = SimpleNamespace(initial_mass_sol=mass * ratio)
                limit = wide_binary.wide_pair_limit_au(system.system_config, system.primary_star, secondary)
                axis = wide_binary.sample_wide_binary_separation_au(limit, mass, mass * (1.0 + ratio),
                                                                    system._wide_inner_boost)
                if axis < companion_classes.INNER_LIMIT_AU:
                    inside += 1
                else:
                    outer += 1
    return Bin(name, mass, fe_h, draws, inside, close, outer)


def run_check(draws=DRAWS):
    """Every bin: `{(name, fe_h): Bin}`."""
    return {(name, fe_h): draw_bin(name, mass, letter, fe_h, draws)
            for name, mass, letter in MASS_BINS for fe_h in FE_H_BINS}


def _ratio(numerator, denominator):
    """`(ratio, standard error)` of two bins' shares (the delta method)."""
    a, b = numerator.share, denominator.share
    ratio = a / b
    error = ratio * math.sqrt((numerator.error() / a) ** 2 + (denominator.error() / b) ** 2)
    return ratio, error


def flags(bins, sigma=SIGMA):
    """The `Flag`s of every row of the check, empty when every bin is inside its range."""
    out = []
    zero = bins[("G dwarf", 0.0)]
    for fe_h, value, error in G_LINE:
        ratio, ratio_error = _ratio(bins[("G dwarf", fe_h)], zero)
        low, high = (value - error) / G_LINE_AT_ZERO - sigma * ratio_error, (value + error) / G_LINE_AT_ZERO + sigma * ratio_error
        if not low <= ratio <= high:
            out.append(Flag("G dwarf ratio to [Fe/H] 0", f"[Fe/H] {fe_h:+.1f}", ratio, low, high))
    for fe_h in (-0.2, 0.0):
        found = bins[("mid B", fe_h)]
        if found.share < MID_B_FLOOR - sigma * found.error():
            out.append(Flag("mid B share inside 10 AU", f"[Fe/H] {fe_h:+.1f}", found.share, MID_B_FLOOR, 1.0))
        found = bins[("O star", fe_h)]
        low, high = O_RANGE[0] - sigma * found.error(), O_RANGE[1] + sigma * found.error()
        if not low <= found.share <= high:
            out.append(Flag("O share inside 10 AU", f"[Fe/H] {fe_h:+.1f}", found.share, low, high))
    for name in FLATNESS:
        ratio, error = _ratio(bins[(name, -0.7)], bins[(name, 0.1)])
        low, high = 1.0 / FLATNESS_EDGE - sigma * error, FLATNESS_EDGE + sigma * error
        if not low <= ratio <= high:
            out.append(Flag(f"{name} flatness, [Fe/H] -0.7 over +0.1", "-0.7 / +0.1", ratio, low, high))
    for name in OUTER_FLAT:
        shares = [bins[(name, fe_h)] for fe_h in FE_H_BINS]
        for one in shares:
            for other in shares:
                if one.fe_h < other.fe_h:
                    difference = abs(one.outer_share - other.outer_share)
                    allowed = sigma * math.hypot(one.error(one.outer), other.error(other.outer))
                    if difference > allowed:
                        out.append(Flag(f"{name} outer share", f"[Fe/H] {one.fe_h:+.1f} vs {other.fe_h:+.1f}",
                                        difference, 0.0, allowed))
    return out


def format_report(bins, found=None):
    """The audit report: every bin's share inside 10 AU (and the outer share), then the flags."""
    found = flags(bins) if found is None else found
    lines = ["Close-binary check (TEST.196): share of primaries with a companion inside 10 AU, outer share in brackets",
             "", "  [Fe/H]  " + "  ".join(f"{name:>14}" for name, _m, _l in MASS_BINS)]
    for fe_h in FE_H_BINS:
        cells = [f"{bins[(name, fe_h)].share:6.3f} ({bins[(name, fe_h)].outer_share:5.3f})" for name, _m, _l in MASS_BINS]
        lines.append(f"  {fe_h:+5.1f}   " + "  ".join(f"{cell:>14}" for cell in cells))
    zero = bins[("G dwarf", 0.0)]
    lines += ["", f"  G dwarf level at [Fe/H] 0: {zero.share:.3f} against Table 13's {G_LEVEL[0]:.2f} +-{G_LEVEL[1]:.2f} "
                  "(reported, never fails)", ""]
    lines += ([f"  OUT OF RANGE: {flag.label()}" for flag in found] if found else ["  Every bin with a range is inside it."])
    return "\n".join(lines) + "\n"
