"""
Generation directives for a sector (GEN.96): minimums a run asks of a
sector, such as "at least 5 systems", "at least 2 habitable worlds" or
"at least 3 G-type stars".

A sector is drawn as usual; when it misses a directive it is thrown away
and drawn again from a fresh stream, up to `DEFAULT_ATTEMPTS` times, and
only the winner is saved (docs/design/sampling-backfill-and-resume.md,
section 8). Attempt 0 uses the sector's plain seed, so a sector that
meets its directives first time is the sector a run with none would have
made. Attempt `a` >= 1 uses a seed from the galaxy seed, the sector's
address, the directives' digest and `a`, so the same galaxy, directives
and sector always give the same sector.
"""
import hashlib
import json
import re
from collections import Counter

from planetgen.galaxy import seed as galaxySeed
from planetgen.util import draw

DEFAULT_ATTEMPTS = 200

MET_NATURALLY = "met_naturally"
MET_AFTER_ATTEMPTS = "met_after_k_attempts"
UNMET = "unmet"

SYSTEMS = "systems"
HABITABLE = "habitable"
TYPE_PREFIX = "type:"
STAR_LETTERS = "OBAFGKM"
SPECIAL_TYPES = {"wd": "white dwarf", "bh": "black hole", "ns": "neutron star"}

_PATTERN = re.compile(r"^\s*([A-Za-z:]+)\s*>=\s*(\d+)\s*$")


class DirectiveError(ValueError):
    """A directive that can't be read, or can never be met."""


class Directive:
    """A set of minimums: `minimums` maps a key (`systems`, `habitable`,
    `type:G`, `type:wd` ...) to the least count wanted."""

    def __init__(self, minimums=None):
        self.minimums = dict(minimums or {})

    def __bool__(self):
        return any(count > 0 for count in self.minimums.values())

    def canonical(self):
        """The directives as JSON with sorted keys: what the digest covers."""
        return json.dumps(self.minimums, sort_keys=True, separators=(",", ":"))

    def digest(self):
        return hashlib.sha256(self.canonical().encode("utf-8")).hexdigest()[:16]

    def describe(self):
        return ", ".join(f"{key} >= {count}" for key, count in sorted(self.minimums.items()))


def parse(texts):
    """
    A `Directive` from `--directive` strings: `systems>=N`, `habitable>=N`
    or `type:X>=N` (X one of O B A F G K M wd bh ns). The same key given
    twice keeps the larger count.

    Raises:
        DirectiveError: Anything else.
    """
    minimums = {}
    for text in texts or ():
        match = _PATTERN.match(text)
        if not match:
            raise DirectiveError(f"{text!r}: expected systems>=N, habitable>=N or type:X>=N")
        key, count = match.group(1), int(match.group(2))
        if key.lower() in (SYSTEMS, HABITABLE):
            key = key.lower()
        elif key.lower().startswith(TYPE_PREFIX):
            kind = key[len(TYPE_PREFIX):]
            kind = kind.upper() if kind.upper() in STAR_LETTERS and len(kind) == 1 else kind.lower()
            if kind not in SPECIAL_TYPES and kind not in STAR_LETTERS:
                raise DirectiveError(f"{text!r}: star type must be one of {', '.join(STAR_LETTERS)}, "
                                     f"{', '.join(SPECIAL_TYPES)}")
            key = TYPE_PREFIX + kind
        else:
            raise DirectiveError(f"{text!r}: unknown directive {key!r}")
        minimums[key] = max(count, minimums.get(key, 0))
    return Directive(minimums)


def _star_key(star):
    yerkes = getattr(star, "yerkes_class", None)
    if yerkes in ("D", "VII"):
        return TYPE_PREFIX + "wd"
    if yerkes == "BH":
        return TYPE_PREFIX + "bh"
    if yerkes == "NS":
        return TYPE_PREFIX + "ns"
    return TYPE_PREFIX + ((getattr(star, "type", None) or "?")[0])


def measure(sector):
    """What `sector` holds, as counts keyed like a directive's minimums."""
    counts = Counter()
    for entry in sector.entries:
        system = entry.star_system
        counts[SYSTEMS] += 1
        counts[HABITABLE] += getattr(system, "hab_count", 0) or 0
        counts[_star_key(system.star)] += 1
    return counts


def shortfall(directive, sector):
    """`{key: (wanted, got)}` for every minimum `sector` misses."""
    got = measure(sector)
    return {key: (want, got[key]) for key, want in directive.minimums.items() if got[key] < want}


class Outcome:
    """How a directed sector came out: `status`, `attempts` made, and
    `missed` (`shortfall`'s result, empty when met)."""

    def __init__(self, status, attempts, missed):
        self.status = status
        self.attempts = attempts
        self.missed = missed

    def describe(self):
        if self.status == MET_NATURALLY:
            return "directives met on the first draw"
        if self.status == MET_AFTER_ATTEMPTS:
            return f"directives met after {self.attempts} draws"
        misses = ", ".join(f"{key} {got}/{want}" for key, (want, got) in sorted(self.missed.items()))
        return f"directives not met after {self.attempts} draws ({misses}); the closest draw was kept"


def _score(missed):
    """Total shortfall; the least wins when no draw meets everything."""
    return sum(want - got for want, got in missed.values())


def generate(directive, build, galaxy_seed=None, address=None, max_attempts=DEFAULT_ATTEMPTS):
    """
    Draws a sector with `build()` until it meets `directive`.

    Args:
        directive (Directive): The minimums wanted; empty means one draw.
        build (callable): Builds one sector from the ambient random stream.
        galaxy_seed (bytes, optional): With it, attempt `a` >= 1 binds its
            own seeded stream; without it the ambient stream just goes on.
        address: The sector's address, for the attempt seeds.
        max_attempts (int): Draws to make at most.

    Returns:
        tuple: `(sector, Outcome)`; when nothing met the directive, the
            draw nearest to it.
    """
    sector = build()
    if not directive:
        return sector, Outcome(MET_NATURALLY, 1, {})
    missed = shortfall(directive, sector)
    if not missed:
        return sector, Outcome(MET_NATURALLY, 1, {})
    best, best_missed = sector, missed
    for attempt in range(1, max(1, max_attempts)):
        if galaxy_seed is None:
            candidate = build()
        else:
            where = f"{galaxySeed.address_text(address)}/{directive.digest()}/{attempt}"
            with draw.bound(galaxySeed.unit_seed(galaxy_seed, "sector-directive", where)):
                candidate = build()
        missed = shortfall(directive, candidate)
        if not missed:
            return candidate, Outcome(MET_AFTER_ATTEMPTS, attempt + 1, {})
        if _score(missed) < _score(best_missed):
            best, best_missed = candidate, missed
    return best, Outcome(UNMET, max(1, max_attempts), best_missed)
