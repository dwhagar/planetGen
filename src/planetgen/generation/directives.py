"""
Generation directives for a sector (GEN.96): minimums a run asks of a
sector, such as "at least 5 systems", "at least 2 habitable worlds" or
"at least 3 G-type stars".

A sector is drawn as usual; when it misses a directive it is thrown away
and drawn again from a fresh stream, up to `DEFAULT_ATTEMPTS` (K) times, and
only the winner is saved (docs/design/generation/sampling-backfill-and-resume.md,
section 8). A habitable-worlds minimum still missed after K draws takes the
fill (GEN.180): the sector is built once more with that many systems each
asking for a habitable world (the odds raised, nothing forced, GEN.1130), and is marked `met_forced`. A request that
cannot occur at all is refused before any drawing (`impossible`). Attempt 0 uses the sector's plain seed, so a sector that
meets its directives first time is the sector a run with none would have
made. Attempt `a` >= 1 uses a seed from the galaxy seed, the sector's
address, the directives' digest and `a`, so the same galaxy, directives
and sector always give the same sector.
"""
import hashlib
import json
import math
import re
from collections import Counter

from planetgen import tuning
from planetgen.galaxy import seed as galaxySeed
from planetgen.util import draw

DEFAULT_ATTEMPTS = tuning.DIRECTIVE_ATTEMPTS

MET_NATURALLY = "met_naturally"
MET_AFTER_ATTEMPTS = "met_after_k_attempts"
MET_FORCED = "met_forced"
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
        counts[HABITABLE] += getattr(system, "habitable_class_count", 0) or 0
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
        if self.status == MET_FORCED:
            return (f"directives not met in {self.attempts - 1} draws, so a fill asked for a habitable world in each of "
                    f"the systems it needed, which raises the odds and forces none (met_forced)")
        misses = ", ".join(f"{key} {got}/{want}" for key, (want, got) in sorted(self.missed.items()))
        return f"directives not met after {self.attempts} draws ({misses}); the closest draw was kept"


def _score(missed):
    """Total shortfall; the least wins when no draw meets everything."""
    return sum(want - got for want, got in missed.values())


HABITABLE_COUNT_PMF = (0.7551, 0.1896, 0.0380, 0.0120, 0.0040, 0.0010)
"""tuple[float]: The share of star systems holding 0, 1, 2, 3, 4 and 5 habitable worlds, from 6,000 systems of
`StarSystem(SystemConfig())` (sampling-backfill-and-resume.md section 8.1), normalised. Tier: verified by simulation of this
generator, which is the thing the directive redraws, so a change to the habitable draw must remeasure it."""


def _poisson_tail(mean, need):
    """P(X >= `need`) for X Poisson(`mean`), summed from the top so a tiny tail keeps its digits."""
    if need <= 0:
        return 1.0
    if mean <= 0:
        return 0.0
    log_term = -mean + need * math.log(mean) - math.lgamma(need + 1)
    if log_term < -700.0 and need > mean:       # the first term underflows: the whole tail is below 1e-300
        return 0.0
    term, total, k = math.exp(log_term), 0.0, need
    while True:
        total += term
        k += 1
        term *= mean / k
        if term < total * 1e-17 and k > mean:
            return min(total, 1.0)


def _compound_poisson_tail(mean, need, pmf=HABITABLE_COUNT_PMF):
    """P(S >= `need`) where S sums the counts (`pmf`: index = count) of a Poisson(`mean`) number of systems, by Panjer's
    recursion; exact."""
    if need <= 0:
        return 1.0
    rates = [mean * share for share in pmf]
    mass = [math.exp(-mean * (1.0 - pmf[0]))]
    for s in range(1, need):
        mass.append(sum(k * rates[k] * mass[s - k] for k in range(1, len(pmf)) if s - k >= 0) / s)
    return max(0.0, 1.0 - sum(mass))


def _fixed_tail(systems, need, pmf=HABITABLE_COUNT_PMF):
    """P(S >= `need`) for exactly `systems` systems, by convolution capped at `need`; 1 when `systems` is too many to
    convolve cheaply (never a false refusal)."""
    if need <= 0:
        return 1.0
    if need > systems * (len(pmf) - 1):
        return 0.0
    if systems > 2000:
        return 1.0
    dist = [1.0]
    for _ in range(systems):
        step = [0.0] * min(len(dist) + len(pmf) - 1, need + 1)
        for have, mass in enumerate(dist):
            for k, share in enumerate(pmf):
                step[min(have + k, need)] += mass * share
        dist = step
    return dist[need] if len(dist) > need else 0.0


def impossible(directive, systems=None, mean=None):
    """
    Why `directive` cannot be met in a sector of `systems` systems (a count the run fixes) or of about `mean` (the mean
    of the Poisson draw of its count), or `None` when a sector can reach every minimum. Judged before any drawing from the
    compound-Poisson law of the habitable worlds (`HABITABLE_COUNT_PMF`) and the Poisson law of the system count: a minimum
    is impossible when a sector reaches it with a probability under `tuning.DIRECTIVE_IMPOSSIBLE_BELOW` and the forced
    fill (which needs that many systems) cannot reach it either. With neither `systems` nor `mean` nothing is judged.
    """
    if systems is None and mean is None:
        return None
    floor = tuning.DIRECTIVE_IMPOSSIBLE_BELOW
    for key, want in sorted(directive.minimums.items()):
        if want <= 0:
            continue
        holds = (1.0 if systems >= want else 0.0) if systems is not None else _poisson_tail(mean, want)
        if holds >= floor:
            continue
        if key == HABITABLE:
            natural = _fixed_tail(systems, want) if systems is not None else _compound_poisson_tail(mean, want)
            if natural >= floor:
                continue
        label = {SYSTEMS: "systems", HABITABLE: "habitable worlds"}.get(key, key.removeprefix(TYPE_PREFIX) + " stars")
        basis = (f"a sector of {systems} system{'s' if systems != 1 else ''}" if systems is not None
                 else f"a sector with about {mean:.3g} systems")
        reason = f"{directive_text(key, want)} cannot be met in {basis}"
        if key == HABITABLE:
            reason += (f": it would need at least {want} habitable worlds, and even forcing one into every system takes "
                       f"{want} systems")
        else:
            reason += f": it needs {want} {label} and the sector cannot hold that many"
        return reason + " (raise the number of systems or the density, or lower the minimum)."
    return None


def directive_text(key, want):
    """`"habitable >= 5"`."""
    return f"{key} >= {want}"


def generate(directive, build, galaxy_seed=None, address=None, max_attempts=DEFAULT_ATTEMPTS, build_forced=None):
    """
    Draws a sector with `build()` until it meets `directive`.

    Args:
        directive (Directive): The minimums wanted; empty means one draw.
        build (callable): Builds one sector from the ambient random stream.
        galaxy_seed (bytes, optional): With it, attempt `a` >= 1 binds its
            own seeded stream; without it the ambient stream just goes on.
        address: The sector's address, for the attempt seeds.
        max_attempts (int): Draws to make at most.
        build_forced (callable, optional): `build_forced(count)` builds a sector in which `count` systems are forced to
            hold a habitable world (GEN.180). Tried once, after `max_attempts` draws, when a habitable-worlds minimum is
            still missed; the sector it builds is kept when it meets every minimum (`MET_FORCED`) or comes nearer.

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
    attempts = max(1, max_attempts)
    if build_forced is not None and HABITABLE in best_missed:
        where = f"{galaxySeed.address_text(address)}/{directive.digest()}/forced" if galaxy_seed is not None else None
        count = directive.minimums[HABITABLE]
        if where is None:
            forced = build_forced(count)
        else:
            with draw.bound(galaxySeed.unit_seed(galaxy_seed, "sector-directive-forced", where)):
                forced = build_forced(count)
        missed = shortfall(directive, forced)
        if not missed:
            return forced, Outcome(MET_FORCED, attempts + 1, {})
        if _score(missed) < _score(best_missed):
            best, best_missed = forced, missed
    return best, Outcome(UNMET, attempts, best_missed)
