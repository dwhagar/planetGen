# planetgen/util/draw.py

"""
Draw
====

Every random draw generation makes goes through here (GEN.56; see
docs/design/generation/reproducible-galaxies.md, section 5).

Python promises the same numbers from one seed across versions only for
`random.Random.random()` (and the seeding itself), not for `choice`,
`uniform`, `gauss`, `shuffle`, `sample` or `randrange`, whose recipes have
changed between releases. So a `Stream` keeps a `random.Random` for its
`random()` alone and builds every other draw on it here, where the recipe
is ours and stays put.

Draws come from the stream bound to the running unit of work: a sector
fill, a run of scatter shells, a backfill block (`planetgen.galaxy.seed.seeded`
binds one from the unit's seed with `bound`). The module-level functions
(`draw.uniform(...)`, `draw.choice(...)`) read that stream, so generators
need no extra argument, and the binding is per thread (a `ContextVar`),
so the web site's threads never share a unit's stream. Outside any unit
draws come from the process's run stream, which `planetgen` seeds from its
run seed (`set_run_seed`) and which otherwise starts from the operating
system's random source.

The module itself has the same draw functions as a `Stream`, so a function
taking an `rng` can default to `rng=draw`.
"""

import bisect
import contextlib
import contextvars
import math
import random as _stdlib_random
import secrets

import numpy as np
from numpy.random import Generator, MT19937

_TWO_PI = 2.0 * math.pi

BULK_MIN = 1000
"""int: The fewest draws `Stream.random_array` and `Stream.randbits_list` hand to NumPy: moving the generator's state
there and back costs about as much as a thousand `random()` calls (PERF.91)."""


class Stream:
    """
    One seeded stream of draws. `seed` is anything `random.Random` takes
    (an int, a str or bytes); a str or bytes seed is hashed with SHA-512
    by `random.Random`, so it doesn't depend on `PYTHONHASHSEED`.
    """

    __slots__ = ("_random",)

    def __init__(self, seed):
        self._random = _stdlib_random.Random(seed)

    def random(self):
        """A float in [0, 1)."""
        return self._random.random()

    def uniform(self, a, b):
        """A float between `a` and `b` (either order)."""
        return a + (b - a) * self._random.random()

    def randrange(self, start, stop=None):
        """An int in [start, stop), or [0, start) with one argument."""
        if stop is None:
            start, stop = 0, start
        span = int(stop) - int(start)
        if span <= 0:
            raise ValueError(f"empty range for randrange({start}, {stop})")
        return int(start) + min(int(self._random.random() * span), span - 1)

    def randint(self, a, b):
        """An int in [a, b], both ends included."""
        return self.randrange(a, int(b) + 1)

    def choice(self, seq):
        """One element of the non-empty sequence `seq`."""
        if not len(seq):
            raise IndexError("cannot choose from an empty sequence")
        return seq[self.randrange(len(seq))]

    def choices(self, population, weights=None, k=1):
        """`k` elements of `population` drawn with replacement, by
        `weights` when given (relative, not cumulative)."""
        n = len(population)
        if not n:
            raise IndexError("cannot choose from an empty population")
        if weights is None:
            return [population[self.randrange(n)] for _ in range(k)]
        if len(weights) != n:
            raise ValueError("the number of weights does not match the population")
        cumulative = []
        total = 0.0
        for weight in weights:
            total += weight
            cumulative.append(total)
        if not total > 0.0:
            raise ValueError("total of weights must be greater than zero")
        last = n - 1
        return [population[min(bisect.bisect_right(cumulative, self._random.random() * total, 0, last), last)]
                for _ in range(k)]

    def shuffle(self, items):
        """Shuffles the list `items` in place (Fisher-Yates)."""
        for i in range(len(items) - 1, 0, -1):
            j = self.randrange(i + 1)
            items[i], items[j] = items[j], items[i]

    def sample(self, population, k):
        """`k` distinct elements of `population`, in draw order."""
        pool = list(population)
        if not 0 <= k <= len(pool):
            raise ValueError("sample larger than population or is negative")
        for i in range(k):
            j = i + self.randrange(len(pool) - i)
            pool[i], pool[j] = pool[j], pool[i]
        return pool[:k]

    def gauss(self, mu=0.0, sigma=1.0):
        """A normal draw (Box-Muller, two `random()` calls each time)."""
        u1 = 1.0 - self._random.random()  # (0, 1], so the log is finite
        u2 = self._random.random()
        return mu + sigma * math.sqrt(-2.0 * math.log(u1)) * math.cos(_TWO_PI * u2)

    normal = gauss

    def gamma(self, shape):
        """A gamma draw of the given `shape` and scale 1 (a shape under 1 is raised by one and
        scaled by a uniform to the power 1 / shape). Built on `random()` and `gauss` here, so it is ours and stays put."""
        if not shape > 0.0:
            raise ValueError("the shape of a gamma draw must be positive")
        if shape < 1.0:
            u = 1.0 - self._random.random()     # (0, 1]
            return self.gamma(shape + 1.0) * u ** (1.0 / shape)
        d = shape - 1.0 / 3.0
        c = 1.0 / math.sqrt(9.0 * d)
        while True:
            x = self.gauss()
            v = 1.0 + c * x
            if v <= 0.0:
                continue
            v = v * v * v
            u = 1.0 - self._random.random()     # (0, 1], so the log is finite
            if u < 1.0 - 0.0331 * x ** 4 or math.log(u) < 0.5 * x * x + d * (1.0 - v + math.log(v)):
                return d * v

    def dirichlet(self, alphas):
        """A list of proportions that adds up to 1 within 1e-15, Dirichlet distributed with the concentrations `alphas`
        (one gamma draw each, divided by their sum). Where every draw underflows to zero (all alphas tiny) the whole
        goes to the component with the largest alpha."""
        if not len(alphas):
            raise ValueError("a Dirichlet draw needs at least one concentration")
        gammas = [self.gamma(a) for a in alphas]
        total = math.fsum(gammas)
        if not total > 0.0:
            best = max(range(len(alphas)), key=lambda i: alphas[i])
            return [1.0 if i == best else 0.0 for i in range(len(alphas))]
        return [g / total for g in gammas]

    def getrandbits(self, k):
        """A `k`-bit non-negative int, 32 bits per `random()` call."""
        value = 0
        bits = 0
        while bits < k:
            value = (value << 32) | int(self._random.random() * 4294967296.0)
            bits += 32
        return value >> (bits - k)

    def random_array(self, count):
        """
        `count` floats in [0, 1) as an array: exactly the numbers `count` calls of `random()` give, and the stream
        is left in exactly the state they leave it in (PERF.91). From `BULK_MIN` draws on the generator's state goes
        to NumPy's `MT19937`, which is the same generator with the same 53-bit float, a batch is drawn in C and the
        state comes back.
        """
        if count < BULK_MIN:
            random = self._random.random
            return np.fromiter((random() for _ in range(count)), dtype=float, count=count)
        version, internal, gauss_next = self._random.getstate()
        bit_generator = MT19937()
        bit_generator.state = {"bit_generator": "MT19937",
                               "state": {"key": np.array(internal[:624], dtype=np.uint32), "pos": internal[624]}}
        values = Generator(bit_generator).random(count)
        after = bit_generator.state["state"]
        self._random.setstate((version, tuple(int(word) for word in after["key"]) + (int(after["pos"]),), gauss_next))
        return values

    def randbits_list(self, k, count):
        """
        `count` values of `getrandbits(k)` as a list, from the same draws and leaving the same state as `count`
        calls of it (PERF.91): in bulk from `BULK_MIN` draws on, for `k` up to 64.
        """
        words = -(-k // 32)
        if count * words < BULK_MIN or words > 2:
            return [self.getrandbits(k) for _ in range(count)]
        parts = (self.random_array(count * words).reshape(count, words) * 4294967296.0).astype(np.uint64)
        values = parts[:, 0]
        for column in range(1, words):
            values = (values << np.uint64(32)) | parts[:, column]
        return (values >> np.uint64(words * 32 - k)).tolist()

    def getstate(self):
        """The stream's state, for `setstate` (a save retried after a
        deadlock replays the same draws)."""
        return self._random.getstate()

    def setstate(self, state):
        self._random.setstate(state)


_run_stream = Stream(secrets.randbits(128))
"""Stream: Draws made outside any unit. `set_run_seed` replaces it."""

_bound = contextvars.ContextVar("planetgen_draw_stream", default=None)


def set_run_seed(seed):
    """Seeds the run stream (`planetgen`'s run seed, a work queue's)."""
    global _run_stream
    _run_stream = Stream(seed)


def current():
    """The stream draws come from: the unit's, else the run's."""
    stream = _bound.get()
    return _run_stream if stream is None else stream


@contextlib.contextmanager
def bound(stream):
    """Draws come from `stream` (a `Stream`, or a seed for one) for the
    length of the block, in this thread only."""
    if not isinstance(stream, Stream):
        stream = Stream(stream)
    token = _bound.set(stream)
    try:
        yield stream
    finally:
        _bound.reset(token)


def random():
    return current().random()


def uniform(a, b):
    return current().uniform(a, b)


def randrange(start, stop=None):
    return current().randrange(start, stop)


def randint(a, b):
    return current().randint(a, b)


def choice(seq):
    return current().choice(seq)


def choices(population, weights=None, k=1):
    return current().choices(population, weights=weights, k=k)


def shuffle(items):
    current().shuffle(items)


def sample(population, k):
    return current().sample(population, k)


def gauss(mu=0.0, sigma=1.0):
    return current().gauss(mu, sigma)


normal = gauss


def gamma(shape):
    return current().gamma(shape)


def dirichlet(alphas):
    return current().dirichlet(alphas)


def getrandbits(k):
    return current().getrandbits(k)


def getstate():
    return current().getstate()


def setstate(state):
    current().setstate(state)


_tally = contextvars.ContextVar("planetgen_redraw_tally", default=None)


def redrawn(kind, count=1):
    """Notes that a drawn item was thrown away and drawn again (or moved, or cut) to make it fit: `count` redraws of `kind`
    (a name of `planetgen.generation.redraw_kinds.KINDS`). It reads no stream and draws nothing, so what is generated is the
    same with or without a tally; outside `counting` it does nothing (TEST.476)."""
    tally = _tally.get()
    if tally is not None:
        tally[kind] = tally.get(kind, 0) + count


@contextlib.contextmanager
def counting():
    """Counts every `redrawn` call made in this thread for the length of the block, as a dict `{kind: count}` that fills as the
    block runs. Blocks nest: the inner one counts its own and the outer one sees them too."""
    tally = {}
    token = _tally.set(tally)
    try:
        yield tally
    finally:
        _tally.reset(token)
        if token.old_value is not contextvars.Token.MISSING and token.old_value is not None:
            for kind, count in tally.items():
                token.old_value[kind] = token.old_value.get(kind, 0) + count
