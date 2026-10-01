# tests/fuzz_support.py

"""
Shared settings for the `test_fuzz_*.py` files: property-based,
brute-force tests (via `hypothesis`) that throw generated inputs at a
function until something breaks, then shrink the failure down to the
smallest input that still breaks it.

They sit beside the seeded `test_bughunt_*.py` files rather than
replacing them: those replay a fixed seed list through whole generation
runs, while these attack one function at a time with every edge the
strategy can reach (zero, negative, NaN, infinity, huge, boundary
values, empty and hostile strings).

Two profiles, picked with `PLANETGEN_FUZZ_PROFILE`:

- `ci` (default): a bounded, derandomized run -- the same examples every
  time, so a normal `pytest` stays fast and a CI failure always
  reproduces.
- `deep`: many more examples from a fresh random seed each run, for
  hunting (`PLANETGEN_FUZZ_PROFILE=deep pytest src/tests/test_fuzz_*.py`,
  and the weekly `Deep fuzz` workflow). A failure prints the
  `@reproduce_failure` line to paste into the test.

`PLANETGEN_FUZZ_EXAMPLES` overrides the example count of either profile.
"""

import contextlib
import math
import os
import random
import secrets
from unittest import mock

from hypothesis import HealthCheck, settings
from hypothesis import strategies as st

_SUPPRESSED = [HealthCheck.too_slow, HealthCheck.data_too_large, HealthCheck.filter_too_much]

settings.register_profile(
    "ci",
    max_examples=int(os.environ.get("PLANETGEN_FUZZ_EXAMPLES") or 60),
    deadline=None,
    derandomize=True,
    suppress_health_check=_SUPPRESSED,
    print_blob=True,
)
settings.register_profile(
    "deep",
    max_examples=int(os.environ.get("PLANETGEN_FUZZ_EXAMPLES") or 3000),
    deadline=None,
    suppress_health_check=_SUPPRESSED,
    print_blob=True,
)
settings.load_profile(os.environ.get("PLANETGEN_FUZZ_PROFILE") or "ci")


@contextlib.contextmanager
def deterministic_entropy(seed):
    """Seeds the global `random` module AND routes every `secrets` call the
    generators make through one seeded `random.Random`, so a whole system
    (planet classes, moon coin-flips, `reseed_rng()` reseeds, flavor text)
    is a pure function of `seed` for the duration of the block."""
    rng = random.Random(seed)
    with mock.patch.object(secrets, "randbits", rng.getrandbits), \
            mock.patch.object(secrets, "randbelow", rng.randrange), \
            mock.patch.object(secrets, "choice", rng.choice):
        random.seed(seed)
        yield rng


def scaled(n):
    """`n` examples under `ci`, scaled up with the active profile -- for a
    test whose single example is expensive (a whole system or sector), so
    it can ask for fewer than the profile default without losing the
    deep run's extra depth."""
    base = settings.get_profile("ci").max_examples
    return max(1, round(n * settings().max_examples / base))


# Finite floats with no NaN/infinity, in a range physics code can take.
finite = st.floats(allow_nan=False, allow_infinity=False, min_value=-1e12, max_value=1e12)

# Every float hypothesis can make, including NaN, +/-inf, -0.0 and
# subnormals -- for "garbage in must not crash uncontrolled" tests.
any_float = st.floats(allow_nan=True, allow_infinity=True)

# Non-finite floats only.
non_finite = st.sampled_from([math.nan, math.inf, -math.inf])

# Hostile text: control characters, surrogates excluded (not encodable),
# markup, SQL and path fragments mixed in with arbitrary unicode.
_HOSTILE_FRAGMENTS = [
    "", " ", "\x00", "\n", "\r\n", "\t", "'", '"', "`", "\\", "%", "_", "%00",
    "<script>alert(1)</script>", "</textarea>", "{{7*7}}", "{% raw %}",
    "' OR '1'='1", "1; DROP TABLE sectors;--", "../../etc/passwd", "..\\..\\boot.ini",
    "javascript:alert(1)", "data:text/html,x", "‮", "﻿", "ｘ", "𝔘𝔫𝔦𝔠𝔬𝔡𝔢",
    "-1", "0", "1e309", "NaN", "inf", "-inf", "9" * 40, "0x10", "1_000",
]
hostile_text = st.one_of(
    st.sampled_from(_HOSTILE_FRAGMENTS),
    st.text(alphabet=st.characters(blacklist_categories=("Cs",)), max_size=80),
    st.lists(st.sampled_from(_HOSTILE_FRAGMENTS), min_size=1, max_size=4).map("".join),
)
