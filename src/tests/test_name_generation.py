"""
Regression tests for `stellarObjects.utils.generate_phoneme_salad_name`.

A base name containing a literal apostrophe (e.g. `PLANET_NAMES`'s
"Hi'iaka") splits into a syllable that *starts* with "'" (see
`split_into_syllables`). `UNIVERSAL_PHONEMES` separately contains chunks
that *end* with "'" (e.g. "ch'", from the Andean/glottalized set). When
`generate_phoneme_salad_name` shuffles syllables and splices in a universal
phoneme, these two independently-apostrophed fragments can land adjacent to
each other, producing a literal "''" in the assembled name. The final
per-word-capitalization step used to split on "'" and blindly index
`part[0]` on every resulting piece, which raised `IndexError: string index
out of range` on the empty piece between the two apostrophes.

This was rare enough (needs a specific base name, a specific phoneme, and a
specific shuffle outcome to all coincide) that a plain sweep of a few
thousand calls usually misses it, so this test pins down the exact
previously-crashing scenario deterministically via monkeypatched
random/secrets calls, in addition to a large-iteration sweep as a
belt-and-suspenders check.

Run with: pytest tests/test_name_generation.py
"""
import random

import pytest

from stellarObjects.names import (
    PLANET_NAMES, PLANET_PREFIXES, PLANET_SUFFIXES,
    STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES,
    UNIVERSAL_PHONEMES,
)
from stellarObjects.utils import generate_phoneme_salad_name, split_into_syllables


def test_apostrophe_base_name_splits_into_a_leading_apostrophe_syllable():
    """Pins down the precondition the crash scenario below relies on."""
    syllables = split_into_syllables("Hi'iaka")
    assert any(s.startswith("'") for s in syllables)


def test_universal_phonemes_contains_a_trailing_apostrophe_chunk():
    """Pins down the other half of the precondition: a phoneme ending in "'"."""
    assert any(p.endswith("'") for p in UNIVERSAL_PHONEMES)


def test_adjacent_apostrophes_from_base_name_and_phoneme_do_not_crash(monkeypatch):
    """
    Reproduces the exact previously-crashing case: base name "Hi'iaka"
    (-> syllables ['Hi', "'ia", 'ka']), shuffled to ['ka', 'Hi'], with the
    universal phoneme "ch'" spliced in right before the "'ia" syllable --
    assembling to ".../..."kaHi" + "ch'" + "'ia" == "...ch''ia...", the
    literal double-apostrophe that used to raise IndexError.
    """
    choice_calls = {"n": 0}

    def fake_secrets_choice(seq):
        choice_calls["n"] += 1
        seq = list(seq)
        if seq == PLANET_NAMES:
            return "Hi'iaka"
        if seq == UNIVERSAL_PHONEMES:
            return "ch'"
        if seq == PLANET_PREFIXES:
            return "Wald"
        if seq == PLANET_SUFFIXES:
            return "os"
        return seq[0]

    random_calls = {"n": 0}

    def fake_random_random():
        # First call gates the UNIVERSAL_PHONEME_CHANCE splice -- must be
        # below the 0.4 threshold so the phoneme is actually inserted.
        return 0.0

    def fake_shuffle(lst):
        # Reproduce the specific shuffle outcome: ['Hi', "'ia", 'ka'] -> ['ka', 'Hi'].
        # (shuffle runs before the phoneme is spliced in, so it only ever
        # sees the two syllables split from "Hi'iaka".)
        lst[:] = ["ka", "Hi"]

    def fake_randint(a, b):
        # Insert position for the phoneme: index 2 in ['ka', 'Hi'] places it
        # at the end -> ['ka', 'Hi', "ch'"], immediately before "'ia" is
        # appended by the join, landing the two apostrophes adjacent.
        return 2

    monkeypatch.setattr("stellarObjects.utils.secrets.choice", fake_secrets_choice)
    monkeypatch.setattr("stellarObjects.utils.random.random", fake_random_random)
    monkeypatch.setattr("stellarObjects.utils.random.shuffle", fake_shuffle)
    monkeypatch.setattr("stellarObjects.utils.random.randint", fake_randint)

    # Should not raise IndexError.
    name = generate_phoneme_salad_name(PLANET_NAMES, PLANET_PREFIXES, PLANET_SUFFIXES)
    assert isinstance(name, str) and name


@pytest.mark.parametrize("name_list,prefixes,suffixes", [
    (PLANET_NAMES, PLANET_PREFIXES, PLANET_SUFFIXES),
    (STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES),
])
def test_generate_phoneme_salad_name_sweep_does_not_crash(name_list, prefixes, suffixes):
    """
    Fallback sweep: the adjacency that crashes generation depends on
    unseedable `secrets` rolls (base name, phoneme, prefix, suffix) as well
    as seedable `random` rolls (shuffle order, splice position), so it can't
    be pinned to a single deterministic seed the way the test above does.
    Running many iterations makes it very likely to hit the same class of
    bug again if a future change reintroduces it.
    """
    for _ in range(20000):
        name = generate_phoneme_salad_name(name_list, prefixes, suffixes)
        assert isinstance(name, str) and name
