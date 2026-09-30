# tests/test_fuzz_names.py

"""
Property-based / brute-force tests for name generation and name-uniqueness
decoration: `stellarObjects/names.py` (the word lists themselves),
`stellarObjects/nameUniqueness.py` (`resolve_greek_roman_collision`,
`resolve_diminutive`, `resolve_companion`, `strip_decoration`) and the
name-generation helpers in `stellarObjects/utils.py` that consume them
(`split_into_syllables`, `is_name_valid`, `split_long_word`,
`generate_phoneme_salad_name`, `generate_sector_name`).

Goes past `test_name_uniqueness.py`/`test_name_generation.py`/
`test_bughunt_name_exhaustion.py`: huge collision counts, arbitrary base
names, long simulated insert runs checking uniqueness, the decoration
stacks `_db.py` actually builds, hostile text into every pure helper, and
a bounded sample of real generated names checked against the offensive-
word list. See `fuzz_support.py` for the `ci`/`deep` profiles.
"""

import itertools
import string
from unittest import mock

import pytest
from hypothesis import assume, example, given, settings
from hypothesis import strategies as st

from stellarObjects import names, utils
from stellarObjects import nameUniqueness as nu
from tests.fuzz_support import hostile_text, scaled

DECORATION_TOKENS = (
    set(names.GREEK_LETTERS)
    | set(names.DIMINUTIVE_PREFIXES)
    | set(names.COMPANION_SUFFIXES)
    | set(names.ROMAN_NUMERALS_BY_VALUE.values())
)

# A plain base-name word: letters (plus an embedded apostrophe, as real
# names have), never itself one of the decoration tokens.
plain_word = st.text(alphabet=string.ascii_letters + "'", min_size=1, max_size=12).filter(
    lambda w: w not in DECORATION_TOKENS
)
plain_base = st.lists(plain_word, min_size=1, max_size=3).map(" ".join)

NAME_SOURCES = {
    "star": (names.STAR_NAMES, names.STAR_PREFIXES, names.STAR_SUFFIXES),
    "planet": (names.PLANET_NAMES, names.PLANET_PREFIXES, names.PLANET_SUFFIXES),
    "moon": (names.MOON_NAMES, names.MOON_PREFIXES, names.MOON_SUFFIXES),
}


# ---------------------------------------------------------------------------
# names.py -- the word lists every generator and decorator depends on
# ---------------------------------------------------------------------------

def test_decoration_vocabularies_are_pairwise_disjoint_single_tokens():
    vocabularies = {
        "greek": names.GREEK_LETTERS,
        "diminutive": names.DIMINUTIVE_PREFIXES,
        "companion": names.COMPANION_SUFFIXES,
        "roman": list(names.ROMAN_NUMERALS_BY_VALUE.values()),
    }
    for label, words in vocabularies.items():
        assert words, label
        assert len(set(words)) == len(words), f"{label} has duplicates"
        for word in words:
            assert word and word == word.strip() and " " not in word, (label, word)
    for (a, wa), (b, wb) in itertools.combinations(vocabularies.items(), 2):
        assert not set(wa) & set(wb), (a, b, set(wa) & set(wb))


def test_greek_and_roman_tables_are_consistent():
    assert len(names.GREEK_LETTERS) == 24
    assert list(names.ROMAN_NUMERALS_BY_VALUE) == names.ROMAN_NUMERAL_VALUES
    assert names.ROMAN_NUMERAL_VALUES == sorted(set(names.ROMAN_NUMERAL_VALUES))
    assert all(v > 0 for v in names.ROMAN_NUMERAL_VALUES)
    assert nu.GREEK_ROMAN_CAPACITY == len(names.GREEK_LETTERS) + len(names.ROMAN_NUMERAL_VALUES)


def test_offensive_word_list_is_usable_by_the_lowercase_substring_filter():
    # `is_name_valid` lowercases the candidate and does `word in name`: an
    # empty entry would reject every name (infinite retry), an uppercase or
    # padded entry could never match anything.
    assert names.NSFW_WORDS
    for word in names.NSFW_WORDS:
        assert word, "empty offensive word would reject every name"
        assert word == word.strip() and word == word.lower(), repr(word)


def test_filter_word_lists_are_sane():
    assert names.DICTIONARY_WORDS
    assert isinstance(names.WORD_SIZE_MEAN, int) and 1 <= names.WORD_SIZE_MEAN <= 30
    assert set(names.VOWELS) == set("aeiou")
    for cluster in names.BAD_CONSONANTS:
        assert cluster and cluster == cluster.lower() and not set(cluster) & set(names.VOWELS), cluster


@pytest.mark.parametrize("source", ["star", "planet", "moon", "sector"])
def test_base_lists_have_no_blank_entries(source):
    lists = NAME_SOURCES.get(source) or (names.SECTOR_NAMES, names.SECTOR_PREFIXES, names.SECTOR_SUFFIXES)
    for words in lists:
        assert words
        for word in words:
            assert isinstance(word, str) and word.strip() == word and word, (source, word)


# ---------------------------------------------------------------------------
# resolve_greek_roman_collision
# ---------------------------------------------------------------------------

@given(base=plain_base, count=st.one_of(st.integers(0, 60), st.integers(0, 10**40)))
def test_greek_roman_result_shape_for_any_count(base, count):
    new_name, rename = nu.resolve_greek_roman_collision(base, count)
    if count >= nu.GREEK_ROMAN_CAPACITY - 1:
        assert (new_name, rename) == (None, None)
        return
    assert isinstance(new_name, str) and base in new_name
    assert nu.strip_decoration(new_name) == base
    if count == 0:
        assert new_name == base and rename is None
    if rename is not None:
        old, renamed = rename
        assert count in (1, len(names.GREEK_LETTERS))
        assert nu.strip_decoration(old) == base == nu.strip_decoration(renamed)
        assert old != renamed and renamed != new_name


def test_greek_roman_full_sequence_for_one_base_is_unique_and_ordered():
    base = "Voranthis"
    live = []
    for count in itertools.count():
        new_name, rename = nu.resolve_greek_roman_collision(base, count)
        if new_name is None:
            break
        if rename is not None:
            old, renamed = rename
            assert live.count(old) == 1, (count, old, live)
            live[live.index(old)] = renamed
        live.append(new_name)
        assert len(set(live)) == len(live), live
    assert len(live) == nu.GREEK_ROMAN_CAPACITY - 1
    assert base not in live  # the bare row was always renamed away


@settings(max_examples=scaled(25))
@given(inserts=st.lists(st.sampled_from(["Aa", "Bb", "Cc", "Dd Ee", "Ff'g"]), min_size=1, max_size=400))
def test_same_level_simulated_insert_run_stays_unique(inserts):
    """A long run of sector (or system) inserts drawing from a tiny base
    pool, replaying `_db.reserve_sector_name`'s own bookkeeping: rename the
    old row when told to, draw a fresh base once a base is exhausted."""
    rows = []                  # current names, in insertion order
    base_of = {}               # current name -> base
    counts = {}                # base -> occurrence_count
    fresh = (f"Fresh{i}" for i in itertools.count())
    for base in inserts:
        while True:
            new_name, rename = nu.resolve_greek_roman_collision(base, counts.get(base, 0))
            if new_name is not None:
                break
            base = next(fresh)
        if rename is not None:
            old, renamed = rename
            assert old in base_of, (old, rows)
            rows[rows.index(old)] = renamed
            base_of[renamed] = base_of.pop(old)
        assert new_name not in base_of, (new_name, rows)
        rows.append(new_name)
        base_of[new_name] = base
        counts[base] = counts.get(base, 0) + 1
    assert len(set(rows)) == len(rows)
    for name in rows:
        assert nu.strip_decoration(name) == base_of[name], name


# ---------------------------------------------------------------------------
# resolve_diminutive / resolve_companion
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("resolver, words", [
    (nu.resolve_diminutive, names.DIMINUTIVE_PREFIXES),
    (nu.resolve_companion, names.COMPANION_SUFFIXES),
])
def test_index_resolvers_walk_every_word_once_then_stop(resolver, words):
    seen, index = [], None
    for _ in range(len(words) + 5):
        word, index = resolver(index)
        if word is None:
            assert index is None
            break
        seen.append(word)
    assert seen == list(words)
    # Exhaustion is sticky: resolving from None again restarts, but from the
    # last valid index it stays exhausted.
    assert resolver(len(words) - 1) == (None, None)


@given(index=st.one_of(st.none(), st.integers(0, 10**40)))
def test_index_resolvers_for_any_non_negative_index(index):
    for resolver, words in ((nu.resolve_diminutive, names.DIMINUTIVE_PREFIXES),
                            (nu.resolve_companion, names.COMPANION_SUFFIXES)):
        word, nxt = resolver(index)
        expected_next = 0 if index is None else index + 1
        if expected_next >= len(words):
            assert (word, nxt) == (None, None)
        else:
            assert nxt == expected_next and word == words[nxt]


# ---------------------------------------------------------------------------
# strip_decoration
# ---------------------------------------------------------------------------

@given(name=hostile_text)
def test_strip_decoration_never_raises_and_only_removes_whole_edge_tokens(name):
    out = nu.strip_decoration(name)
    assert isinstance(out, str)
    tokens = name.split(" ")
    kept = out.split(" ") if out else []
    # The result is a contiguous run of the original tokens.
    if kept:
        n = len(kept)
        assert any(tokens[i:i + n] == kept for i in range(len(tokens) - n + 1))
    if kept:
        start = next(i for i in range(len(tokens) - len(kept) + 1) if tokens[i:i + len(kept)] == kept)
        removed = tokens[:start] + tokens[start + len(kept):]
    else:
        removed = tokens
    assert all(t in DECORATION_TOKENS for t in removed) or not kept and name.strip(" ") == "" or \
        all(t in DECORATION_TOKENS or t == "" for t in removed), (name, out)
    if not (set(tokens[:2]) | set(tokens[-2:])) & DECORATION_TOKENS:
        assert out == name


@given(base=plain_base,
       greek=st.sampled_from(names.GREEK_LETTERS),
       roman=st.sampled_from(list(names.ROMAN_NUMERALS_BY_VALUE.values())),
       companion=st.sampled_from(names.COMPANION_SUFFIXES),
       diminutive=st.sampled_from(names.DIMINUTIVE_PREFIXES))
def test_strip_decoration_inverts_each_single_decoration(base, greek, roman, companion, diminutive):
    for decorated in (f"{greek} {base}", f"{greek} {base} {roman}", f"{base} {companion}",
                      f"{diminutive} {base}", f"{greek} {diminutive} {base}"):
        assert nu.strip_decoration(decorated) == base, decorated


@given(base=plain_base, count=st.integers(1, nu.GREEK_ROMAN_CAPACITY - 2),
       diminutives=st.integers(1, len(names.DIMINUTIVE_PREFIXES)),
       companions=st.integers(1, len(names.COMPANION_SUFFIXES)))
@example(base="Voranthis", count=1, diminutives=2, companions=2)   # "Little Beta ...", "Petit Little ..."
@example(base="Voranthis", count=30, diminutives=1, companions=1)  # "Little Alpha Voranthis <roman>"
def test_strip_decoration_inverts_the_stacks_db_actually_builds(base, count, diminutives, companions):
    """Regression: `strip_decoration("Little Beta Foo")` used to give
    "Beta Foo" (and "Petit Little Foo" -> "Little Foo", "Foo Kin Ami" ->
    "Foo Kin"). These are the exact stacks `_db.py` builds:
    `reserve_system_name` prefixes the diminutive onto the Greek/Roman
    name, and `_rename_existing_system_for_diminutive`/
    `_rename_existing_body_for_companion` decorate the same row again on
    every later collision."""
    greek_name, _ = nu.resolve_greek_roman_collision(base, count)
    prefixes, index = [], None
    for _ in range(diminutives):
        prefix, index = nu.resolve_diminutive(index)
        prefixes.insert(0, prefix)          # each rename goes on the outside
    suffixes, index = [], None
    for _ in range(companions):
        suffix, index = nu.resolve_companion(index)
        suffixes.append(suffix)
    stacks = [
        f"{prefixes[-1]} {greek_name}",                  # reserve_system_name after a sector hit
        " ".join(prefixes + [greek_name]),               # ...then renamed again by later sectors
        " ".join(prefixes + [base]),                     # _rename_existing_system_for_diminutive, repeatedly
        " ".join([base] + suffixes),                     # _rename_existing_body_for_companion, repeatedly
    ]
    for stack in stacks:
        assert nu.strip_decoration(stack) == base, stack


@given(token=st.sampled_from(sorted(DECORATION_TOKENS)), other=st.sampled_from(["gughoe", "kelmo", "tavira"]))
@example(token="Liten", other="gughoe")   # brute force produced the plain sector name "Liten Gughoe"
@example(token="Ohana", other="kelmo")
def test_generated_words_are_never_decoration_tokens(token, other):
    """Regression: `is_name_valid` accepted decoration words ("liten",
    "klein", "ohana", "amigo", ...), so a plain generated sector name like
    "Liten Gughoe" was stripped to "Gughoe" and grouped with another base."""
    assert utils.is_name_valid(token.lower()) is False
    assert utils.is_name_valid(f"{token.lower()} {other}") is False
    assert utils.is_name_valid(other) is True  # the rejection is the token's, not the filler's


def test_split_halves_of_generated_names_are_never_decoration_tokens():
    """A long name split in two can land a decoration word as either half
    ("Kavel Ohana"); the generator must draw again rather than return it."""
    splits = iter(["kavel Ohana", "kavel Onami"])
    with mock.patch.object(utils, "is_name_valid", return_value=True), \
            mock.patch.object(utils, "split_long_word", side_effect=lambda _name: next(splits)):
        name = utils.generate_phoneme_salad_name(*NAME_SOURCES["star"])
    assert name == "Kavel Onami"


def test_resolvers_reject_negative_indices():
    with pytest.raises(ValueError):
        nu.resolve_greek_roman_collision("X", -1)
    with pytest.raises(ValueError):
        nu.resolve_diminutive(-1)
    with pytest.raises(ValueError):
        nu.resolve_companion(-100)


@given(n=st.integers(max_value=-1))
def test_resolvers_reject_any_negative_index(n):
    for call in (lambda: nu.resolve_greek_roman_collision("X", n), lambda: nu.resolve_diminutive(n),
                 lambda: nu.resolve_companion(n)):
        with pytest.raises(ValueError):
            call()


# ---------------------------------------------------------------------------
# utils name helpers
# ---------------------------------------------------------------------------

@given(name=hostile_text)
def test_split_into_syllables_is_a_partition(name):
    syllables = utils.split_into_syllables(name)
    assert "".join(syllables) == name
    assert all(syllables)
    assert (syllables == []) == (name == "")


@given(name=hostile_text)
def test_is_name_valid_never_raises_and_rejects_every_offensive_substring(name):
    result = utils.is_name_valid(name)
    assert isinstance(result, bool)
    lower = name.lower()
    if lower in names.DICTIONARY_WORDS or any(w in lower for w in names.NSFW_WORDS):
        assert result is False


@given(prefix=st.text(alphabet=string.ascii_lowercase, max_size=4),
       word=st.sampled_from(sorted(w for w in names.NSFW_WORDS if w.isalpha())),
       suffix=st.text(alphabet=string.ascii_lowercase, max_size=4))
def test_is_name_valid_rejects_embedded_offensive_words(prefix, word, suffix):
    assert utils.is_name_valid(prefix + word + suffix) is False


@given(run=st.sampled_from(["aei", "ouu", "bcd", "rst", "xyzq"]), pad=st.text(alphabet="ab", max_size=4))
def test_is_name_valid_rejects_three_in_a_row(run, pad):
    assert utils.is_name_valid(pad + run + pad) is False


@given(prefix=st.text(alphabet="bdgklmnprt", max_size=3),
       word=st.sampled_from(sorted(w for w in names.NSFW_WORDS if w.isalpha() and len(w) >= 3)),
       cut=st.integers(1, 10), suffix=st.text(alphabet="aeiou", max_size=2))
@example(prefix="ge", word="paki", cut=3, suffix="")   # generated "Gepak'I Conio"
@example(prefix="aln", word="kike", cut=1, suffix="m")  # generated "Alnik'Ikem" (via "i" + "k'ike")
def test_is_name_valid_sees_through_embedded_apostrophes(prefix, word, cut, suffix):
    """Regression: an apostrophe spliced in by a phoneme chunk ("k'",
    "pak'") hid an offensive word from the substring filter."""
    cut = min(cut, len(word) - 1)
    hidden = prefix + word[:cut] + "'" + word[cut:] + suffix
    assert utils.is_name_valid(hidden) is False
    assert utils.is_name_valid(hidden.replace("'", " ")) is False


@given(name=st.text(alphabet=string.ascii_lowercase + "'", min_size=0, max_size=40))
def test_split_long_word_properties(name):
    out = utils.split_long_word(name)
    if len(name) <= names.WORD_SIZE_MEAN or out == name:
        assert out == name
        return
    assert out.count(" ") == 1
    first, second = out.split(" ")
    assert first and second
    assert (first + second).lower() == name
    assert second[0] == second[0].upper()
    if "'" not in name:
        assert abs(len(first) - len(second)) <= 1


@given(words=st.lists(st.text(alphabet=string.ascii_lowercase, min_size=1, max_size=9), min_size=1, max_size=3))
@example(words=["ethel", "nacyon"])   # base "El Nath" -> generated "Ethel  Nacyon"
def test_split_long_word_never_doubles_an_existing_space(words):
    name = " ".join(words)
    out = utils.split_long_word(name)
    assert "  " not in out
    assert out.count(" ") <= max(1, name.count(" "))


@settings(max_examples=scaled(30))
@given(base_list=st.lists(st.text(alphabet=string.ascii_letters + "'", min_size=1, max_size=15),
                          min_size=1, max_size=5),
       prefixes=st.lists(st.text(alphabet=string.ascii_letters, min_size=1, max_size=4), min_size=1, max_size=3),
       suffixes=st.lists(st.text(alphabet=string.ascii_letters, min_size=1, max_size=4), min_size=1, max_size=3),
       allow_split=st.booleans(),
       fraction=st.floats(0.01, 1.0),
       max_length=st.one_of(st.none(), st.integers(1, 12)))
def test_generate_phoneme_salad_name_terminates_on_arbitrary_inputs(base_list, prefixes, suffixes,
                                                                      allow_split, fraction, max_length):
    """Any word lists either yield a valid name or a clean RuntimeError --
    never another exception, never a hang (attempt cap lowered here so
    even an always-rejecting input is quick)."""
    with mock.patch.object(utils, "MAX_NAME_GENERATION_ATTEMPTS", 60):
        try:
            name = utils.generate_phoneme_salad_name(base_list, prefixes, suffixes, allow_split=allow_split,
                                                     syllable_fraction=fraction, max_length=max_length)
        except RuntimeError:
            return
    assert name and name[0] == name[0].upper()
    assert utils.is_name_valid(name.replace(" ", "").lower()) or " " in name
    if max_length is not None:
        assert len(name.replace(" ", "")) <= max_length
    if not allow_split:
        assert " " not in name


def _check_generated(name, sector=False):
    lower = name.lower()
    # (A double space is possible today -- see
    # test_split_long_word_never_doubles_an_existing_space -- so words are
    # split on any whitespace run here.)
    assert name and name == name.strip()
    assert name[0].isupper(), name
    assert not any(word in lower for word in names.NSFW_WORDS), name
    for word in name.split():
        assert word and not word.startswith("'") and not word.endswith("'"), name
    if sector:
        words = name.split(" ")
        assert len(words) == 2, name
        assert all(len(w) <= 7 for w in words), name


@pytest.mark.parametrize("kind", sorted(NAME_SOURCES))
def test_bulk_generated_names_are_clean(kind):
    source = NAME_SOURCES[kind]
    for _ in range(scaled(60) * 25):
        _check_generated(utils.generate_phoneme_salad_name(*source))


def test_bulk_generated_sector_names_are_clean():
    for _ in range(scaled(60) * 25):
        _check_generated(utils.generate_sector_name(), sector=True)


def test_generation_is_bounded_even_when_every_candidate_is_rejected():
    with mock.patch.object(utils, "is_name_valid", return_value=False) as validator:
        with pytest.raises(RuntimeError):
            utils.generate_phoneme_salad_name(*NAME_SOURCES["star"])
    assert validator.call_count == utils.MAX_NAME_GENERATION_ATTEMPTS


@given(huge=st.integers(10**6, 10**40))
def test_huge_counts_never_index_out_of_range(huge):
    assume(huge >= nu.GREEK_ROMAN_CAPACITY)
    assert nu.resolve_greek_roman_collision("X", huge) == (None, None)
    assert nu.resolve_diminutive(huge) == (None, None)
    assert nu.resolve_companion(huge) == (None, None)
