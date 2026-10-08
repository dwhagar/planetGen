"""planetgen.names.gated_phoneme_codec (GEN.120): hex IDs to pronounceable
words and back. Names are seed-reproducible, so the golden words pin the
codec: a change that moves them renames objects and needs a naming key."""

import itertools
import random

import pytest

from planetgen.names.gated_phoneme_codec import GatedPhonemeCodec, PhonemeDecodeError

codec = GatedPhonemeCodec()


def _random_hex(rng, length):
    return "".join(rng.choice("0123456789ABCDEF") for _ in range(length))


# Fixed (ID, domain) pairs and the words they must keep giving.
GOLDEN = [
    ("0123456789ABCDEF012", "rogue_planet", "vinparkuf bunsanmee hibunsan"),
    ("0123456789ABCDEF012", "nebula", "favinpark hibunsan meehibun"),
    ("3A7F00B2", "rogue_planet", "jeru fogay"),
    ("10D00B002C7", "comet", "visujai largeek"),
    ("FFFFFFFFFFFFFFFFFFF", "constellation", "rintunzad fohaikoo morpenson"),
    ("0000000000000000000", "quasar", "digujeel morpenson vaybefo"),
]


@pytest.mark.parametrize("identifier, domain, words", GOLDEN)
def test_fixed_ids_give_pinned_words(identifier, domain, words):
    assert codec.encode(identifier, domain) == words
    assert codec.decode(words, domain, length=len(identifier)) == identifier


def test_wordiness_presets_pin_the_word_count():
    assert codec.encode("3A7F00B2C9", "comet") == "muvooj keefaiz"
    assert codec.encode("3A7F00B2C9", "comet", wordiness="balanced") == "muvooj keefaiz"
    assert codec.encode("3A7F00B2C9", "comet", wordiness="verbose") == "muvoo jeek faiz"
    assert GatedPhonemeCodec("verbose").encode("3A7F00B2C9", "comet") == "muvoo jeek faiz"
    with pytest.raises(ValueError, match="Unknown wordiness"):
        codec.encode("3A7F00B2C9", wordiness="chatty")


@pytest.mark.parametrize("wordiness", ["compact", "balanced", "verbose", 5])
@pytest.mark.parametrize("style", ["auto", "words", "sentences", "paragraphs"])
def test_encode_then_decode_returns_the_same_hex(wordiness, style):
    rng = random.Random(120)
    for _ in range(600):
        length = rng.randint(1, 60)
        identifier = _random_hex(rng, length)
        domain = rng.choice(["rogue_planet", "nebula", "star", 7])
        phrase = codec.encode(identifier, domain, wordiness=wordiness, style=style)
        assert codec.decode(phrase, domain, wordiness=wordiness, length=length) == identifier, (
            identifier, domain, phrase)


def test_the_ids_the_game_names_round_trip():
    """GEN.64's IDs are 19 hex digits."""
    rng = random.Random(64)
    for _ in range(3000):
        identifier = "%019X" % rng.getrandbits(76)
        phrase = codec.encode(identifier, "quasar")
        assert codec.decode_safe(phrase, "quasar", length=19) == (identifier, None)


def test_short_ids_need_no_length_and_a_word_that_reads_two_ways_asks_for_it():
    for identifier in ("A", "3F", "0BC"):
        assert codec.decode(codec.encode(identifier, "star"), "star") == identifier
    # "bar" is the digits 08 or 00B, so the length settles which.
    assert codec._word_readings("bar").keys() == {2, 3}
    ambiguous = next(
        h for h in (f"{n:09X}" for n in itertools.count(0x1234567))
        if codec.decode_safe(codec.encode(h, "d"), "d")[0] is None
    )
    assert ambiguous == "001234568"
    phrase = codec.encode(ambiguous, "d")
    result, error = codec.decode_safe(phrase, "d")
    assert result is None and "pass the identifier's length" in error
    assert codec.decode(phrase, "d", length=9) == ambiguous


def test_two_ids_of_the_same_length_never_share_a_name():
    for length in (1, 2, 3, 4):
        names = {codec.encode("".join(digits), "star") for digits in
                 itertools.product("0123456789ABCDEF", repeat=length)}
        assert len(names) == 16 ** length


def test_domains_give_the_same_id_different_words():
    rng = random.Random(7)
    domains = ["rogue_planet", "nebula", "quasar", "comet", "asteroid_field", "constellation"]
    differing = 0
    for _ in range(200):
        identifier = _random_hex(rng, 19)
        words = {domain: codec.encode(identifier, domain) for domain in domains}
        differing += len(set(words.values())) > 1
        for domain, phrase in words.items():
            assert codec.decode(phrase, domain, length=19) == identifier
    assert differing == 200
    # The same words under the wrong domain are another ID.
    phrase = codec.encode("0123456789ABCDEF012", "nebula")
    assert codec.decode(phrase, "quasar", length=19) != "0123456789ABCDEF012"


def test_a_domain_is_a_name_or_a_number_and_the_permutation_is_stable():
    assert codec._get_domain_params("rogue_planet") == (13, 14)
    assert codec._get_domain_params(7) == (1, 2)
    assert codec._get_domain_params(7) == codec._get_domain_params("7")


def test_raw_text_round_trips_as_utf8():
    phrase = codec.encode("Hello, Wörld", "label", is_raw_text=True)
    assert codec.decode(phrase, "label", as_raw_text=True, length=2 * len("Hello, Wörld".encode())) == "Hello, Wörld"
    assert codec.encode("Hi", "t", is_raw_text=True) == "se mar"


def test_words_follow_the_gating_grammar():
    rng = random.Random(3)
    vowels = set("aeiou")
    for _ in range(500):
        for word in codec.encode(_random_hex(rng, 19), "star", style="words").split():
            assert word[0] not in vowels, word
            assert "aa" not in word and "ii" not in word and "uu" not in word


def test_long_phrases_are_split_into_sentences_and_paragraphs():
    identifier = "0123456789ABCDEF" * 3
    sentences = codec.encode(identifier, "d", style="sentences")
    assert sentences == "Seebanhun miseeban hunmisee banhunmi seebanhun miseeban. Hunmisee banhunmi."
    assert codec.decode(sentences, "d", length=48) == identifier
    many = codec.encode("0123456789ABCDEF" * 20, "d", style="paragraphs")
    assert many.count("\n\n") == 2
    assert codec.decode(many, "d", length=320) == "0123456789ABCDEF" * 20


@pytest.mark.parametrize("phrase", ["zzzq", "qua", "bar!!x1"])
def test_bad_words_fail_cleanly(phrase):
    message = codec.decode(phrase, "d")
    assert message.startswith("Decode Error")
    assert codec.decode_safe(phrase, "d") == (None, message)
    with pytest.raises(PhonemeDecodeError):
        codec.decode(phrase, "d", strict=True)


def test_a_phrase_that_is_no_encoding_is_refused_with_the_length_it_was_given():
    with pytest.raises(PhonemeDecodeError, match="not the encoding of any identifier of 7 digits"):
        codec.decode("jeru fogay", "rogue_planet", strict=True, length=7)


def test_empty_input_is_empty():
    assert codec.encode("  ", "d") == ""
    assert codec.decode("  ", "d") == ""
