"""
sectorGen.generate_sector_name regression tests.

Covers the two-word invariant: `generate_sector_name` joins two independent
`generate_phoneme_salad_name` calls into one name, and that inner function
can itself split a long result into two words (`split_long_word`) -- left
unchecked, that silently produced 3-4 word sector names instead of 2.

Run with: pytest tests/test_sector_gen.py
"""
from sectorGen import generate_sector_name

TRIALS = 500


def test_generate_sector_name_is_always_two_words():
    for _ in range(TRIALS):
        name = generate_sector_name()
        words = name.split(" ")
        assert len(words) == 2, f"expected exactly 2 words, got {len(words)}: {name!r}"
        assert all(word for word in words), f"unexpected empty word in {name!r}"
