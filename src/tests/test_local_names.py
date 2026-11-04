# tests/test_local_names.py

"""
Local names from the Object ID (GEN.1255, `planetgen/names/naming_key.py`): the serial and body number of an Object ID spelled alone
as two words through a keyed bijective scramble of the serial's count bits, with the birth sector shown beside the name.
No database is needed.
"""

import random

import pytest

from planetgen.galaxy import object_uid
from planetgen.names import naming_key

KEY = "00AB12CD"
OTHER_KEY = "5E11C0DE"


def _uid(shell, band, wedge, serial, body=0, kind=object_uid.SERIAL_GENERATED):
    return object_uid.pack(shell, band, wedge, kind, serial, body)


def _first_word(name):
    return name.split()[0]


def test_the_codec_version_counts_local_names():
    assert naming_key.CODEC_VERSION == 2


def test_golden_local_names():
    """The words under one key; they change only with `CODEC_VERSION`."""
    golden = [
        (_uid(0, 0, 0, 0), "Soojork Dusonv"),
        (_uid(0, 0, 0, 1), "Jiteh Hoosonv"),
        (_uid(0, 0, 0, 2), "Zanforf Lansonv"),
        (_uid(5, 3, 1000, 12345, 7), "Zaidonr Jaisonn"),
        (_uid(65535, -102900, 411000, 67108863, 4095, kind=object_uid.SERIAL_RUNTIME), "Bingonr Bezad"),
        (_uid(100, -5, 77, 1, 3, kind=object_uid.SERIAL_FIELD), "Lileez Noosonj"),
        (_uid(100, -5, 77, 1, 3, kind=object_uid.SERIAL_REGION), "Rileez Noosonj"),
    ]
    for uid, words in golden:
        assert naming_key.local_name(uid, KEY) == words, object_uid.format_id(uid)


def test_a_local_name_is_two_words_and_carries_the_sector_beside_it():
    uid = _uid(5, 3, 1000, 12345, 7)
    name = naming_key.local_name(uid, KEY)
    assert len(name.split()) == 2 and name == name.title()
    assert naming_key.sector_text(uid) == "0000B0001803E8" == object_uid.format_id(uid).split("-")[0]
    assert naming_key.local_name_with_sector(uid, KEY) == f"{name}, sector 0000B0001803E8"


@pytest.mark.parametrize("key", [KEY, OTHER_KEY, "00000000", "FFFFFFFF"])
def test_the_scramble_is_a_bijection_and_reverses(key):
    rng = random.Random(1)
    for sector in (0, 1, 0x0000B0001803E8, (1 << 56) - 1, rng.getrandbits(56)):
        a, b = naming_key._scramble_params(key, sector)
        assert a % 2 == 1
        counts = [rng.getrandbits(26) for _ in range(300)] + [0, 1, (1 << 26) - 1]
        for count in counts:
            scrambled = naming_key.scramble(count, sector, key)
            assert 0 <= scrambled < 1 << 26
            assert naming_key.unscramble(scrambled, sector, key) == count
        # A bijection on a sample: no two of 20,000 consecutive counts land on the same value.
        assert len({naming_key.scramble(count, sector, key) for count in range(20000)}) == 20000


def test_serials_of_a_sector_do_not_share_their_first_word():
    """Serials 0, 1, 2 ... spelled directly share their first word; scrambled, the first `SPREAD_SERIALS` of every sector have their own."""
    sectors = [(0, 0, 0), (1, -1, 5), (19, 200, 3000), (65535, -102900, 411000), (4000, 12, 70000)]
    for key in (KEY, OTHER_KEY):
        for shell, band, wedge in sectors:
            names = [naming_key.local_name(_uid(shell, band, wedge, serial), key) for serial in range(naming_key.SPREAD_SERIALS)]
            assert len({_first_word(name) for name in names}) == naming_key.SPREAD_SERIALS, (key, shell, band, wedge)
            assert len(set(names)) == naming_key.SPREAD_SERIALS


def test_local_names_are_unique_inside_a_sector_across_kinds_and_bodies():
    seen = {}
    for kind in object_uid.SERIAL_KINDS:
        for serial in range(0, 700):
            for body in (0, 1, 2, 3, 4095):
                uid = _uid(7, -3, 90, serial, body, kind)
                name = naming_key.local_name(uid, KEY)
                assert name not in seen, (name, object_uid.format_id(uid), object_uid.format_id(seen[name]))
                seen[name] = uid
    assert len(seen) == 4 * 700 * 5


def test_the_bodies_of_a_system_share_the_first_word_and_differ_in_the_second():
    names = [naming_key.local_name(_uid(7, -3, 90, 41, body), KEY) for body in range(0, 40)]
    assert len({_first_word(name) for name in names}) == 1
    assert len(set(names)) == 40


@pytest.mark.parametrize("key", [KEY, OTHER_KEY])
def test_a_name_and_its_sector_give_back_the_object_id(key):
    rng = random.Random(7)
    cases = [(0, 0, 0, 0, 0, 0), (65535, -102900, 411000, (1 << 26) - 1, 4095, 3), (19, 5, 5, 1, 0, 1), (3, -1, 4, 77, 12, 2)]
    for _ in range(300):
        cases.append((rng.randrange(1 << 16), rng.randrange(-100000, 100000), rng.randrange(400000), rng.randrange(1 << 26),
                      rng.randrange(1 << 12), rng.randrange(4)))
    for shell, band, wedge, serial, body, kind in cases:
        uid = _uid(shell, band, wedge, serial, body, kind)
        name = naming_key.local_name(uid, key)
        sector = naming_key.sector_text(uid)
        assert naming_key.uid_of_local_name(name, sector, key) == uid
        assert naming_key.uid_of_local_name(name.lower(), int(sector, 16), key) == uid
        assert naming_key.uid_of_local_name(f"  {name}  ", f"sector {sector}", key) == uid


def test_a_name_in_another_sector_is_another_object():
    uid = _uid(7, -3, 90, 1234, 5)
    name = naming_key.local_name(uid, KEY)
    assert naming_key.uid_of_local_name(name, naming_key.sector_text(uid), KEY) == uid
    assert naming_key.uid_of_local_name(name, naming_key.sector_text(_uid(7, -3, 91, 0)), KEY) != uid


def test_a_name_does_not_carry_the_words_of_its_sector():
    """The object's neighbours in space do not share leading words: the first 200 serials of two neighbouring sectors name apart."""
    first = {naming_key.local_name(_uid(7, -3, 90, serial), KEY) for serial in range(200)}
    second = {naming_key.local_name(_uid(7, -3, 91, serial), KEY) for serial in range(200)}
    assert not first & second
    assert len({_first_word(name) for name in first | second}) > 390


def test_a_key_change_renames_every_local_name():
    uid = _uid(7, -3, 90, 1234, 5)
    names = {naming_key.local_name(uid, f"{number:08X}") for number in range(40)}
    assert len(names) > 30
    # The Object ID is the identity: the new key reads the new name back to the same ID.
    for number in (3, 17):
        key = f"{number:08X}"
        assert naming_key.uid_of_local_name(naming_key.local_name(uid, key), naming_key.sector_text(uid), key) == uid


def test_what_cannot_be_read_is_a_value_error():
    with pytest.raises(ValueError):
        naming_key.uid_of_local_name("Soojork Dusonv", "00FA0FFFE9E24", KEY)  # a sector is 14 hex digits
    with pytest.raises(ValueError):
        naming_key.uid_of_local_name("Soojork Dusonv", 1 << 56, KEY)
    with pytest.raises(ValueError):
        naming_key.uid_of_local_name("xyzzy", "00001000000000", KEY)
    with pytest.raises(ValueError):
        naming_key.local_name(_uid(0, 0, 0, 0), "not a key")
