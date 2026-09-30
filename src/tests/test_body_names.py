# tests/test_body_names.py

"""
Tests for `stellarObjects/bodyNames.py` -- stars, planets and moons named
from their system (`Voranthis I`, `Voranthis IIa`, a binary's
`Voranthis Kelmoor`) -- and for `StarSystem` applying it at generation.

Run with: pytest tests/test_body_names.py
"""

import pytest

from stellarObjects.bodyNames import moon_letters, rename_prefix, to_roman
from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem


def test_to_roman():
    assert [to_roman(n) for n in (1, 2, 3, 4, 5, 9, 10, 14, 40, 49, 90, 400, 1994)] == [
        "I", "II", "III", "IV", "V", "IX", "X", "XIV", "XL", "XLIX", "XC", "CD", "MCMXCIV",
    ]
    with pytest.raises(ValueError):
        to_roman(0)


def test_moon_letters():
    assert [moon_letters(i) for i in (0, 1, 25, 26, 27, 51, 52)] == ["a", "b", "z", "aa", "ab", "az", "ba"]


def test_rename_prefix():
    assert rename_prefix("Voranthis", "Voranthis", "Sol") == "Sol"
    assert rename_prefix("Voranthis IIa", "Voranthis", "Sol") == "Sol IIa"
    assert rename_prefix("Voranthis Kelmoor I", "Voranthis", "Alpha Voranthis") == "Alpha Voranthis Kelmoor I"
    assert rename_prefix("Voranthisa II", "Voranthis", "Sol") is None
    assert rename_prefix("Terra", "Voranthis", "Sol") is None
    assert rename_prefix(None, "Voranthis", "Sol") is None


def _system(binary=False, wide=None):
    """A generated system with at least one real planet on each star that
    can have one, retried since generation is random."""
    for _ in range(60):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = binary
        cfg.WIDE_BINARY = wide
        system = StarSystem(system_config=cfg)
        lists = [system.planets] + ([system.secondary_planets] if system.binary_type == "wide" else [])
        if all(any(p.body_type != "a" for p in bodies) for bodies in lists):
            if any(p.moons for bodies in lists for p in bodies if p.body_type != "a"):
                return system
    raise AssertionError("could not generate a system with planets and moons")


def _assert_numbered(prefix, bodies):
    planets = [p for p in bodies if p.body_type != "a"]
    for number, planet in enumerate(planets, start=1):
        assert planet.name == f"{prefix} {to_roman(number)}"
        assert [m.name for m in planet.moons] == [
            f"{prefix} {to_roman(number)}{moon_letters(i)}" for i in range(len(planet.moons))
        ]
    assert all(not hasattr(p, "name") or p.name is None for p in bodies if p.body_type == "a")


def _all_names(system):
    names = [star.name for star in system.stars]
    for planet in system.planets + system.secondary_planets:
        if planet.body_type != "a":
            names.append(planet.name)
            names.extend(m.name for m in planet.moons)
    return names


def test_single_star_system_numbers_planets_and_lettered_moons():
    system = _system()
    assert system.binary_type is None
    assert system.star.name == system.name
    _assert_numbered(system.name, system.planets)
    assert str(system).startswith(f"= {system.name} =")


def test_close_binary_stars_follow_the_system_name_and_planets_orbit_both():
    system = _system(binary=True, wide=False)
    assert system.binary_type == "close"
    primary, secondary = system.primary_star.name, system.secondary_star.name
    assert primary.startswith(system.name + " ") and secondary.startswith(system.name + " ")
    assert primary != secondary
    assert system.star.name == system.name  # the pair's proxy titles the page
    _assert_numbered(system.name, system.planets)


def test_wide_binary_planets_are_numbered_after_their_own_star():
    system = _system(binary=True, wide=True)
    assert system.binary_type == "wide"
    assert system.primary_star.name.startswith(system.name + " ")
    assert system.secondary_star.name.startswith(system.name + " ")
    _assert_numbered(system.primary_star.name, system.planets)
    _assert_numbered(system.secondary_star.name, system.secondary_planets)
    # The page is titled with the system, not the primary star.
    assert str(system).startswith(f"= {system.name} =")


@pytest.mark.parametrize("binary,wide", [(False, None), (True, False), (True, True)])
def test_no_name_uses_prime_or_a_letter_suffix(binary, wide):
    system = _system(binary=binary, wide=wide)
    names = _all_names(system)
    assert len(names) == len(set(names))
    for name in names:
        assert "Prime" not in name.split()
        assert name.split()[-1] not in ("A", "B")


def test_assign_names_follows_a_new_system_name():
    system = _system(binary=True, wide=True)
    words = [star.name.split()[-1] for star in system.stars]
    system.assign_names("Beta Rigel")
    assert system.name == "Beta Rigel"
    assert [star.name for star in system.stars] == [f"Beta Rigel {word}" for word in words]
    _assert_numbered(f"Beta Rigel {words[0]}", system.planets)


def test_serialization_round_trip_keeps_the_system_name():
    system = _system(binary=True, wide=True)
    restored = StarSystem.from_dict(system.to_dict())
    assert restored.name == system.name
    assert restored.primary_star.name == system.primary_star.name
    assert [getattr(p, "name", None) for p in restored.planets] == [getattr(p, "name", None) for p in system.planets]


def test_star_words_are_always_one_word(monkeypatch):
    """Regression: STAR_NAMES' "El Nath" could come through as a two-word
    star word, so a star's name no longer ended in its own word."""
    from stellarObjects import bodyNames

    draws = iter(["El Nath", "Vega"])
    monkeypatch.setattr(bodyNames, "generate_phoneme_salad_name", lambda *args, **kwargs: next(draws))
    assert bodyNames.generate_star_word() == "Vega"
