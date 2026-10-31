"""
The condition-based planet and moon classes (GEN.645, GEN.647; C1, C2 and C5 of docs/design/habitability/condition-based-planet-classes.md):
`physics/planet_classes.py` and its rule file `physics/data/planet_classes.json`.

* the rule file is well formed: every letter once, the free letters kept apart, the last row of each kind takes the rest;
* every row has a body that gets that letter and nothing else (a test per rule), and the order of the rows decides;
* the function is pure and total: any stored conditions give one letter, the same every time, from an object or a mapping;
* the class sweep finds every stored letter is the derived one;
* no module but the generator's current writers assigns a body's class (C5: the list shrinks at the cut-over, GEN.646).
"""

import os
import random
import re

import pytest

from planetgen import tuning
from planetgen.physics import atmosphere, constants, gases, planet_classes as pc

EARTH = {
    "body_type": "t", "mass_kg": constants.EARTH_MASS_TO_KG, "gravity_g": 1.0, "surface_temperature_k": 288.0,
    "atmospheric_pressure_pa": 101325.0, "hydrosphere": "surface ocean", "ocean_fraction": 0.7,
    "water_mass_fraction": 2.3e-4, "flux_variation": 0.0,
    atmosphere.GAS_FIELD: {"O2": 21200.0, "N2": 78800.0, "Ar": 900.0, "CO2": 40.0},
}


def put_kpa(values, gas, kpa):
    """Sets one species' partial pressure (kPa) in `values`' own copy of the air."""
    air = dict(values[atmosphere.GAS_FIELD])
    air[gases.stored(gas).formula] = kpa * 1000.0
    values[atmosphere.GAS_FIELD] = air


def body(**changes):
    """Earth's stored conditions with `changes`; a `mass_earth` change sets the mass, `share_<gas>` a gas's partial pressure."""
    values = dict(EARTH)
    for name, value in changes.items():
        if name == "mass_earth":
            values["mass_kg"] = value * constants.EARTH_MASS_TO_KG
        elif name.startswith("share_"):
            put_kpa(values, name[len("share_"):], value * values["atmospheric_pressure_pa"] / 1000.0)
        elif name.startswith("p_") and name.endswith("_kpa"):
            put_kpa(values, name[2:-4], value)
        else:
            values[name] = value
    return values


def airless(**changes):
    values = body(atmospheric_pressure_pa=0.0, hydrosphere="dry", ocean_fraction=0.0, **changes)
    values[atmosphere.GAS_FIELD] = {}
    return values


def no_oxygen(**changes):
    values = body(**changes)
    put_kpa(values, "o2", 0.0)
    return values


CASES = {
    "T": body(body_type="g", mass_earth=6.0),
    "I": body(body_type="g", mass_earth=20.0),
    "J": body(body_type="g", mass_earth=300.0),
    "D": airless(water_mass_fraction=0.3, mass_earth=0.01),
    "S": airless(mass_earth=3.0, gravity_g=1.4),
    "C": airless(mass_earth=0.05, gravity_g=0.3),
    "A": body(share_so2=0.2, p_o2_kpa=0.0),
    "N": body(atmospheric_pressure_pa=9e6, surface_temperature_k=700.0, hydrosphere="vapour", ocean_fraction=0.0,
              p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=8800.0),
    "V": body(gravity_g=2.5, mass_earth=8.0),
    "L": body(p_co2_kpa=3.0, p_o2_kpa=10.0, gravity_g=0.9),
    "K": body(atmospheric_pressure_pa=500.0, surface_temperature_k=300.0, hydrosphere="dry", gravity_g=0.4,
              p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=0.5),
    "P": body(surface_temperature_k=250.0, hydrosphere="ice", ocean_fraction=0.4, atmospheric_pressure_pa=6000.0),
    "B": body(gravity_g=0.3, mass_earth=0.04, surface_temperature_k=450.0, hydrosphere="dry", ocean_fraction=0.0),
    "E": body(surface_temperature_k=200.0, hydrosphere="ice", ocean_fraction=0.5, atmospheric_pressure_pa=6000.0),
    "F": no_oxygen(surface_temperature_k=300.0, share_ch4=0.2, share_co2=0.5, p_n2_kpa=10.0, p_ar_kpa=0.0),
    "G": body(share_co2=0.5, p_o2_kpa=0.0),
    "O": body(ocean_fraction=0.97),
    "H": body(ocean_fraction=0.02, hydrosphere="dry"),
    "M": body(),
}


def test_class_q_is_retired_and_its_letter_reserved():
    """GEN.1135 (Boss, 6:07 PM EDT 2026-10-30): the swing of the starlight is the stored flux_variation and a flag, not a class."""
    assert "Q" not in {c.letter for c in pc.CLASSES} and "Q" not in tuning.PLANET_CLASSES and "Q" in pc.RESERVE
    assert "Q" not in tuning.HABITABLE_PLANET_CLASSES
    assert all(rule.get("quantity") != "flux_variation" for c in pc.CLASSES for rule in c.when)
    assert pc.classify(body(flux_variation=0.9)) == "M"      # a very eccentric orbit alone makes no class


def test_the_rule_file_is_well_formed():
    letters = [c.letter for c in pc.CLASSES]
    assert len({(c.letter, c.kind) for c in pc.CLASSES}) == len(letters)      # a letter once in a kind (T has a giant and a terrestrial row)
    assert set(letters) == set(tuning.PLANET_CLASSES), "a letter was lost or added: ask Boss first"
    assert not set(pc.RESERVE) & set(letters)
    assert set(pc.RESERVE) | set(letters) == set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assert set(pc.CORE_LETTERS) == set("HPOJMICDN") and set(pc.CORE_LETTERS) <= set(letters)
    for kind in (pc.TERRESTRIAL, pc.GIANT):
        rows = [c for c in pc.CLASSES if c.kind == kind]
        assert not rows[-1].when and all(c.when for c in rows[:-1])
    for cls in pc.CLASSES:
        assert cls.name and cls.sentence and cls.basis.startswith("design")


def test_every_rule_has_a_body_that_gets_its_letter():
    assert set(CASES) == set(pc.LETTERS)
    for letter, values in CASES.items():
        assert pc.classify(values) == letter, f"{letter}: {pc.rule_words(pc.describe(letter))}"


def test_order_decides_between_two_matching_rows():
    # An airless super-Earth is also a dead world: the earlier row (S) takes it.
    assert pc.classify(airless(mass_earth=3.0)) == "S"
    assert pc.classify(airless(mass_earth=0.5)) == "C"
    # A hot small world with sulfur dioxide air is volcanic, not a small hot world.
    assert pc.classify(body(gravity_g=0.3, share_so2=0.1)) == "A"
    # Giants are never read by the terrestrial rows, and the other way round.
    assert pc.classify(body(body_type="g", mass_earth=100.0, gravity_g=2.4)) == "J"


def test_thresholds_are_edges_of_their_row():
    assert pc.classify(airless()) == "C"
    assert pc.classify(body(atmospheric_pressure_pa=9.99, p_o2_kpa=0, p_n2_kpa=0, p_ar_kpa=0, p_co2_kpa=0,
                            hydrosphere="dry", ocean_fraction=0.0, surface_temperature_k=300.0)) == "C"
    assert pc.classify(body(ocean_fraction=0.8)) == "O"
    assert pc.classify(body(ocean_fraction=0.7999)) == "M"
    assert pc.classify(body(ocean_fraction=0.2)) == "M"
    assert pc.classify(body(ocean_fraction=0.1999)) == "H"
    assert pc.classify(body(surface_temperature_k=273.15)) == "M" and pc.classify(body(surface_temperature_k=273.1)) == "P"


def test_the_function_is_pure_total_and_the_same_from_an_object():
    rng = random.Random(7)
    for _ in range(3000):
        values = body(
            body_type=rng.choice("tg"), mass_earth=10 ** rng.uniform(-4, 3.6), gravity_g=rng.uniform(0.01, 5),
            surface_temperature_k=rng.uniform(2.7, 2500), atmospheric_pressure_pa=rng.choice([0.0, 10 ** rng.uniform(-3, 8)]),
            ocean_fraction=rng.random(), water_mass_fraction=10 ** rng.uniform(-8, -0.3),
            hydrosphere=rng.choice(["dry", "vapour", "ice", "ice-covered ocean", "surface ocean", "brine seep", "hycean", None]),
            flux_variation=rng.choice([0.0, rng.random()]))
        for gas in ("o2", "co2", "n2", "ar", "h2", "h2o", "ch4", "h2s", "so2", "co"):
            put_kpa(values, gas, rng.choice([0.0, values["atmospheric_pressure_pa"] / 1000.0 * rng.random() / 6]))
        first = pc.classify(values)
        assert first in pc.LETTERS and first == pc.classify(dict(values))
        assert (pc.describe(first).kind == values["body_type"])

    class Planet:
        def __init__(self):
            self.body_type, self.mass, self.gravity, self.surface_temperature = "t", EARTH["mass_kg"], 1.0, 288.0
            self.atmospheric_pressure, self.hydrosphere, self.ocean_fraction = 101325.0, "surface ocean", 0.7
            self.water_mass_fraction, self.flux_variation = 2.3e-4, 0.0
            self.gas_partials_pa = {"O2": 21200.0, "N2": 78800.0}

    assert pc.classify(Planet()) == "M"


def test_a_body_without_what_the_rules_read_is_refused():
    with pytest.raises(pc.ClassError):
        pc.classify({"body_type": "t"})
    with pytest.raises(pc.ClassError):
        pc.classify({"body_type": "x", "mass_kg": 1.0})
    values = dict(EARTH, surface_temperature_k=None)
    with pytest.raises(pc.ClassError):
        pc.classify(values)  # cold air worlds read the temperature
    with pytest.raises(pc.ClassError):
        pc.describe("R")


def test_a_class_follows_its_conditions():
    values = body()
    assert pc.classify(values) == "M"
    values["surface_temperature_k"] = 250.0
    assert pc.classify(values) == "P"
    values.update(atmospheric_pressure_pa=0.0, **{atmosphere.GAS_FIELD: {}})
    assert pc.classify(values) == "C"


def test_rules_read_in_plain_words():
    assert pc.rule_words(pc.describe("S")) == ["surface pressure below 10 Pa", "mass at least 1.5 Earth masses"]
    assert pc.rule_words(pc.describe("N"))[0] == "surface pressure at least 1,000,000 Pa"
    assert any(word.startswith("oxygen share of the air at least 3 percent") for word in pc.rule_words(pc.describe("L")))
    assert pc.rule_words(pc.describe("M")) == []


def test_the_class_sweep_finds_every_stored_letter_is_the_derived_one():
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    import class_sweep_cli

    bodies = class_sweep_cli.sample_bodies(12, seed=3, habitable_share=0.4)
    assert sum(not is_moon for _, is_moon in bodies) > 20 and any(is_moon for _, is_moon in bodies)
    # GEN.646: the letter is derived from the stored conditions once the system is final, so classifying again changes nothing.
    assert class_sweep_cli.changed_letters(bodies) == []
    report = class_sweep_cli.format_report(bodies, 12, 3)
    assert "Share of each letter" in report and "Planets by mass" in report


# C5: only the generator's own writers assign a class.
CLASS_WRITERS = {
    "generation/phenomena/rogue.py", "generation/planet.py", "physics/planets.py",
}
"""set: The modules that assign `.planet_class`: a planet or moon is classified by `physics/planets.apply_class` (and cleared there and
in `Planet.__init__` while it is built), a rogue by `rogue.rogue_class`. A new writer fails this test."""


def test_only_the_known_writers_assign_a_class():
    root = os.path.join(os.path.dirname(__file__), "..", "planetgen")
    found = set()
    for folder, _, files in os.walk(root):
        for name in files:
            if name.endswith(".py"):
                path = os.path.join(folder, name)
                with open(path, encoding="utf-8") as handle:
                    if re.search(r"(?<!\})\.planet_class\s*=[^=]", handle.read()):
                        found.add(os.path.relpath(path, root).replace(os.sep, "/"))
    assert found == CLASS_WRITERS, f"new or removed writers of a body's class: {sorted(found ^ CLASS_WRITERS)}"


def test_gen_1142_lines():
    """GEN.1142: the lines of the train 6 epoch, each on both sides of its edge (registry 6.121)."""
    hydrogen = body(mass_earth=5.0, atmospheric_pressure_pa=1.0e6, surface_temperature_k=600.0)
    put_kpa(hydrogen, "h2", 900.0)
    assert pc.classify(hydrogen) == "T"                                  # a hydrogen-rich terrestrial-kind body under 10 Earth masses
    assert pc.classify(dict(hydrogen, mass_kg=12.0 * constants.EARTH_MASS_TO_KG)) != "T"      # not under 10 Earth masses
    assert pc.classify(dict(hydrogen, mass_kg=12.0 * constants.EARTH_MASS_TO_KG)) != "N"      # and hydrogen-dominated, so not a Venus-like air
    assert pc.classify(CASES["N"]) == "N" and pc.classify(dict(CASES["N"], surface_temperature_k=399.0)) != "N"      # hot
    assert pc.describe("T").kind == pc.GIANT and {c.kind for c in pc.CLASSES if c.letter == "T"} == {pc.GIANT, pc.TERRESTRIAL}
    assert "E" in {c.letter for c in pc.CLASSES} and "E" not in pc.RESERVE                    # E is the deep-frozen world now
    assert pc.classify(body(surface_temperature_k=233.14)) == "E" and pc.classify(body(surface_temperature_k=233.15)) == "P"
    thin = body(atmospheric_pressure_pa=2999.0, surface_temperature_k=300.0, p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=0.2)
    assert pc.classify(thin) == "K"
    assert pc.classify(dict(thin, surface_temperature_k=272.0)) == "P"          # the cold part of the old K is glaciated now
    assert pc.classify(dict(thin, surface_temperature_k=321.0)) != "K"
    mask = body(p_co2_kpa=3.0, p_o2_kpa=10.0)
    assert pc.classify(mask) == "L"
    assert pc.classify(dict(mask, atmospheric_pressure_pa=29000.0)) != "L" and pc.classify(dict(mask, atmospheric_pressure_pa=310000.0)) != "L"
    put_kpa(mask, "o2", 20.0)
    assert pc.classify(mask) != "L"                                      # too much oxygen to need a mask
