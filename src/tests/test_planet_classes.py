"""
The condition-based planet and moon classes (GEN.645, GEN.647; C1, C2 and C5 of docs/design/habitability/condition-based-planet-classes.md):
`physics/planet_classes.py` and its rule file `physics/data/planet_classes.json`.

* the rule file is well formed: every letter once, the free letters kept apart, the last row of each kind takes the rest;
* every row has a body that gets that letter and nothing else (a test per rule), and the order of the rows decides;
* the function is pure and total: any stored conditions give one letter, the same every time, from an object or a mapping;
* the shadow sweep keeps most of today's drawn letters (the rules stay as close to today's classes as they can);
* no module but the generator's current writers assigns a body's class (C5: the list shrinks at the cut-over, GEN.646).
"""

import os
import random
import re

import pytest

from planetgen import tuning
from planetgen.physics import constants, planet_classes as pc

EARTH = {
    "body_type": "t", "mass_kg": constants.EARTH_MASS_TO_KG, "gravity_g": 1.0, "surface_temperature_k": 288.0,
    "atmospheric_pressure_pa": 101325.0, "hydrosphere": "surface ocean", "ocean_fraction": 0.7,
    "water_mass_fraction": 2.3e-4, "flux_variation": 0.0,
    "p_o2_kpa": 21.2, "p_n2_kpa": 78.8, "p_ar_kpa": 0.9, "p_co2_kpa": 0.04,
}


def body(**changes):
    """Earth's stored conditions with `changes`; a `mass_earth` change sets the mass, `share_<gas>` a gas's partial pressure."""
    values = dict(EARTH)
    for name, value in changes.items():
        if name == "mass_earth":
            values["mass_kg"] = value * constants.EARTH_MASS_TO_KG
        elif name.startswith("share_"):
            values[f"p_{name[len('share_'):]}_kpa"] = value * values["atmospheric_pressure_pa"] / 1000.0
        else:
            values[name] = value
    return values


def airless(**changes):
    values = body(atmospheric_pressure_pa=0.0, hydrosphere="dry", ocean_fraction=0.0, **changes)
    for field in [k for k in values if k.startswith("p_")]:
        values[field] = 0.0
    return values


def no_oxygen(**changes):
    values = body(**changes)
    values["p_o2_kpa"] = 0.0
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
    "Q": body(flux_variation=0.5),
    "V": body(gravity_g=1.5, mass_earth=4.0),
    "L": body(atmospheric_pressure_pa=1500.0, share_o2=0.25, gravity_g=0.8, surface_temperature_k=230.0,
              hydrosphere="ice-covered ocean", p_n2_kpa=0.0, p_ar_kpa=0.0),
    "K": body(atmospheric_pressure_pa=500.0, surface_temperature_k=190.0, hydrosphere="ice", gravity_g=0.4,
              p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=0.5),
    "P": body(surface_temperature_k=175.0, hydrosphere="ice", ocean_fraction=0.4, atmospheric_pressure_pa=6000.0),
    "B": body(gravity_g=0.3, mass_earth=0.04, surface_temperature_k=450.0, hydrosphere="dry", ocean_fraction=0.0),
    "E": no_oxygen(surface_temperature_k=340.0, share_ch4=0.25, share_h2o=0.4, p_n2_kpa=5.0, p_ar_kpa=0.0),
    "F": no_oxygen(surface_temperature_k=300.0, share_ch4=0.2, share_co2=0.5, p_n2_kpa=10.0, p_ar_kpa=0.0),
    "G": body(share_co2=0.5, p_o2_kpa=10.0, surface_temperature_k=260.0),
    "O": body(ocean_fraction=0.97),
    "H": body(ocean_fraction=0.02, hydrosphere="dry"),
    "M": body(),
}


def test_the_rule_file_is_well_formed():
    letters = [c.letter for c in pc.CLASSES]
    assert len(letters) == len(set(letters))
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
    assert pc.classify(body(ocean_fraction=0.9)) == "O"
    assert pc.classify(body(ocean_fraction=0.8999)) == "M"
    assert pc.classify(body(ocean_fraction=0.15)) == "M"
    assert pc.classify(body(ocean_fraction=0.1499)) == "H"


def test_the_function_is_pure_total_and_the_same_from_an_object():
    rng = random.Random(7)
    for _ in range(3000):
        values = body(
            body_type=rng.choice("tg"), mass_earth=10 ** rng.uniform(-4, 3.6), gravity_g=rng.uniform(0.01, 5),
            surface_temperature_k=rng.uniform(2.7, 2500), atmospheric_pressure_pa=rng.choice([0.0, 10 ** rng.uniform(-3, 8)]),
            ocean_fraction=rng.random(), water_mass_fraction=10 ** rng.uniform(-8, -0.3),
            hydrosphere=rng.choice(["dry", "vapour", "ice", "ice-covered ocean", "surface ocean", "hycean", None]),
            flux_variation=rng.choice([0.0, rng.random()]))
        for gas in ("o2", "co2", "n2", "ar", "h2", "h2o", "ch4", "h2s", "so2", "co"):
            values[f"p_{gas}_kpa"] = rng.choice([0.0, values["atmospheric_pressure_pa"] / 1000.0 * rng.random() / 6])
        first = pc.classify(values)
        assert first in pc.LETTERS and first == pc.classify(dict(values))
        assert (pc.describe(first).kind == values["body_type"])

    class Planet:
        def __init__(self):
            self.body_type, self.mass, self.gravity, self.surface_temperature = "t", EARTH["mass_kg"], 1.0, 288.0
            self.atmospheric_pressure, self.hydrosphere, self.ocean_fraction = 101325.0, "surface ocean", 0.7
            self.water_mass_fraction, self.flux_variation = 2.3e-4, 0.0
            self.p_o2_kpa, self.p_n2_kpa = 21.2, 78.8

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
    values["surface_temperature_k"] = 150.0
    assert pc.classify(values) == "P"
    values.update(atmospheric_pressure_pa=0.0, p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=0.0)
    assert pc.classify(values) == "C"


def test_rules_read_in_plain_words():
    assert pc.rule_words(pc.describe("S")) == ["surface pressure below 10 Pa", "mass at least 1.5 Earth masses"]
    assert pc.rule_words(pc.describe("N")) == ["surface pressure at least 1,000,000 Pa"]
    assert pc.rule_words(pc.describe("L"))[1] == "oxygen share of the air at least 10 percent"
    assert pc.rule_words(pc.describe("M")) == []


def test_the_shadow_sweep_keeps_most_drawn_letters():
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    import class_shadow_cli

    bodies = class_shadow_cli.sample_bodies(25, seed=3, habitable_share=0.4)
    pairs = class_shadow_cli.compare(bodies)
    total = sum(pairs.values())
    same = sum(n for (drawn, derived), n in pairs.items() if drawn == derived)
    assert total > 100 and same / total > 0.9
    report = class_shadow_cli.format_report(pairs, 25, 3)
    assert "keep their letter" in report and "Confusion table" in report


# C5: only the generator's own writers assign a class.
CLASS_WRITERS = {
    "admin/edits.py", "generation/phenomena/rogue.py", "generation/planet.py", "physics/planets.py",
}
"""set: The modules that assign `.planet_class` today. Cut-over (GEN.646 to GEN.647) shrinks it to `physics/planet_classes.py`'s
caller; a new writer fails this test."""


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
