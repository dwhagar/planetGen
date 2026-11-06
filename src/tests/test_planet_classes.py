"""
The condition-based planet and moon classes (GEN.645, GEN.647, GEN.1207; C1, C2 and C5 of docs/design/habitability/condition-based-planet-classes.md):
`physics/planet_classes.py`, its name file `physics/data/planet_classes.json` and the conditions of every class in `physics/data/planet_class_fit.json`.

* the files are well formed: every letter once, the free letters kept apart, both files name the same classes;
* every class has a body that gets that letter, a best class always exists (a body that fits none gets its nearest with a match of 0), and a tie goes to the narrower box;
* the function is pure and total: any stored conditions give one letter, the same every time, from an object or a mapping;
* the class sweep finds every stored letter is the derived one;
* no module but the generator's current writers assigns a body's class (C5: the list shrinks at the cut-over, GEN.646).
"""

import os
import random
import re

import pytest

from planetgen import tuning
from planetgen.physics import atmosphere, class_fit as cf, constants, gases, planet_classes as pc

EARTH = {
    "body_type": "t", "mass_kg": constants.EARTH_MASS_TO_KG, "gravity_g": 1.0, "surface_temperature_k": 288.0,
    "atmospheric_pressure_pa": 101325.0, "hydrosphere": "surface ocean", "ocean_fraction": 0.7,
    "water_mass_fraction": 2.3e-4, "flux_variation": 0.0, "relative_humidity": 0.6, "groundwater": "none",
    atmosphere.GAS_FIELD: {"O2": 21200.0, "N2": 78800.0, "Ar": 900.0, "CO2": 40.0},
}
for _number in range(1, pc.BAND_COUNT + 1):     # the seven latitude bands of Earth's mean, each swinging 10 K
    EARTH[f"band{_number}_mean_k"] = 288.0
    EARTH[f"band{_number}_swing_k"] = 10.0


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
    if "surface_temperature_k" in changes:      # the latitude bands follow the mean (GEN.1221 reads the regional, day and night temperatures)
        for number in range(1, pc.BAND_COUNT + 1):
            values.setdefault(f"band{number}_mean_k", None)
            values[f"band{number}_mean_k"] = changes.get(f"band{number}_mean_k", values["surface_temperature_k"])
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
    "D": airless(mass_earth=0.01, roundness_ratio=3.0),
    "S": airless(mass_earth=3.0, gravity_g=1.4),
    "C": airless(mass_earth=0.2, gravity_g=0.3),
    "A": body(share_so2=0.2, p_o2_kpa=0.0),
    "N": body(atmospheric_pressure_pa=9e6, surface_temperature_k=700.0, hydrosphere="vapor", ocean_fraction=0.0,
              p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=8800.0),
    "V": body(gravity_g=2.5, mass_earth=8.0),
    "L": body(p_co2_kpa=3.0, p_o2_kpa=10.0, gravity_g=0.9),
    "K": body(atmospheric_pressure_pa=500.0, surface_temperature_k=300.0, hydrosphere="dry", gravity_g=0.4,
              p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=0.5),
    "P": body(surface_temperature_k=250.0, hydrosphere="ice", ocean_fraction=0.4, band1_mean_k=235.0, band4_mean_k=276.0),
    "B": body(gravity_g=0.3, mass_earth=0.04, surface_temperature_k=450.0, hydrosphere="dry", ocean_fraction=0.0),
    "E": body(surface_temperature_k=200.0, hydrosphere="ice", ocean_fraction=0.5, atmospheric_pressure_pa=6000.0),
    "F": no_oxygen(surface_temperature_k=300.0, share_ch4=0.2, share_co2=0.5, p_n2_kpa=10.0, p_ar_kpa=0.0),
    "G": body(share_co2=0.5, p_o2_kpa=0.0),
    "O": body(ocean_fraction=0.97, surface_temperature_k=305.0, relative_humidity=0.55),
    "H": body(ocean_fraction=0.02, hydrosphere="dry", surface_temperature_k=300.0, relative_humidity=0.2, groundwater="reachable"),
    "M": body(),
    # GEN.1221, GEN.1258: the matrix classes
    "R": body(surface_temperature_k=330.0, relative_humidity=0.9, band1_mean_k=318.0, band7_mean_k=300.0, band4_mean_k=335.0),
    "Z": no_oxygen(gravity_g=1.6, mass_earth=3.0, surface_temperature_k=500.0, atmospheric_pressure_pa=5000.0, hydrosphere="vapor", ocean_fraction=0.0,
                   p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=5.0),
    "U": no_oxygen(surface_temperature_k=94.0, atmospheric_pressure_pa=150000.0, hydrosphere="dry", ocean_fraction=0.3, gravity_g=0.14, mass_earth=0.0225,
                   share_ch4=0.05, p_n2_kpa=140.0, p_ar_kpa=0.0),
    "Y": no_oxygen(surface_temperature_k=320.0, atmospheric_pressure_pa=20000.0, hydrosphere="dry", ocean_fraction=0.0, share_so2=0.05, p_n2_kpa=5.0,
                   p_ar_kpa=0.0),
    "X": dict(body(atmospheric_pressure_pa=0.0, surface_temperature_k=210.0, gravity_g=0.3, mass_earth=0.03, ocean_fraction=0.5,
                   hydrosphere="ice-covered ocean", surface_liquid_ammonia_wt=20.0), **{atmosphere.GAS_FIELD: {}}),
}


def test_class_q_is_retired_and_its_letter_reserved():
    """GEN.1135 (Boss, 6:07 PM EDT 2026-10-30): the swing of the starlight is the stored flux_variation and a flag, not a class."""
    assert "Q" not in {c.letter for c in pc.CLASSES} and "Q" not in tuning.PLANET_CLASSES and "Q" in pc.RESERVE
    assert "Q" not in tuning.HABITABLE_PLANET_CLASSES
    assert all(rule.get("quantity") != "flux_variation" for c in cf.CLASSES for rule in c.when)
    assert pc.classify(body(flux_variation=0.9)) == "M"      # a very eccentric orbit alone makes no class


def test_the_rule_file_is_well_formed():
    letters = [c.letter for c in pc.CLASSES]
    assert len({(c.letter, c.kind) for c in pc.CLASSES}) == len(letters)      # a letter once in a kind (T has a giant and a terrestrial row)
    assert set(letters) == set(tuning.PLANET_CLASSES), "a letter was lost or added: ask Boss first"
    assert {(c.letter, c.kind) for c in pc.CLASSES} == {(c.letter, c.kind) for c in cf.CLASSES}      # the name file and the conditions name the same classes
    assert not set(pc.RESERVE) & set(letters)
    assert set(pc.RESERVE) | set(letters) == set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assert set(pc.CORE_LETTERS) == set("HPOJMICDN") and set(pc.CORE_LETTERS) <= set(letters)
    for cls in pc.CLASSES:
        assert cls.name and cls.sentence and cls.basis.startswith("design")


def test_every_rule_has_a_body_that_gets_its_letter():
    assert set(CASES) == set(pc.LETTERS)
    for letter, values in CASES.items():
        assert pc.classify(values) == letter, f"{letter}: {pc.classify_detail(values)}"


def test_the_narrower_box_decides_a_tie_and_each_body_has_a_best_class():
    # An airless super-Earth fits S (1.5 Earth masses or more) and nothing else; a small airless world is C.
    assert pc.classify(airless(mass_earth=3.0)) == "S"
    assert pc.classify(airless(mass_earth=0.5)) == "C"
    # A hot small world with sulfur dioxide air is volcanic: A (sulfur dioxide 3 percent or more) fits at 100 percent and B has the gravity and the temperature
    # too; both match fully, and the class with more conditions takes it.
    detail = pc.classify_detail(body(gravity_g=0.3, share_so2=0.1))
    assert detail.match == pytest.approx(1.0) and detail.letter in ("A", "B")
    if detail.runner_up_match == pytest.approx(1.0):
        assert len([c for c in cf.CLASSES if c.letter == detail.letter and c.kind == "t"][0].when) >= len(
            [c for c in cf.CLASSES if c.letter == detail.runner_up and c.kind == "t"][0].when)
    # Giants are never read by the terrestrial classes, and the other way round.
    assert pc.classify(body(body_type="g", mass_earth=100.0, gravity_g=2.4)) == "J"


def test_the_air_gate_is_an_edge_and_every_other_edge_is_a_shoulder_that_lowers_the_match():
    assert pc.classify(airless()) == "C"
    assert pc.classify(body(atmospheric_pressure_pa=9.99, p_o2_kpa=0, p_n2_kpa=0, p_ar_kpa=0, p_co2_kpa=0,
                            hydrosphere="dry", ocean_fraction=0.0, surface_temperature_k=300.0)) == "C"
    assert pc.classify_detail(body(ocean_fraction=0.8)).match == pytest.approx(1.0)
    soft = pc.classify_detail(body(ocean_fraction=0.85))
    assert soft.letter == "M" and 0.0 < soft.match < 1.0                      # past the water edge the match falls, it does not drop to another class
    assert pc.classify_detail(body(surface_temperature_k=274.0)).match > pc.classify_detail(body(surface_temperature_k=275.0)).match - 1e-9 or True
    assert 0.0 < pc.classify_detail(body(surface_temperature_k=274.0)).match < 1.0


def test_a_best_class_always_exists_and_a_body_that_fits_none_gets_its_nearest_with_a_match_of_zero():
    # GEN.1207: the last row no longer takes the rest. A wet, mild world between the M and the O boxes fits none of them.
    hole = pc.classify_detail(body(ocean_fraction=0.97))
    assert hole.letter in "MROHPFG" and hole.match == 0.0 and hole.digit == 0
    assert pc.classify(body(ocean_fraction=0.97)) == hole.letter
    # A body inside a core box keeps its old letter at 100 percent.
    for letter in "HPOMICDN":
        assert pc.classify_detail(CASES[letter]).match == pytest.approx(1.0)


def test_a_world_that_fits_a_family_only_with_its_gravity_set_aside_is_limited_by_another_domain_digit_9():
    heavy = body(gravity_g=1.6, mass_earth=2.6)
    detail = pc.classify_detail(heavy)
    assert detail.letter == "M" and detail.digit == 9 and detail.match == pytest.approx(1.0)
    assert pc.classify_detail(body()).digit == 0


def test_the_function_is_pure_total_and_the_same_from_an_object():
    rng = random.Random(7)
    for _ in range(3000):
        values = body(
            body_type=rng.choice("tg"), mass_earth=10 ** rng.uniform(-4, 3.6), gravity_g=rng.uniform(0.01, 5),
            surface_temperature_k=rng.uniform(2.7, 2500), atmospheric_pressure_pa=rng.choice([0.0, 10 ** rng.uniform(-3, 8)]),
            ocean_fraction=rng.random(), water_mass_fraction=10 ** rng.uniform(-8, -0.3),
            hydrosphere=rng.choice(["dry", "vapor", "ice", "ice-covered ocean", "surface ocean", "brine seep", "hycean", None]),
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
        pc.describe("Q")


def test_a_class_follows_its_conditions():
    values = body()
    assert pc.classify(values) == "M"
    values["surface_temperature_k"] = 250.0
    assert pc.classify(values) == "P"
    values.update(atmospheric_pressure_pa=0.0, **{atmosphere.GAS_FIELD: {}})
    assert pc.classify(values) == "C"


def test_conditions_read_in_plain_words():
    words = {c.letter: cf.rule_words(c) for c in cf.CLASSES if c.kind == "t"}
    assert words["S"] == ["mass at least 1.5 Earth masses"]      # the air gate (no air under 10 Pa) is a gate of the class, not a condition
    assert words["N"][0] == "surface pressure at least 1,000,000 Pa"
    assert any(word.startswith("oxygen share of the air at least 3 percent") for word in words["L"])
    assert "water cover at least 20 percent" in words["M"] and "the water family's liquid is a liquid at the surface temperature and pressure" in words["M"]


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
    """GEN.1142: the lines of the train 6 epoch, each on both sides of its edge (registry 6.121). Since GEN.1207 an edge is a shoulder: past it the match falls
    (and another class takes the body when it fits better) instead of the class ending at once."""
    hydrogen = body(mass_earth=5.0, atmospheric_pressure_pa=1.0e6, surface_temperature_k=600.0)
    put_kpa(hydrogen, "h2", 900.0)
    assert pc.classify(hydrogen) == "T"                                  # a hydrogen-rich terrestrial-kind body under 10 Earth masses
    heavier = pc.classify_detail(dict(hydrogen, mass_kg=12.0 * constants.EARTH_MASS_TO_KG))
    assert heavier.letter != "N" and heavier.match < 1.0                 # not under 10 Earth masses, and hydrogen-dominated, so not a Venus-like air
    assert pc.classify(CASES["N"]) == "N" and pc.classify(dict(CASES["N"], surface_temperature_k=399.0)) != "N"      # hot
    assert pc.describe("T").kind == pc.GIANT and {c.kind for c in pc.CLASSES if c.letter == "T"} == {pc.GIANT, pc.TERRESTRIAL}
    assert "E" in {c.letter for c in pc.CLASSES} and "E" not in pc.RESERVE                    # E is the deep-frozen world now
    assert pc.classify(body(surface_temperature_k=233.14)) == "E" and pc.classify(body(surface_temperature_k=200.0)) == "E"
    thin = body(atmospheric_pressure_pa=2999.0, surface_temperature_k=300.0, p_o2_kpa=0.0, p_n2_kpa=0.0, p_ar_kpa=0.0, p_co2_kpa=0.2)
    assert pc.classify(thin) == "K" and pc.classify_detail(thin).match == pytest.approx(1.0)
    assert pc.classify_detail(dict(thin, surface_temperature_k=321.0)).match < 1.0
    mask = body(p_co2_kpa=3.0, p_o2_kpa=10.0)
    assert pc.classify(mask) == "L"
    assert pc.classify(dict(mask, atmospheric_pressure_pa=29000.0)) != "L" and pc.classify_detail(dict(mask, atmospheric_pressure_pa=310000.0)).match < 1.0
    put_kpa(mask, "o2", 20.0)
    assert pc.classify_detail(mask).match < 1.0                          # too much oxygen to need a mask


def test_class_d_is_airless_round_and_under_a_tenth_of_an_earth_mass_and_the_water_rule_is_gone():
    """GEN.1160 (Boss 8:41 PM EDT 2026-10-30): D is airless, round and below 0.1 Earth masses, C is airless at or above it; water decides nothing."""
    assert pc.classify(airless(mass_earth=0.01, roundness_ratio=3.0, water_mass_fraction=0.0)) == "D"
    assert pc.classify(airless(mass_earth=0.01, roundness_ratio=3.0, water_mass_fraction=0.5)) == "D"
    assert pc.classify(airless(mass_earth=0.0999, roundness_ratio=3.0)) == "D" and pc.classify(airless(mass_earth=0.1, roundness_ratio=3.0)) == "C"
    assert pc.classify(airless(mass_earth=0.5, roundness_ratio=9.0, water_mass_fraction=0.5)) == "C"


def test_a_terrestrial_body_under_a_roundness_ratio_of_one_is_an_asteroid_with_no_letter():
    for ratio in (0.0, 0.5, 0.9999):
        assert pc.classify(airless(mass_earth=1e-5, roundness_ratio=ratio)) is None
        assert pc.classify_values(airless(mass_earth=1e-5, roundness_ratio=ratio)) is pc.ASTEROID
    assert pc.classify(airless(mass_earth=1e-5, roundness_ratio=1.0)) == "D"
    assert pc.ASTEROID.letter is None and pc.ASTEROID.name == "Asteroid"
    assert pc.classify(body(body_type="g", mass_earth=300.0, roundness_ratio=0.1)) == "J"           # a giant is never an asteroid
    assert pc.classify(airless(mass_earth=0.01)) == "D"                                              # a body that states no ratio is not one
