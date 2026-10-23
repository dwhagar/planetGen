# tests/test_atmosphere_mixture_draw.py

"""`draw_mixture` (GEN.479): the family each recipe draw calls for, the inputs, and that every body gets a mixture of gas-table gases."""

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import gases, mixture_draw as md, mixture_giant, species_thermo
from planetgen.physics.atmosphere_recipes import INPUTS, RecipeDraw
from planetgen.physics.mixture_base import MixtureError
from planetgen.util import draw
from tests.test_atmosphere_mixture_giant import JUPITER
from tests.test_atmosphere_mixture_rocky import EARTH


class Row:
    def __init__(self, name):
        self.name = name


def drawn(**groups):
    return RecipeDraw(groups={group: Row(name) for group, name in groups.items()})


def test_the_family_follows_the_groups_the_recipe_drew():
    cold = Subject("A", {"ts_k": 300.0, "tp_k": 300.0})
    assert md.family_of(drawn(g="g2", f="f3"), cold) == "brown-dwarf"
    assert md.family_of(drawn(f="f4", h="h2"), cold) == "hot-jupiter"
    assert md.family_of(drawn(f="f1"), cold) == "giant"
    assert md.family_of(drawn(e="e2", f="f2"), cold) == "giant"           # a class of giants stays a giant
    assert md.family_of(drawn(e="e1", f="f2"), cold) == "giant"
    assert md.family_of(drawn(e="e2"), cold) == "envelope"
    for rocky in (drawn(e="e1"), drawn(e="e3"), drawn(e="e4"), drawn(a="a6", b="b1", c="c3")):
        assert md.family_of(rocky, cold) == "rocky"


def test_the_family_of_the_evolved_cold_coma_and_molten_rows():
    cold, hot = Subject("A", {"ts_k": 300.0, "tp_k": 300.0}), Subject("A", {"ts_k": 2000.0, "tp_k": 2000.0})
    for row, family in (("j1", "venus"), ("j2", "mars"), ("j3", "titan"), ("j4", "frost"), ("j5", "abiotic-oxygen")):
        assert md.family_of(drawn(a="a6", j=row), cold) == family
    assert md.family_of(drawn(k="k1"), cold) == "coma"
    assert md.family_of(drawn(l="l2"), hot) == "airless"
    assert md.family_of(drawn(l="l1"), hot) == "lava"
    assert md.family_of(drawn(l="l1"), cold) == "airless"                  # a surface below the melting temperature has no vapour
    assert md.family_of(drawn(a="a6", i="i2"), hot) == "lava"
    assert md.family_of(drawn(a="a6", i="i1"), cold) == "rocky"           # a hot-zone body that is not molten outgasses


def test_every_gas_the_solver_can_make_of_the_budget_elements_is_in_the_gas_table():
    allowed = set(mixture_giant.METALS) | {"H", "He"}
    missing = [name for name in species_thermo.GASES
               if set(species_thermo.SPECIES[name].elements) <= allowed and not gases.has(species_thermo.SPECIES[name].formula)]
    assert missing == []


def test_the_required_inputs_are_stored_quantities_and_what_the_chemistry_reads():
    needed = md.required_inputs("M")
    assert {"ts_k", "pressure_bar", "mass_earth", "age_gy", "snow_ratio", "ocean_present", "long_lived_geology"} <= set(needed)
    assert "radius_km" not in needed
    assert "redox_diw" not in needed and "o2_pal" not in needed          # drawn by the tables; the life stage comes after the mixture
    assert all(name in INPUTS for name in needed)
    giants = md.required_inputs("J")
    assert {"t_layer_k", "radius_km", "mass_jupiter", "grain_um"} <= set(giants)


def test_a_subject_without_a_quantity_a_row_reads_is_refused_by_name():
    subject = Subject("A", {k: v for k, v in EARTH.items() if k != "snow_ratio"})
    with draw.bound(0), pytest.raises(MixtureError, match="snow_ratio"):
        md.draw_mixture(subject)


@pytest.mark.parametrize("body_class, changes, family", [
    ("A", {}, "rocky"),
    ("M", {"ts_k": 600.0, "pressure_bar": 30.0}, "rocky"),
    ("J", {"t_layer_k": 120.0, "teq_k": 120.0}, "giant"),
    ("J", {"t_layer_k": 1400.0, "teq_k": 1400.0, "tday_k": 1500.0, "zone": "h", "insolation": 800.0}, "hot-jupiter"),
    ("brown-dwarf", {"mass_jupiter": 30.0, "mass_earth": 9500.0, "teff_k": 900.0, "t_layer_k": 900.0, "radius_km": 70000.0}, "brown-dwarf"),
])
def test_a_body_gets_a_mixture_of_gas_table_gases_and_the_draw_it_followed(body_class, changes, family):
    subject = Subject(body_class, {**EARTH, **JUPITER, **changes} if body_class != "A" and body_class != "M" else {**EARTH, **changes})
    with draw.bound(3):
        result = md.draw_mixture(subject)
    assert result.family == family
    assert sum(result.mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
    assert all(gases.has(formula) for formula in result.mole_fractions)
    assert result.recipe is not None and "d" in result.recipe.groups or family != "rocky"
    assert result.drawn == {group: row.name for group, row in result.recipe.groups.items()}
    assert (result.redox_diw is None) == (family != "rocky")


def test_the_retention_row_and_the_mixture_come_from_one_draw():
    with draw.bound(8):
        first = md.draw_mixture(Subject("A", EARTH))
    with draw.bound(8):
        second = md.draw_mixture(Subject("A", EARTH))
    assert first.mole_fractions == second.mole_fractions
    assert first.recipe.groups["d"].name == second.recipe.groups["d"].name


def test_a_sweep_from_a_cold_giant_to_an_ultra_hot_one_never_refuses():
    for temperature in (45.0, 110.0, 250.0, 600.0, 1100.0, 1800.0, 2600.0, 3400.0):
        for seed in range(3):
            changes = {"t_layer_k": temperature, "tday_k": temperature, "teq_k": temperature, "teff_k": temperature,
                       "zone": "h" if temperature > 1300 else "c", "insolation": (temperature / 255.0) ** 4}
            with draw.bound(seed):
                result = md.draw_mixture(Subject("J", {**JUPITER, **changes}))
            assert sum(result.mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
            assert result.mole_fractions["H2"] + result.mole_fractions.get("H", 0.0) > 0.5


def test_a_random_body_of_any_rocky_class_gets_a_mixture_from_the_inputs_the_function_names():
    import random
    rng = random.Random(7)
    base = dict(EARTH)
    families = set()
    for index in range(150):
        body_class = rng.choice("ABCDEFGHKLMNOPQSVUYZ")
        mass, temperature = 10 ** rng.uniform(-3, 1.3), 10 ** rng.uniform(1.4, 3.5)
        quantities = dict(base, mass_earth=mass, ts_k=temperature, tp_k=temperature * rng.uniform(1.0, 1.3),
                          pressure_bar=10 ** rng.uniform(-6, 2), age_gy=rng.uniform(0.1, 10), long_lived_geology=rng.random() < 0.5,
                          ocean_present=rng.random() < 0.3, snow_ratio=rng.uniform(0.2, 3), zone=rng.choice("hecr"),
                          host_type=rng.choice(["G", "K", "M", "F", "WD", "PSR"]), is_moon=rng.random() < 0.2,
                          insolation=10 ** rng.uniform(-4, 3), luminosity_solar=10 ** rng.uniform(-3, 1), water_lost=rng.random() < 0.2,
                          runaway_ratio=rng.uniform(0.2, 2), teq_k=temperature * 0.8, xuv_index=10 ** rng.uniform(-1, 2),
                          v_esc_kms=10 ** rng.uniform(0, 1.6), radius_km=10 ** rng.uniform(2.5, 4.3), t_layer_k=temperature * 0.9)
        wanted = set(md.required_inputs(body_class))
        subject = Subject(body_class, {name: value for name, value in quantities.items() if name in wanted})
        with draw.bound(index):
            result = md.draw_mixture(subject)
        families.add(result.family)
        assert sum(result.mole_fractions.values()) == pytest.approx(1.0, abs=1e-12) or result.family == "airless"
        assert all(gases.has(formula) for formula in result.mole_fractions)
    assert {"rocky", "airless", "lava"} <= families
