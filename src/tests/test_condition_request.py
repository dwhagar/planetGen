"""A request that names conditions, not a class letter (GEN.688, generation/condition_request.py)."""

import pytest

from planetgen import tuning
from planetgen.generation import condition_request as cr
from planetgen.physics import class_fit, planet_classes


def _earth(**changes):
    values = dict(body_type="t", mass_earth=1.0, pressure_pa=101325.0, temperature_k=288.0, gravity_g=1.0, share_h2=0.0, share_so2=0.0,
                  o2_partial_pa=21000.0, co2_partial_pa=40.0, share_co2=0.0004, reducing_share=0.0)
    values.update(changes)
    return values


def test_a_letter_is_sugar_for_the_conditions_of_its_class():
    request = cr.from_letter("G")
    cls = next(c for c in class_fit.CLASSES if c.letter == "G")
    (option,) = request.options
    assert option.kind == "t" and option.air is True and option.when == cls.when and request.letter == "G"
    assert request.words() == "class G"


def test_a_letter_two_kinds_share_is_either_of_them():
    request = cr.from_letter("T")
    assert {option.kind for option in request.options} == {"g", "t"}
    assert cr.meets(request, _earth(body_type="g", mass_earth=5.0)) and not cr.meets(request, _earth(body_type="g", mass_earth=50.0))


def test_a_letter_that_is_no_class_is_refused():
    with pytest.raises(cr.RequestError):
        cr.from_letter("?")


def test_conditions_become_the_fit_tables_own_shape():
    request = cr.from_conditions({"kind": "t", "air": True, "temperature_k": [273.15, 323.15], "gravity_g": [0.5, None], "phase_water": "surface ocean"})
    (option,) = request.options
    assert {"quantity": "temperature_k", "min": 273.15} in option.when and {"quantity": "temperature_k", "lt": 323.15} in option.when
    assert {"quantity": "gravity_g", "min": 0.5} in option.when and {"quantity": "gravity_g", "lt": None} not in option.when
    assert {"quantity": "phase_water", "equals": "surface ocean"} in option.when
    assert "temperature_k at least 273.15" in request.words()


@pytest.mark.parametrize("spec", [{"temperature": [1, 2]}, {"kind": "x"}, {"air": "yes"}, {"gravity_g": 1.0}, {"gravity_g": [2, 1]}, {"gravity_g": [None, None]},
                                  {"gravity_g": ["a", 2]}, {"gravity_g": [True, 2]}, [1, 2]])
def test_conditions_that_cannot_be_read_are_refused(spec):
    with pytest.raises(cr.RequestError):
        cr.from_conditions(spec)


def test_a_slot_names_a_class_or_conditions_not_both_and_may_name_neither():
    assert cr.from_spec() is None
    assert cr.from_spec(planet_class="M").letter == "M"
    assert cr.from_spec(conditions={"mass_earth": [0.5, 5.0]}).letter is None
    with pytest.raises(cr.RequestError):
        cr.from_spec(planet_class="M", conditions={"mass_earth": [0.5, 5.0]})


def test_a_finished_body_meets_a_request_by_its_conditions_and_not_its_letter():
    request = cr.from_conditions({"temperature_k": [273.15, 323.15], "pressure_pa": [10000.0, None]})
    assert cr.meets(request, _earth())
    assert not cr.meets(request, _earth(temperature_k=330.0))
    assert not cr.meets(request, _earth(pressure_pa=5000.0))
    assert not cr.meets(cr.from_conditions({"kind": "g", "mass_earth": [35.0, None]}), _earth())


def test_the_edges_are_hard_a_little_outside_does_not_meet():
    request = cr.from_conditions({"gravity_g": [0.75, 2.0]})
    assert cr.meets(request, _earth(gravity_g=0.75)) and not cr.meets(request, _earth(gravity_g=0.74)) and not cr.meets(request, _earth(gravity_g=2.0))


def test_the_air_gate_is_the_class_tables():
    assert cr.meets(cr.from_conditions({"air": True}), _earth())
    assert not cr.meets(cr.from_conditions({"air": False}), _earth())
    assert cr.meets(cr.from_conditions({"air": False}), _earth(pressure_pa=0.0, temperature_k=200.0))


def test_a_body_with_no_kind_or_an_unreadable_one_meets_nothing():
    request = cr.from_conditions({})
    assert not cr.meets(request, type("Belt", (), {"body_type": "a"})())
    assert not cr.meets(request, {"body_type": "t"})


def test_a_request_for_a_warm_but_not_hot_surface_is_placed_in_the_ecosphere_and_no_other():
    assert cr.wants_ecosphere(cr.from_letter("M")) and cr.wants_ecosphere(cr.from_letter("G"))
    assert cr.wants_ecosphere(cr.from_conditions({"temperature_k": [275.0, 330.0]}))
    assert not cr.wants_ecosphere(cr.from_letter("E"))                                 # a deep-frozen world keeps its slot's place
    assert not cr.wants_ecosphere(cr.from_letter("T"))                                 # a letter a giant class also uses keeps it too
    assert not cr.wants_ecosphere(cr.from_conditions({"temperature_k": [273.15, 500.0]}))
    assert not cr.wants_ecosphere(cr.from_conditions({"temperature_k": [200.0, 300.0]}))
    assert not cr.wants_ecosphere(cr.from_conditions({"temperature_k": [300.0, None]}))
    assert not cr.wants_ecosphere(cr.from_conditions({"kind": "g", "temperature_k": [300.0, 350.0]}))
    assert not cr.wants_ecosphere(cr.from_conditions({"mass_earth": [1.0, 5.0]}))
    assert tuning.REQUEST_ECOSPHERE_TEMPERATURE_K == (273.15, 400.0)


def test_every_letter_of_the_table_expands_to_a_request():
    for letter in {c.letter for c in class_fit.CLASSES}:
        request = cr.from_letter(letter)
        assert request.options and all(option.when and option.kind in cr.KINDS for option in request.options)
