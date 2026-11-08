# planetgen/generation/condition_request.py

"""
A request that names conditions, not a class letter (GEN.688)
===============================================================

A slot of `SystemConfig.SLOTS` and the admin and web forms ask for a planet by what it is like: its kind (rocky or giant), whether it holds an air, and a
range for any of the quantities the classes read (`physics/data/planet_class_fit.json`: `mass_earth`, `gravity_g`, `temperature_k`, `pressure_pa`, the
partial pressures of oxygen and carbon dioxide, the water's phase, the ocean fraction and so on). That is the shape of `habitable_request.meets`, which
asks the same of a habitable world; here the caller says which conditions.

A class letter is only sugar. `from_letter("M")` reads the conditions of the class from the fit table, the same list the classifier grades, and the request
holds those conditions. The planet that is built is never selected, steered or seeded by the letter: a finished body meets a request when every one of
its conditions holds, whatever letter the classifier gives it (it may give another, as a body at the edge of two classes).

A request is a `Request`; `meets(request, body)` judges a finished planet or moon. A request that asks for a surface warm enough for liquid water and
cool enough to keep it can only be met in the ecosphere, and a slot with one is placed there (`wants_ecosphere`).
"""

import dataclasses
from collections.abc import Mapping

from planetgen import tuning
from planetgen.physics import class_fit, planet_classes

KINDS = (planet_classes.TERRESTRIAL, planet_classes.GIANT)
"""tuple: The two kinds of a body a request can name, `"t"` (rocky) and `"g"` (giant)."""

EQUALS_ONLY = ("phase_water", "groundwater_present", "phase_methane_ethane", "phase_sulfuric_acid")
"""tuple: The quantities of the fit table that are a state and not a number: a request names the value it must have."""

QUANTITIES = tuple(sorted({condition["quantity"] for cls in class_fit.CLASSES for condition in cls.when}))
"""tuple: The quantities a request may condition, those the fit table's classes read."""


class RequestError(ValueError):
    """A request that cannot be read."""


@dataclasses.dataclass(frozen=True)
class Option:
    """One way to meet a request: a `kind`, the `air` gate (`True` for an air of 10 Pa or more, `False` for none, `None` for either) and each condition in the fit table's own shape."""

    kind: str
    air: object
    when: tuple


@dataclasses.dataclass(frozen=True)
class Request:
    """The conditions a planet must meet: one `Option`, or several when a class letter is shared by a rocky and a giant class (T), which a body meets by meeting any one."""

    options: tuple
    letter: object = None
    """str or None: The class letter this request is the sugar of, when it came from one (it is only a name for the report; nothing judges by it)."""

    def words(self):
        """The request as a phrase, for the system's report of an unmet request."""
        if self.letter:
            return f"class {self.letter}"
        return " or ".join(
            ("these conditions: " + ", ".join(_condition_words(c) for c in option.when)) if option.when else f"a body of kind {option.kind}" for option in self.options)


def _condition_words(condition):
    quantity = condition["quantity"]
    if "equals" in condition:
        return f"{quantity} is {condition['equals']}"
    if "min" in condition:
        return f"{quantity} at least {condition['min']:g}"
    return f"{quantity} under {condition['lt']:g}"


def from_letter(letter):
    """
    The request a class letter stands for: the kind, air gate and conditions of that class in the fit table (a letter two kinds share, as T, gives both: a body meets
    it by meeting either).

    Raises:
        RequestError: The letter is no class.
    """
    rows = [cls for cls in class_fit.CLASSES if cls.letter == letter]
    if not rows:
        raise RequestError(f"{letter!r} is not a planet class.")
    return Request(tuple(Option(cls.kind, cls.air, cls.when) for cls in rows), letter)


def from_conditions(spec):
    """
    The request a mapping of conditions names. Keys: `"kind"` (`"t"` or `"g"`, default `"t"`), `"air"` (`true`, `false` or leave out), and any quantity of `QUANTITIES`, to
    `[low, high]` (low reached, high not; `null` leaves an end open) or, for a state (`EQUALS_ONLY`), the value it must have.

    Raises:
        RequestError: A key that is no quantity, a range that is not two numbers, or a state without a value.
    """
    if not isinstance(spec, dict):
        raise RequestError("the conditions must be an object of quantities.")
    kind = spec.get("kind", planet_classes.TERRESTRIAL)
    if kind not in KINDS:
        raise RequestError(f"kind must be 't' or 'g', not {kind!r}.")
    air = spec.get("air")
    if air is not None and not isinstance(air, bool):
        raise RequestError("air must be true, false or left out.")
    when = []
    for quantity, value in spec.items():
        if quantity in ("kind", "air"):
            continue
        if quantity not in QUANTITIES:
            raise RequestError(f"{quantity!r} is not a quantity a request can condition; the quantities are {', '.join(QUANTITIES)}.")
        if quantity in EQUALS_ONLY:
            when.append({"quantity": quantity, "equals": value})
            continue
        if not (isinstance(value, (list, tuple)) and len(value) == 2):
            raise RequestError(f"{quantity} must be a range [low, high] (null for an open end).")
        low, high = value
        for end in (low, high):
            if end is not None and (isinstance(end, bool) or not isinstance(end, (int, float))):
                raise RequestError(f"{quantity} has an end that is not a number: {end!r}.")
        if low is None and high is None:
            raise RequestError(f"{quantity} names no end.")
        if low is not None and high is not None and not low < high:
            raise RequestError(f"{quantity} must have its low end under its high end.")
        if low is not None:
            when.append({"quantity": quantity, "min": float(low)})
        if high is not None:
            when.append({"quantity": quantity, "lt": float(high)})
    return Request((Option(kind, air, tuple(when)),))


def from_spec(planet_class=None, conditions=None):
    """The request of a slot or a form: a class letter (sugar) or conditions, never both; `None` when it names neither."""
    if planet_class is not None and conditions is not None:
        raise RequestError("name a class or conditions, not both.")
    if planet_class is not None:
        return from_letter(planet_class)
    if conditions is not None:
        return from_conditions(conditions)
    return None


def _holds(condition, found):
    value = found.get(condition["quantity"])
    if value is None:
        return False
    if "equals" in condition:
        return value == condition["equals"]
    if "min" in condition:
        return value >= condition["min"]
    return value < condition["lt"]


def meets(request, body):
    """Whether the finished planet or moon `body` meets every condition of one of the options of `request` (a low end is reached and a high end is not: no shoulder, a body outside by a little does not meet it)."""
    body_type = body.get("body_type") if isinstance(body, Mapping) else getattr(body, "body_type", None)
    if body_type not in KINDS:
        return False
    try:
        found = class_fit.found_quantities(body)
    except planet_classes.ClassError:
        return False
    return any(_option_met(option, found) for option in request.options)


def _option_met(option, found):
    if found["body_type"] != option.kind:
        return False
    if option.air is not None and option.kind == planet_classes.TERRESTRIAL and (found["pressure_pa"] >= class_fit.GATE_PRESSURE_PA) != option.air:
        return False
    return all(_holds(condition, found) for condition in option.when)


def wants_ecosphere(request):
    """
    Whether only a world in the ecosphere can meet `request`: every option is a rocky body that must have a mean surface temperature of at least the freezing point and under
    `tuning.REQUEST_ECOSPHERE_TEMPERATURE_K`'s high end. A request with no such temperature (a cold world, a hot one, a giant, a letter a giant also uses) keeps the place its slot has.
    """
    low_floor, high_ceiling = tuning.REQUEST_ECOSPHERE_TEMPERATURE_K

    def temperate(option):
        lows = [c["min"] for c in option.when if c["quantity"] == "temperature_k" and "min" in c]
        highs = [c["lt"] for c in option.when if c["quantity"] == "temperature_k" and "lt" in c]
        return option.kind == planet_classes.TERRESTRIAL and bool(lows) and bool(highs) and max(lows) >= low_floor and min(highs) <= high_ceiling

    return all(temperate(option) for option in request.options)
