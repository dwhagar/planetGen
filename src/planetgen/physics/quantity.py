# planetgen/physics/quantity.py

"""
The Quantity value class (GEN.569; docs/design/architecture/unified-unit-quantity.md, section 4.3)
=====================================================================================

A number in one unit of the registry (`unit_registry`), kept as given and answered in every other unit of its kind. It is
for the edges, where a value is stored, loaded, shown or handed to the API; the generator's kernels stay on plain floats
in a field's declared unit (Q2), so nothing here is on a hot path.

* **The truth is never converted.** `Quantity("power", 3.828e26, "W")` holds 3.828e26 and "W" for good until `set`
  changes them. Every other form is the truth scaled by the two factors, never one derived form from another, so no
  drift builds up, and a conversion between units whose ratio is a power of ten is one correctly rounded step
  (`unit_registry.scale`): 1.5 km is exactly 1500.0 m.
* **Forms are lazy.** `quantity.kW`, `quantity.Lsun` and `quantity.as_("erg_s")` are worked out on first read, kept, and
  dropped by `set`, so every form is always in step with the truth.
* **Relations** (`unit_registry.RELATIONS`) answer a quantity of another kind from this one: a mass's gravitational
  parameter, a radius's diameter, a period's frequency (`related`), always in the partner's base unit.
* **Checks.** A value is finite and inside its kind's bounds (`unit_registry.check`) when a quantity is made or set; a
  unit of another kind is refused with a message naming both.
* **A quantity compares equal** to another of the same kind with the same truth in the same unit, or in another unit
  with the same value in base units, within a relative 1e-12.

Dimensionless values (fractions, ratios, counts) are not quantities.
"""

import math

from planetgen.physics import unit_registry as units

_UNKNOWN_ATTRIBUTE = "'Quantity' object has no attribute {!r}"


class Quantity:
    """A number of one kind in one unit, and every other unit of that kind on demand."""

    __slots__ = ("kind", "value", "unit", "_forms")

    def __init__(self, kind, value, unit):
        """
        Args:
            kind (str): A kind of `unit_registry.KINDS`.
            value (float): The number, in `unit`; finite and inside the kind's bounds.
            unit (str): A key of `unit_registry.UNITS` whose kind is `kind`.

        Raises:
            ValueError: For a kind or unit the registry does not hold, a unit of another kind, or a value outside the
                kind's bounds.
        """
        units.kind(kind)
        row = units.unit(unit)
        if row.kind != kind:
            raise ValueError(f"{row.symbol} ({unit}) measures {row.kind}, not {kind}")
        self.kind = kind
        self.unit = unit
        self.value = units.check(float(value), unit)
        self._forms = None

    # --- Forms ---------------------------------------------------------------------------------------------------

    def as_(self, unit):
        """This quantity's number in `unit` (a key of this kind): the truth itself for its own unit, otherwise the truth
        scaled by the two units' factors, kept after the first read."""
        if unit == self.unit:
            return self.value
        forms = self._forms
        if forms is None:
            forms = self._forms = {}
        elif unit in forms:
            return forms[unit]
        form = forms[unit] = units.convert(self.value, self.unit, unit)
        return form

    def __getattr__(self, key):
        # Reached only for a name that is not a slot: a unit key of this quantity's kind (`quantity.kW`); a key that is no
        # Python name ("m-3") is read with `as_`. A unit of another kind is an AttributeError (so `hasattr` is False).
        if key.startswith("_") or key not in units.UNITS:
            raise AttributeError(_UNKNOWN_ATTRIBUTE.format(key))
        row = units.UNITS[key]
        if row.kind != self.kind:
            raise AttributeError(f"a {self.kind} has no {row.symbol} ({key}), which measures {row.kind}")
        return self.as_(key)

    def set(self, value, unit=None):
        """Changes the truth to `value` of `unit` (this quantity's own when omitted); every form follows."""
        unit = self.unit if unit is None else unit
        row = units.unit(unit)
        if row.kind != self.kind:
            raise ValueError(f"{row.symbol} ({unit}) measures {row.kind}, not {self.kind}")
        self.value = units.check(float(value), unit)
        self.unit = unit
        self._forms = None

    def base(self):
        """The quantity in its kind's base unit, as a number (offsets applied)."""
        return units.to_base(self.value, self.unit)

    # --- Relations ---------------------------------------------------------------------------------------------------

    def related(self, relation, role=None):
        """
        The partner of this quantity in a declared relation (`unit_registry.RELATIONS`), as a new `Quantity` in the
        partner kind's base unit.

        Args:
            relation (str): "mass-mu", "radius-diameter", "period-frequency", "rotation_period-angular_rate".
            role (str, optional): Which side of the relation this quantity is ("radius" or "diameter", for two of one
                kind); the side its kind alone names when omitted.

        Raises:
            ValueError: For an unknown relation, or a quantity that is neither side, or whose side is ambiguous.
        """
        declared = units.RELATIONS.get(relation)
        if declared is None:
            raise ValueError(f"unknown relation {relation!r} ({', '.join(units.RELATIONS)})")
        if role is None:
            sides = [(side, kind) for side, kind in ((declared.first_role, declared.first_kind),
                                                      (declared.second_role, declared.second_kind))
                     if kind == self.kind]
            if len(sides) != 1:
                raise ValueError(f"{relation} needs a role: this {self.kind} could be "
                                 f"{' or '.join(side for side, _kind in sides) or 'neither side'}")
            role = sides[0][0]
        if role == declared.first_role:
            own_kind, partner_kind = declared.first_kind, declared.second_kind
        elif role == declared.second_role:
            own_kind, partner_kind = declared.second_kind, declared.first_kind
        else:
            raise ValueError(f"{relation} relates {declared.first_role} and {declared.second_role}, not {role!r}")
        if own_kind != self.kind:
            raise ValueError(f"{relation}'s {role} is a {own_kind}, not a {self.kind}")
        return Quantity(partner_kind, units.apply_relation(relation, self.base(), role), units.KINDS[partner_kind].base)

    # --- Dict, equality, text ------------------------------------------------------------------------------------------

    def to_dict(self):
        """The describe and API shape: `{"kind", "value", "unit"}`, the truth as given."""
        return {"kind": self.kind, "value": self.value, "unit": self.unit}

    @classmethod
    def from_dict(cls, data):
        """A quantity from `to_dict`'s shape; a missing or unknown key is a `ValueError`."""
        try:
            return cls(data["kind"], data["value"], data["unit"])
        except (KeyError, TypeError) as error:
            raise ValueError(f"not a quantity: {data!r}") from error

    def __eq__(self, other):
        if not isinstance(other, Quantity):
            return NotImplemented
        if self.kind != other.kind:
            return False
        if self.unit == other.unit:
            return self.value == other.value
        return math.isclose(self.base(), other.base(), rel_tol=1e-12, abs_tol=0.0)

    __hash__ = None  # mutable (`set`), and equal across units

    def __repr__(self):
        return f"Quantity({self.kind!r}, {self.value!r}, {self.unit!r})"

    def __str__(self):
        return f"{self.value:g} {units.unit(self.unit).symbol}"

    def __copy__(self):
        return Quantity(self.kind, self.value, self.unit)

    def __reduce__(self):
        return (Quantity, (self.kind, self.value, self.unit))
