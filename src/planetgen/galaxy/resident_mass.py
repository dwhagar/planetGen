# planetgen/galaxy/resident_mass.py

"""
What a sector already holds, and what is left of its mass budget (GEN.1124; docs/design/systems/region-object-store.md section 8.5, stage 0)
============================================================================================================================================

Before anything is drawn for a sector, whatever the sector already holds when it is generated comes off the mass it is calculated to have: rogues that arrived from other sectors, bright
stars placed beforehand, scattered phenomena, the parents of belts and extended objects whose centre lies in it. What is left is the budget the phases of the sector build (the mass-backfill
rings, the system draw, phenomena, remnants, belts) draw against. The residents themselves are never moved or changed: this module reads them and returns a new, separate `Budget`.

`remaining_budget(total_kg, residents)` is the pure function:

* A resident counts only when its centre is in the sector (`Resident.centre_inside`); a large object that reaches in from a neighbour is that neighbour's, and is listed as `ignored`.
* The remaining budget is the total less the counted residents, **never negative**; mass the residents hold beyond the total is reported as `over_kg` and the budget is zero (nothing more is drawn).
* A sector without residents returns the total unchanged.
* Each resident may also stand for a number of systems (`Resident.systems`), taken off the expected number of systems the same way (`remaining_systems`).
* The held mass is also split by kind (`held_by_kind`), so each phase can report what it was given.

`shortfall_kg(budget, drawn_kg)` is the stage 3 measure: the budget less everything the sector's build drew, never negative, which the rogue top-up (`physics/rogue_topup`) makes up.

Nothing calls this module yet; the generation phases and the region run wire it (GEN.1120).
"""

import math
from dataclasses import dataclass, field

KINDS = ("rogue", "star", "phenomenon", "belt parent", "extended")
"""tuple: The kinds of resident a sector can hold, for the split of the held mass."""


class ResidentMassError(ValueError):
    """An input the budget cannot be worked out from."""


@dataclass(frozen=True)
class Resident:
    """
    One thing a sector already holds.

    Attributes:
        kind (str): One of `KINDS`.
        mass_kg (float): Its mass, kg, not negative.
        centre_inside (bool): Whether its centre lies in the sector. An extended object reaching in from a neighbour has `False` and is not counted.
        systems (float): How many of the sector's expected systems it stands for, not negative.
    """
    kind: str
    mass_kg: float
    centre_inside: bool = True
    systems: float = 0.0

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ResidentMassError(f"a resident is one of {KINDS}, not {self.kind!r}")
        for name in ("mass_kg", "systems"):
            value = getattr(self, name)
            if not math.isfinite(value) or value < 0.0:
                raise ResidentMassError(f"a resident's {name} must be a finite number, not negative: {value!r}")


@dataclass(frozen=True)
class Budget:
    """
    A sector's budget once its residents are off.

    Attributes:
        total_kg (float): The mass the sector was calculated to hold.
        held_kg (float): The mass of the residents counted.
        remaining_kg (float): `total_kg - held_kg`, never negative.
        over_kg (float): How far the residents pass the total (zero unless they do).
        total_systems (float): The expected systems.
        remaining_systems (float): The expected systems less those the residents stand for, never negative.
        held_by_kind (dict): The counted mass for each kind that has any.
        counted (int): The residents counted.
        ignored (int): The residents left out because their centre lies elsewhere.
    """
    total_kg: float
    held_kg: float
    remaining_kg: float
    over_kg: float
    total_systems: float
    remaining_systems: float
    held_by_kind: dict = field(default_factory=dict)
    counted: int = 0
    ignored: int = 0


def remaining_budget(total_kg, residents=(), total_systems=0.0):
    """
    The mass (and systems) a sector still has to draw once its residents are off.

    Args:
        total_kg (float): The mass the sector is calculated to hold, kg, not negative.
        residents (iterable): The `Resident`s the sector already holds, none changed.
        total_systems (float): The expected systems, not negative.

    Returns:
        Budget: The remaining budget, never negative.

    Raises:
        ResidentMassError: If the total or the systems are negative or not finite, or an entry is not a `Resident`.
    """
    for name, value in (("total_kg", total_kg), ("total_systems", total_systems)):
        if not math.isfinite(value) or value < 0.0:
            raise ResidentMassError(f"{name} must be a finite number, not negative: {value!r}")
    held = 0.0
    systems = 0.0
    by_kind = {}
    counted = ignored = 0
    for resident in residents:
        if not isinstance(resident, Resident):
            raise ResidentMassError(f"a resident is a Resident, not {resident!r}")
        if not resident.centre_inside:
            ignored += 1
            continue
        counted += 1
        held += resident.mass_kg
        systems += resident.systems
        if resident.mass_kg > 0.0:
            by_kind[resident.kind] = by_kind.get(resident.kind, 0.0) + resident.mass_kg
    remaining = max(0.0, total_kg - held)
    return Budget(total_kg=total_kg, held_kg=held, remaining_kg=remaining, over_kg=max(0.0, held - total_kg),
                  total_systems=total_systems, remaining_systems=max(0.0, total_systems - systems),
                  held_by_kind=by_kind, counted=counted, ignored=ignored)


def shortfall_kg(budget, drawn_kg):
    """
    How much of a sector's budget is still unmatched once everything its build drew is counted (stage 3).

    Args:
        budget (Budget): From `remaining_budget`.
        drawn_kg (float): The mass the build drew in every phase, kg, not negative.

    Returns:
        float: `budget.remaining_kg - drawn_kg`, never negative.
    """
    if not math.isfinite(drawn_kg) or drawn_kg < 0.0:
        raise ResidentMassError(f"drawn_kg must be a finite number, not negative: {drawn_kg!r}")
    return max(0.0, budget.remaining_kg - drawn_kg)
