# planetgen/physics/elements.py

"""
The elements and their isotopes (GEN.438, docs/design/composition/composition-model.md section 7.1 step 1)
================================================================================================

The one table of atomic data the gas, materials and stellar tables take their molar masses from, and the one set of
isotope constants the radiogenic heat and the ground's dose read. Atomic physics only: no abundance pattern lives here,
no weight, no draw; the weights that use it are the abundance rows of GEN.439 and the recipes above it. It is a constants
table and not a possibility table (`generation/possibility_table.py`): it has no weights and no class gate.

The data are three files in `physics/data/`, built by `scripts/build_element_tables.py` from published tables and checked in
with the values read from them; each states its sources, citations and retrieval date, and the loader below
refuses a file that does not (`DataError`):

* `elements.json`: atomic number, symbol, name and atomic weight of the 118 elements, from CIAAW's Atomic Weights 2021
  (with its 2024 revisions). `atomic_weight` is the abridged conventional value, the one an element with an interval
  standard weight (H, Li, B, C, N, O, Mg, Si, S, Cl, Ar, Br, Tl) also has; the full standard weight or interval is kept.
* `isotopes.json`: every isotope of natural composition (NIST, from IUPAC and the atomic mass evaluation) and the
  radioactive isotopes the heat, the dose and the short-lived-radionuclide heating need, with the ground-state
  half-life, the decay modes and their branchings (NUBASE2020), and the energy released per decay ([R2768: Table 4]), with [R1203]'s heat constants of the four long-lived parents beside it as a cross-check.
* `decay_series.json`: the decay chains of uranium-238, uranium-235 and thorium-232 down every branch of 0.005 percent or
  more to the stable lead isotope, each step with the energy computed from the mass excesses.

Units: molar masses in kg/mol (the unit the rest of the generator stores, `habitability_world.MOLAR_KG`; the files hold
g/mol, the unit the tables are published in), half-lives in seconds, energies in joules, the atomic mass in unified atomic
mass units. Where two sources disagree the newer verified one is used (Boss, 2026-10-17: when in doubt verified wins, when
two verified sources are in doubt the newest wins); the half-lives are NUBASE2020's (2021), Dye's decay constants (2012) are
kept only as the cross-check the tests use.

Use:

    molar_mass("CaSO4.2H2O")            # kg/mol, from the formula
    element("Fe").atomic_weight         # g/mol
    isotope("U-238").half_life_s
    heat_per_kg_of_isotope("K-40")      # W per kg of the isotope itself (all of its decays, neutrinos left out)
"""

import json
import math
import os
import re
from dataclasses import dataclass

AVOGADRO = 6.02214076e23
"""float: Avogadro's constant, 1/mol (exact in the SI since 2019)."""

DATA_DIRECTORY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
"""str: Where the data files are."""

KG_PER_G = 1e-3


class DataError(ValueError):
    """A data file or a row of it that the loader refuses at load."""


class UnknownElement(KeyError):
    """A symbol, name or atomic number that is not an element."""


@dataclass(frozen=True)
class Element:
    """
    One element.

    Attributes:
        z (int): Atomic number.
        symbol (str): The symbol (`"Fe"`).
        name (str): The name in English (`"Iron"`).
        atomic_weight (float or None): The abridged conventional standard atomic weight, g/mol; `None` for an element that
            has no standard weight (every isotope is radioactive and none is common enough to define one).
        atomic_weight_uncertainty (float or None): The stated uncertainty of `atomic_weight`.
        standard_weight (float or None): The full-precision standard weight where the element has one.
        standard_interval (tuple or None): The standard weight as an interval, for an element whose natural variation
            exceeds the measurement uncertainty.
        has_natural_isotopes (bool): Whether any isotope occurs in nature (the primordial ones and the daughters of their
            chains).
        longest_lived_mass_number (int or None): For an element with no standard weight, the mass number of its
            longest-lived isotope (NUBASE2020).
    """

    z: int
    symbol: str
    name: str
    atomic_weight: float | None
    atomic_weight_uncertainty: float | None
    standard_weight: float | None
    standard_interval: tuple | None
    has_natural_isotopes: bool
    longest_lived_mass_number: int | None

    @property
    def molar_mass(self):
        """The molar mass in kg/mol.

        Raises:
            DataError: The element has no standard atomic weight.
        """
        if self.atomic_weight is None:
            raise DataError(f"{self.name} has no standard atomic weight; use an isotope's mass")
        return self.atomic_weight * KG_PER_G


@dataclass(frozen=True)
class Isotope:
    """
    One isotope.

    Attributes:
        z (int): Atomic number.
        a (int): Mass number.
        name (str): `"U-238"`.
        mass_u (float): The mass of the neutral atom in unified atomic mass units (NIST, from the AME).
        mass_uncertainty_u (float or None): Its uncertainty.
        abundance (float or None): The atom fraction of its element in normal material; `None` where it is not found in
            nature in a standard proportion.
        stable (bool): Whether it is observationally stable.
        half_life_s (float or None): The ground-state half-life in seconds; `None` for a stable or unmeasured one.
        half_life_uncertainty_s (float or None): Its uncertainty.
        modes (tuple): `(mode, percent)` pairs from NUBASE2020 (`"A"` alpha, `"B-"`, `"B+"` beta plus or electron
            capture, `"SF"`, `"EC"`, ...); a branch NUBASE2020 gives as approximately 100 percent has the remainder
            after the others and a third item `"~"`.
        decay_energy_j (float or None): The energy released per decay, joules, without the neutrinos [R2768];
            for U-238, U-235 and Th-232 the whole chain to lead.
        dye_2012 (dict or None): Dye's heat constants where he gives them (U-238, U-235, Th-232, K-40).
        mass_estimated (bool): The mass is NIST's estimate from systematic trends, not a measurement.
    """

    z: int
    a: int
    name: str
    mass_u: float
    mass_uncertainty_u: float | None
    abundance: float | None
    stable: bool
    half_life_s: float | None
    half_life_uncertainty_s: float | None
    modes: tuple
    decay_energy_j: float | None
    dye_2012: dict | None
    mass_estimated: bool = False

    @property
    def decay_constant(self):
        """The decay constant in 1/s.

        Raises:
            DataError: The isotope has no half-life (it is stable or none is measured).
        """
        if not self.half_life_s:
            raise DataError(f"{self.name} has no half-life")
        return math.log(2.0) / self.half_life_s

    @property
    def molar_mass(self):
        """The molar mass in kg/mol: the atomic mass in u times 1 g/mol (the molar mass constant, within 1e-9)."""
        return self.mass_u * KG_PER_G


@dataclass(frozen=True)
class DecayStep:
    """One step of a decay chain: `parent` decays by `mode` (`"alpha"`, `"beta-"`, `"ec"`) with `percent` of its decays to
    `daughter`, releasing `q_kev` (the Q value from the mass excesses; `None` where one mass is missing)."""

    parent: str
    daughter: str
    mode: str
    percent: float
    q_kev: float | None


@dataclass(frozen=True)
class DecaySeries:
    """The chain of a long-lived parent to the stable end: `members` (heaviest first) and every branch as `steps`."""

    parent: str
    end: str
    members: tuple
    steps: tuple

    def daughters(self, name):
        """The steps that leave nuclide `name`."""
        return tuple(step for step in self.steps if step.parent == name)


def _load(filename, rows_key):
    path = os.path.join(DATA_DIRECTORY, filename)
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    for field in ("sources", "retrieved", "title"):
        if not data.get(field):
            raise DataError(f"{filename} states no {field}")
    if not data.get(rows_key):
        raise DataError(f"{filename} holds no {rows_key}")
    return data


def _load_elements():
    rows = []
    for item in _load("elements.json", "elements")["elements"]:
        try:
            interval = item.get("standard_interval")
            rows.append(Element(
                z=int(item["z"]), symbol=item["symbol"], name=item["name"], atomic_weight=item.get("atomic_weight"),
                atomic_weight_uncertainty=item.get("atomic_weight_uncertainty"),
                standard_weight=item.get("standard_weight"), standard_interval=None if interval is None else tuple(interval),
                has_natural_isotopes=bool(item["has_natural_isotopes"]),
                longest_lived_mass_number=item.get("longest_lived_mass_number")))
        except KeyError as exc:
            raise DataError(f"elements.json: an element lacks {exc}") from None
    if [row.z for row in rows] != list(range(1, len(rows) + 1)):
        raise DataError("elements.json: the atomic numbers are not 1 to N in order")
    return tuple(rows)


def _load_isotopes():
    rows = []
    for item in _load("isotopes.json", "isotopes")["isotopes"]:
        try:
            rows.append(Isotope(
                z=int(item["z"]), a=int(item["a"]), name=item["name"], mass_u=float(item["mass_u"]),
                mass_uncertainty_u=item.get("mass_uncertainty_u"), abundance=item.get("abundance"),
                stable=bool(item["stable"]), half_life_s=item.get("half_life_s"),
                half_life_uncertainty_s=item.get("half_life_uncertainty_s"),
                modes=tuple(tuple(mode) for mode in item["modes"]),
                decay_energy_j=item.get("decay_energy_j"), dye_2012=item.get("dye_2012"),
                mass_estimated=bool(item.get("mass_estimated"))))
        except KeyError as exc:
            raise DataError(f"isotopes.json: an isotope lacks {exc}") from None
    return tuple(rows)


def _load_series():
    rows = []
    for item in _load("decay_series.json", "series")["series"]:
        steps = tuple(DecayStep(step["from"], step["to"], step["mode"], step["percent"], step["q_kev"])
                      for step in item["steps"])
        rows.append(DecaySeries(item["parent"], item["end"], tuple(item["members"]), steps))
    return tuple(rows)


ELEMENTS = _load_elements()
"""tuple[Element]: Every element, by atomic number."""
ISOTOPES = _load_isotopes()
"""tuple[Isotope]: The isotopes of natural composition and the radioactive ones the generator needs, by atomic number and
mass number."""
SERIES = _load_series()
"""tuple[DecaySeries]: The chains of uranium-238, uranium-235 and thorium-232."""

_BY_SYMBOL = {element.symbol: element for element in ELEMENTS}
_BY_NAME = {element.name.lower(): element for element in ELEMENTS}
_ISOTOPE_BY_NAME = {isotope.name: isotope for isotope in ISOTOPES}


def element(key):
    """The element with symbol (`"Fe"`), name (`"iron"`, any case) or atomic number `key`.

    Raises:
        UnknownElement: There is none.
    """
    try:
        if isinstance(key, int):
            if not 1 <= key <= len(ELEMENTS):
                raise KeyError(key)
            return ELEMENTS[key - 1]
        return _BY_SYMBOL[key] if key in _BY_SYMBOL else _BY_NAME[key.lower()]
    except KeyError:
        raise UnknownElement(f"no element {key!r}") from None


def isotope(name, mass_number=None):
    """The isotope `"U-238"`, or `isotope("U", 238)`.

    Raises:
        KeyError: The table holds no such isotope.
    """
    if mass_number is not None:
        name = f"{name}-{mass_number}"
    return _ISOTOPE_BY_NAME[name]


def isotopes_of(symbol):
    """The table's isotopes of an element, by mass number."""
    z = element(symbol).z
    return tuple(isotope for isotope in ISOTOPES if isotope.z == z)


def series(parent):
    """The decay chain whose parent is `"U-238"`, `"U-235"` or `"Th-232"`.

    Raises:
        KeyError: There is none.
    """
    for found in SERIES:
        if found.parent == parent:
            return found
    raise KeyError(f"no decay chain for {parent!r}")


_ATOM = re.compile(r"([A-Z][a-z]?)(\d*)|([(\[])|([)\]])(\d*)")


def parse_formula(formula):
    """
    The atoms of a chemical formula as `{symbol: count}`: element symbols with whole-number subscripts, groups in
    parentheses or brackets with a subscript (`"Ca3(PO4)2"`, `"Fe[OH]3"`), and hydrates after a dot (`"CaSO4.2H2O"`, the
    centre dot too).

    Raises:
        ValueError: The text is not a formula of whole numbers of known elements.
    """
    text = formula.strip().replace("·", ".").replace("•", ".")
    if not text:
        raise ValueError("an empty formula")
    total = {}
    for part in text.split("."):
        multiplier = re.match(r"(\d+)", part)
        count = int(multiplier.group(1)) if multiplier else 1
        for symbol, number in _parse_part(part[multiplier.end():] if multiplier else part, formula).items():
            total[symbol] = total.get(symbol, 0) + number * count
    return total


def _parse_part(part, formula):
    position, stack = 0, [{}]
    closers = {")": "(", "]": "["}
    openers = []
    while position < len(part):
        match = _ATOM.match(part, position)
        if not match:
            raise ValueError(f"cannot read {formula!r} at {part[position:]!r}")
        position = match.end()
        if match.group(1):
            if match.group(1) not in _BY_SYMBOL:
                raise ValueError(f"{formula!r}: {match.group(1)!r} is not an element")
            number = int(match.group(2)) if match.group(2) else 1
            if number == 0:
                raise ValueError(f"{formula!r}: a subscript of 0")
            stack[-1][match.group(1)] = stack[-1].get(match.group(1), 0) + number
        elif match.group(3):
            openers.append(match.group(3))
            stack.append({})
        else:
            if not openers or openers.pop() != closers[match.group(4)]:
                raise ValueError(f"{formula!r}: unbalanced brackets")
            inner = stack.pop()
            number = int(match.group(5)) if match.group(5) else 1
            if number == 0:
                raise ValueError(f"{formula!r}: a subscript of 0")
            for symbol, count in inner.items():
                stack[-1][symbol] = stack[-1].get(symbol, 0) + count * number
    if openers:
        raise ValueError(f"{formula!r}: unbalanced brackets")
    return stack[0]


def molar_mass(formula):
    """The molar mass of a compound in kg/mol: the sum of the abridged standard atomic weights of its atoms (`parse_formula`).

    Raises:
        ValueError: The formula cannot be read.
        DataError: An element of it has no standard atomic weight.
    """
    return math.fsum(_BY_SYMBOL[symbol].molar_mass * count for symbol, count in parse_formula(formula).items())


def heat_per_kg_of_isotope(name):
    """
    The power of one kilogram of the isotope itself, W/kg: Avogadro's constant times the decay constant times the energy
    released per decay (`Isotope.decay_energy_j`, the neutrinos left out, the whole chain for a parent that has one) over
    the molar mass ([R1203], equation 6). Multiply by the isotope's mass in the rock; the radiogenic heat module
    (GEN.459) does that with the abundances it draws.

    Raises:
        DataError: The isotope has no half-life or no decay energy.
    """
    found = isotope(name)
    if found.decay_energy_j is None:
        raise DataError(f"{name} has no decay energy in the table")
    return AVOGADRO * found.decay_constant * found.decay_energy_j / found.molar_mass


def specific_activity(name):
    """The activity of one kilogram of the isotope itself, Bq/kg: Avogadro's constant times the decay constant over the
    molar mass."""
    found = isotope(name)
    return AVOGADRO * found.decay_constant / found.molar_mass


def isotopic_mass_fraction(name):
    """The mass fraction of the isotope in its element in normal material: abundance times its mass over the element's
    mean isotopic mass (the sum over the element's natural isotopes in the table).

    Raises:
        DataError: The isotope has no natural abundance.
    """
    found = isotope(name)
    if found.abundance is None:
        raise DataError(f"{name} has no natural abundance")
    mean = math.fsum(item.abundance * item.mass_u for item in isotopes_of(element(found.z).symbol)
                     if item.abundance is not None)
    return found.abundance * found.mass_u / mean
