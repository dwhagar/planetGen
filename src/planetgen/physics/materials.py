# planetgen/physics/materials.py

"""
The shared materials and components table (GEN.427, docs/design/composition/bulk-composition.md section 3, composition-model.md section 4b)
================================================================================================================================

The one table of every element, mineral, ice, organic and alloy a body can contain: its formula, molar mass, grain density,
phase points, latent heats, the band albedos and the equation-of-state constants of the model GEN.426 chose, each number with
the source it was read from. There is no component list owned by one kind: the 84 names of `tuning.ASTEROID_COMPONENTS`, the
taxa of `physics/small_body_composition.py` and the composition texts of
`tuning.PLANET_CLASSES` all resolve to rows of it (`material`, `class_composition`), and a test fails when one does not. A
planet's or moon's bulk and surface draws (GEN.428, GEN.441) read the same rows; this item holds the data only and
changes no draw, so it bumps nothing (Boss, 2:10 PM EDT, card "Data table, no bump").

It is a constants table with the gate of the possibility-table format (`generation/possibility_table.py`, GEN.472): each
row says which classes can hold it (`gate`: `any` for every planet and moon class and every small-body kind, `planetary` for the
planet and moon classes only, which is where the high-pressure silicate phases, the high-pressure ices, iron alloys of a core
and the hydrogen-helium envelopes sit), a source and a source state. The licence of the data is not a column: it is recorded in `DATA-LICENSE.md`
at the repository root (Boss, 2026-10-22). The rows are
registered as the table `materials` (`TABLE`, `generation.possibility_table.SHARED_COMPONENTS`), so the fractions of any other
table refer to them and `REGISTRY.validate` checks it. Its weights are 1 and unused: nothing is drawn from this table.

The data are `physics/data/materials.json`, built by `scripts/build_materials_table.py` from part E's materials table, part B's
SLB phases and iron rows, part C's ice phases and the values in `scripts/materials_data.py` ([R3992: Tables 1 to 3],
the Handbook of Mineralogy, PubChem), and checked in. The loader refuses (`DataError`) a file without its sources, a
row without a source and a state, a number without a source key, a partly verified row without its limit, a formula it
cannot read, a gate it does not know, an equation of state without the constants of its model, an alias that is two rows' name,
a group whose member is missing, and a `labels` entry that names no property of the row. A value with no source is `None` with
the code `N` (not reached); none is kept from memory.

Units: densities in g/cm3 (`rho_grain_g_cm3`, the grain density at ambient conditions, solids at room temperature), temperatures
in K, pressures in bar, latent heats in kJ/mol and kJ/kg, interaction length in g/cm2, albedo as a fraction, molar mass in kg/mol
(computed here from `physics.elements`). An equation of state states its own units in the names of its constants (`rho0_kg_m3`,
`k0_gpa`).

Use:

    material("forsterite").molar_mass                      # kg/mol
    material("plagioclase feldspars").members              # the rows a group name stands for
    material("carbon dioxide ice").prop("rho_grain_g_cm3").value
    material("iron_epsilon").eos_of("vinet")               # the Vinet constants, source and range
    class_composition("silicon, iron, and magnesium")     # the rows a class text names
    MATERIALS                                              # every row, in the file's order
"""

import json
import math
import os
from dataclasses import dataclass

from planetgen import tuning
from planetgen.generation import possibility_table as pt
from planetgen.physics import elements
from planetgen.physics.elements import DataError

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "materials.json")
"""str: The data file."""

ANY = "any"
PLANETARY = "planetary"
SMALL_BODY_KINDS = ("asteroid-field", "belt", "comet")
"""tuple: The small-body kinds a gate of `any` includes beside the planet and moon classes."""

PROPERTIES = {
    "rho_grain_g_cm3": ("g/cm3", "rho"),
    "t_melt_k": ("K", "Tm"),
    "t_boil_or_sub_k": ("K", "Tb"),
    "t_triple_k": ("K", "Tt"),
    "p_triple_bar": ("bar", "Pt"),
    "t_crit_k": ("K", "Tc"),
    "p_crit_bar": ("bar", "Pc"),
    "h_fus_kj_mol": ("kJ/mol", "Hfus"),
    "h_vap_kj_mol": ("kJ/mol", "Hvap"),
    "h_sub_kj_mol": ("kJ/mol", "Hsub"),
    "h_fus_kj_kg": ("kJ/kg", "Hfus_kg"),
    "lambda_i_g_cm2": ("g/cm2", "lamI"),
    "debye_k": ("K", "debye"),
    "alb_040_050": ("1", "alb"),
    "alb_050_060": ("1", "alb"),
    "alb_060_070": ("1", "alb"),
    "rho_gas_stp_g_l": ("g/L", "rhog"),
}
"""dict: Every property a row may carry: key -> (unit, the key of part E's label column). The albedo of the three bands share one label."""

CODES = ("R", "C", "H", "S", "N")
"""tuple: Part E's label codes: R read, C computed, H second-hand, S snippet, N not reached."""
CODE_STATE = {"R": pt.VERIFIED, "C": pt.VERIFIED, "H": pt.PARTLY_VERIFIED, "S": pt.PARTLY_VERIFIED, "N": pt.OPEN_STATE}
"""dict: The source state a code gives a value."""

EOS_MODELS = {
    "vinet": ("rho0_kg_m3", "k0_gpa", "k0p"),
    "bme": ("k0_gpa", "k0p", "order"),
    "holzapfel": ("rho0_kg_m3", "k0_gpa", "k0p"),
    "mgd-slb3": ("v0_cm3_mol", "rho0_g_cm3", "k0_gpa", "k0p", "theta0_k", "gamma0", "q0", "unit_formula", "unit_molar_mass_g_mol"),
    "modified-polytrope": ("rho0_kg_m3", "c_kg_m3_pa_n", "n"),
    "tabulated": ("columns", "points"),
    "reference-points": ("columns", "points"),
    "ideal-gas": (),
}
"""dict: The equation-of-state models and the constants each must give (`bme` also gives `rho0_kg_m3` or `v0_cm3_mol`, and
`k0pp_per_gpa` at order 4)."""

LICENCE_FILE = "DATA-LICENSE.md"
"""str: The file in the repository root where the licence of the data is recorded."""

UNIT_MASS_TOLERANCE_G_MOL = 0.05
"""float: How far an SLB unit's printed molar mass may be from the element table's sum (the research used burnman's atomic weights)."""


@dataclass(frozen=True)
class Prop:
    """
    One number of a row.

    Attributes:
        value (float or None): The number; `None` where none was reached (`code` `"N"`).
        unit (str): Its unit (`PROPERTIES`).
        source (str or None): The key of its source in `SOURCES`.
        code (str): One of `CODES`.
        state (str): The source state the code gives (`CODE_STATE`).
        measured_range (tuple or None): For a density, the measured range the compilation gives.
        basis (str): How the number was taken (D(calc.), the temperature it was measured at).
    """

    value: float | None
    unit: str
    source: str | None
    code: str
    state: str
    measured_range: tuple | None
    basis: str


@dataclass(frozen=True)
class Eos:
    """
    One equation of state of a row.

    Attributes:
        model (str): A key of `EOS_MODELS`.
        params (dict): Its constants, with the unit in each name.
        valid (dict): The range it is valid over (`p_gpa`, `p_mpa`, `p_pa`, `t_k`), each `[low, high]`; empty where the source gives none.
        uncertainty (dict): The stated uncertainty of constants, by name.
        source (str): The key of its source in `SOURCES`.
        code (str): One of `CODES`.
        note (str): What the source says about it.
    """

    model: str
    params: dict
    valid: dict
    uncertainty: dict
    source: str
    code: str
    note: str


@dataclass(frozen=True)
class Material:
    """
    One row.

    Attributes:
        id (str): The row's name, lower case with underscores.
        name (str): What it is called in words.
        family (str): The kind of material (`ice`, `oxide`, `metal`, `mineral group`).
        formula (str or None): A formula `physics.elements` reads; `None` for a solid solution, a group or a mixture.
        composition (str): The composition in words where the formula is not the whole story.
        aliases (tuple): The other names the lists use for it (`tuning.ASTEROID_COMPONENTS` and the like).
        members (tuple or None): For a group, the rows it stands for; a group has no number of its own.
        gas_formula (str or None): The row of the gas table (`physics.gases`) of the same molecule, which owns its gas properties.
        molar_mass (float or None): kg/mol from the formula.
        gate (tuple): The classes that can hold it: `(pt.ANY_CLASS,)` or the planet and moon classes.
        props (dict): Property key -> `Prop`.
        eos (tuple): The `Eos` entries.
        labels (str): Part E's label column for the row (`rho=R;Tm=H`), rebuilt from the codes of its numbers.
        source, state, limit (str): The row's source summary, source state, and for a partly verified row the limit.
        note (str): What was dropped or replaced.
    """

    id: str
    name: str
    family: str
    formula: str | None
    composition: str
    aliases: tuple
    members: tuple | None
    gas_formula: str | None
    molar_mass: float | None
    gate: tuple
    props: dict
    eos: tuple
    labels: str
    source: str
    state: str
    limit: str
    note: str

    def prop(self, key):
        """The `Prop` of `key` (`"rho_grain_g_cm3"`); raises `KeyError` where the row has no such column."""
        try:
            return self.props[key]
        except KeyError:
            raise KeyError(f"the row {self.id!r} has no {key!r}") from None

    def eos_of(self, model):
        """The first `Eos` of `model` on this row; raises `KeyError` where it has none."""
        for entry in self.eos:
            if entry.model == model:
                return entry
        raise KeyError(f"the row {self.id!r} has no {model!r} equation of state")

    @property
    def is_group(self):
        return self.members is not None


class UnknownMaterial(KeyError):
    """A name that is neither a row of the table nor an alias of one."""


def _plain_classes():
    return tuple(sorted(tuning.PLANET_CLASSES))


def gate_classes(gate):
    """The classes a gate group names: `any` is every class (`pt.ANY_CLASS`), `planetary` the planet and moon classes."""
    if gate == ANY:
        return (pt.ANY_CLASS,)
    if gate == PLANETARY:
        return _plain_classes()
    raise DataError(f"materials.json: the gate {gate!r} is neither {ANY!r} nor {PLANETARY!r}")


def _number(label, value, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise DataError(f"materials.json: {label} is not a finite number")
    if positive and not value > 0:
        raise DataError(f"materials.json: {label} is not above zero")


def _prop(row_id, key, item, sources):
    if key not in PROPERTIES:
        raise DataError(f"materials.json: {row_id} has the unknown property {key!r}")
    label = f"{row_id} {key}"
    try:
        value, unit, source, code = item["value"], item["unit"], item["source"], item["code"]
    except KeyError as exc:
        raise DataError(f"materials.json: {label} lacks {exc}") from None
    if code not in CODES:
        raise DataError(f"materials.json: {label} has the code {code!r}, none of {', '.join(CODES)}")
    if unit != PROPERTIES[key][0]:
        raise DataError(f"materials.json: {label} is in {unit!r}, not {PROPERTIES[key][0]!r}")
    if value is None:
        if code != "N" or source is not None:
            raise DataError(f"materials.json: {label} has no value and must be the code N with no source")
    else:
        if code == "N":
            raise DataError(f"materials.json: {label} has a value and the code N")
        _number(label, value)
        if not source or source not in sources:
            raise DataError(f"materials.json: {label} has no source (every number has one; a value with none is null)")
        if key.startswith("alb_") and not 0.0 <= value <= 1.0:
            raise DataError(f"materials.json: {label} is not a fraction")
        if key != "debye_k" and key not in ("rho_gas_stp_g_l",) and not value >= 0:
            raise DataError(f"materials.json: {label} is negative")
    measured = item.get("measured_range")
    if measured is not None:
        low, high = measured
        _number(label + " measured range", low)
        _number(label + " measured range", high)
        if not low <= high:
            raise DataError(f"materials.json: {label} measured range runs backwards")
        measured = (float(low), float(high))
    return Prop(None if value is None else float(value), unit, source, code, CODE_STATE[code], measured, item.get("basis", ""))


def _eos(row_id, item, sources):
    label = f"{row_id} equation of state"
    try:
        model, params, source, code = item["model"], item["params"], item["source"], item["code"]
    except KeyError as exc:
        raise DataError(f"materials.json: {label} lacks {exc}") from None
    if model not in EOS_MODELS:
        raise DataError(f"materials.json: {label} is the unknown model {model!r}")
    if code not in CODES or code == "N":
        raise DataError(f"materials.json: {label} has the code {code!r}")
    if source not in sources:
        raise DataError(f"materials.json: {label} has no source")
    missing = [name for name in EOS_MODELS[model] if name not in params]
    if missing:
        raise DataError(f"materials.json: {label} ({model}) lacks the constants {', '.join(missing)}")
    if model == "bme":
        if params["order"] not in (3, 4):
            raise DataError(f"materials.json: {label}: the Birch-Murnaghan order is 3 or 4")
        if ("rho0_kg_m3" in params) == ("v0_cm3_mol" in params):
            raise DataError(f"materials.json: {label}: give the zero-pressure density or the molar volume, not both or neither")
        if params["order"] == 4 and "k0pp_per_gpa" not in params:
            raise DataError(f"materials.json: {label}: a fourth-order fit names K'' (None where the source gives none)")
    for name, value in params.items():
        if name in ("unit_formula", "columns", "points", "order") or value is None:
            continue
        _number(f"{label} {name}", value)
    if model in ("tabulated", "reference-points"):
        width = len(params["columns"])
        if not params["points"] or any(len(point) != width for point in params["points"]):
            raise DataError(f"materials.json: {label}: a point does not have {width} values")
    if model == "mgd-slb3":
        mass = elements.molar_mass(params["unit_formula"]) * 1000.0
        if abs(mass - params["unit_molar_mass_g_mol"]) > UNIT_MASS_TOLERANCE_G_MOL:
            raise DataError(f"materials.json: {label}: the unit {params['unit_formula']} has molar mass {mass:.3f} g/mol, not "
                            f"{params['unit_molar_mass_g_mol']}")
    for name, span in item.get("valid", {}).items():
        if len(span) != 2 or not span[0] < span[1]:
            raise DataError(f"materials.json: {label}: the validity range {name} is not [low, high]")
    return Eos(model, dict(params), {k: tuple(v) for k, v in item.get("valid", {}).items()}, dict(item.get("uncertainty", {})), source,
               code, item.get("note", ""))


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    for field in ("sources", "retrieved", "title"):
        if not data.get(field):
            raise DataError(f"materials.json states no {field}")
    sources = data["sources"]
    rows, seen, names = [], set(), {}
    for item in data.get("materials") or ():
        row_id = item.get("id")
        if not row_id or row_id in seen:
            raise DataError(f"materials.json: a row has no id, or {row_id!r} appears twice")
        seen.add(row_id)
        for field in ("name", "family", "gate", "source", "state"):
            if not item.get(field):
                raise DataError(f"materials.json: {row_id} has no {field}")
        if item["state"] not in pt.STATES:
            raise DataError(f"materials.json: {row_id} state {item['state']!r} is none of {', '.join(pt.STATES)}")
        if item["state"] in (pt.PARTLY_VERIFIED, pt.OPEN_STATE) and not str(item.get("limit", "")).strip():
            raise DataError(f"materials.json: {row_id} is {item['state']} and names no limit")
        gate = gate_classes(item["gate"])
        formula = item.get("formula")
        mass = None
        if formula:
            try:
                mass = elements.molar_mass(formula)
            except (ValueError, KeyError) as exc:
                raise DataError(f"materials.json: cannot read the formula of {row_id}: {exc}") from None
        members = item.get("members")
        props = {key: _prop(row_id, key, entry, sources) for key, entry in (item.get("props") or {}).items()}
        eos_entries = tuple(_eos(row_id, entry, sources) for entry in item.get("eos") or ())
        if members is not None and (props or eos_entries or formula):
            raise DataError(f"materials.json: {row_id} is a group and carries a number or formula of its own")
        if members is not None and len(members) < 2:
            raise DataError(f"materials.json: the group {row_id} names fewer than two members")
        labels = item.get("labels", "")
        for pair in filter(None, labels.split(";")):
            key, _, code = pair.partition("=")
            columns = [k for k, (_unit, label_key) in PROPERTIES.items() if label_key == key]
            if not columns or not any(k in props and props[k].code == code for k in columns):
                raise DataError(f"materials.json: {row_id} labels {pair!r} name no property of the row with that code")
        for name in [item["name"], row_id, *item.get("aliases", ())]:
            key = name.lower()
            if key in names and names[key] != row_id:
                raise DataError(f"materials.json: the name {name!r} is both {names[key]!r} and {row_id!r}")
            names[key] = row_id
        rows.append(Material(
            id=row_id, name=item["name"], family=item["family"], formula=formula, composition=item.get("composition", ""),
            aliases=tuple(item.get("aliases", ())), members=None if members is None else tuple(members),
            gas_formula=item.get("gas_formula"), molar_mass=mass, gate=gate, props=props, eos=eos_entries, labels=labels,
            source=item["source"], state=item["state"], limit=item.get("limit", ""), note=item.get("note", "")))
    if not rows:
        raise DataError("materials.json holds no materials")
    ids = {row.id for row in rows}
    for row in rows:
        for member in row.members or ():
            if member not in ids:
                raise DataError(f"materials.json: the group {row.id} names {member!r}, which is not a row")
    return data, tuple(rows), names


_DATA, MATERIALS, _NAMES = _load()
"""tuple[Material]: Every row, in the file's order."""
SOURCES = _DATA["sources"]
"""dict: The source keys and their citations."""
_BY_ID = {row.id: row for row in MATERIALS}


def material(name):
    """
    The row a name stands for: its id, its name in words or an alias (case does not matter).

    Raises:
        UnknownMaterial: The table has no such material; add it with its source.
    """
    row_id = _NAMES.get(str(name).lower())
    if row_id is None:
        raise UnknownMaterial(f"{name!r} is not in the materials table (physics/data/materials.json); add it with its source")
    return _BY_ID[row_id]


def has(name):
    """Whether the table has the material, by id, name or alias."""
    return str(name).lower() in _NAMES


def class_composition(text):
    """
    The rows a class composition text of `tuning.PLANET_CLASSES` names (`"silicon, iron, and magnesium"`), in the order the text
    gives them.

    Raises:
        UnknownMaterial: The text is not one the file maps.
    """
    try:
        names = _DATA["legacy_names"]["class_compositions"][text]
    except KeyError:
        raise UnknownMaterial(f"the class composition {text!r} is not mapped in materials.json (legacy_names.class_compositions)") from None
    return tuple(material(name) for name in names)


def expand(name):
    """The rows a name ends in: the row itself, or for a group the rows of its members, recursively."""
    row = material(name)
    if not row.is_group:
        return (row,)
    found = []
    for member in row.members:
        found.extend(expand(member))
    return tuple(found)


def _row_for_table(row):
    return pt.Row(
        name=row.id, classes=row.gate, weight=1.0, source=row.source, state=row.state, limit=row.limit,
        notes=row.note)


TABLE = pt.REGISTRY.register(pt.Table(
    pt.SHARED_COMPONENTS, 1, [_row_for_table(row) for row in MATERIALS],
    description="The shared materials and components table: every element, mineral, ice, organic and alloy a body can "
                "contain, with the class gate, the source of each row. A constants table: its weights are 1 and are not "
                "drawn from."))
"""Table: The rows as a possibility table, registered as `materials` so the fractions of other tables refer to them."""
