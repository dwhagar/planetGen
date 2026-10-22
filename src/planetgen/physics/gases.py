# planetgen/physics/gases.py

"""
The gases and compounds that can be in an atmosphere (GEN.410, docs/design/atmosphere-composition.md section 4)
================================================================================================================

The one constants table of the gases and volatile compounds the atmosphere chain reads: formula, name, CAS number, molar
mass, the flags that decide what a gas does (infrared absorber, near-ultraviolet absorber, volatility class), the lowest
photolysis threshold, triple, boiling and critical points, ideal-gas heat capacity at 100, 300 and 1000 K, latent heats,
and the Rayleigh cross section with the reference temperature its refractivity fit was made at. Every row says where it
comes from (`source`), under what licence (`licence`, `"open"` where none is found) and how far it is checked (`state`:
`verified`, `partly verified` with its `limit`, `open`). It is a constants table and not a possibility table
(`generation/possibility_table.py`): it has no weights and no class gate; the recipes that weight these gases (GEN.411,
GEN.454 to GEN.456) and the class gate belong to them, and cite these rows by `formula`.

The data are `physics/data/gases.json`, built by `scripts/build_gas_table.py` from the tables of the research parts
(`docs/design/atmosphere-composition-parts/`, parts 01 and 02) and checked in with the values read from them. The loader
refuses a file or a row without the fields above (`DataError`), a duplicate formula, a formula it cannot read, and a molar
mass of the research part that is more than 0.01 g/mol from the sum of the element table (`physics.elements`, GEN.438).

Units: molar mass in kg/mol (computed here from `physics.elements`; the research part's own figure is `research_molar_mass`),
temperatures in K, pressures in Pa, heat capacity in J/(mol K), latent heat in kJ/mol, Rayleigh cross section in cm2 per
molecule (the unit the optics tables are published in). A number the research does not give is `None`.

Use:

    gas("CO2").molar_mass            # kg/mol
    gas("CO2").rayleigh.sigma_cm2[550]
    stored("h2o")                    # the row of a stored partial (O2, CO2, CO, N2, Ar, H2, H2O, CH4, H2S, SO2)
    GASES                            # every row, in the research table's order
"""

import json
import math
import os
from dataclasses import dataclass

from planetgen.generation.possibility_table import OPEN_STATE, PARTLY_VERIFIED, STATES
from planetgen.physics import elements
from planetgen.physics.elements import DataError

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "gases.json")
"""str: The data file."""

MOLAR_MASS_TOLERANCE_G_MOL = 0.01
"""float: How far the research part's molar mass may be from the element table's sum (70 of its 72 agree to this)."""

STORED_PARTIALS = ("o2", "co2", "co", "n2", "ar", "h2", "h2o", "ch4", "h2s", "so2")
"""tuple: The ten partial pressures the planet stores (`habitability_world.MOLAR_KG`), lower case."""

_IR = (None, "none", "cia", "weak", "ir")
_VOLATILITY = (None, "UV", "SV", "V", "SmV", "Ref")
_WAVELENGTHS = (450, 550, 650)


@dataclass(frozen=True)
class Rayleigh:
    """
    The Rayleigh scattering of one gas (part 02 table 2.1).

    Attributes:
        sigma_cm2 (dict): Cross section per molecule at 450, 550 and 650 nm, cm2.
        king_factor (float): The King correction factor F_K at 550 nm, (6 + 3 rho)/(6 - 7 rho).
        king_factor_assumed (bool): F_K was not found and is assumed (1.00); the cross section is good without it to a few percent.
        reference_temperature_k (float or None): The temperature the refractivity fit was made at (GEN.408: mixing conventions
            is 1.1 percent off for O2 and about 8 percent for codes that do); `None` where the research gives none (the
            polarizability route has none, and water vapour's is unknown).
        state (str): One of `STATES`.
        limit (str): For a partly verified entry, the limit.
        note (str): How the number was made.
        research_row (str): The row of part 02 table 2.1 the numbers are read from.
    """

    sigma_cm2: dict
    king_factor: float
    king_factor_assumed: bool
    reference_temperature_k: float | None
    state: str
    limit: str
    note: str
    research_row: str


@dataclass(frozen=True)
class Gas:
    """
    One gas or volatile compound.

    Attributes:
        formula (str): As the research table writes it (`"CO2"`, `"(CH3)2S"`, `"3He"`); the key everything else cites.
        name (str): The English name.
        cas (str or None): The CAS registry number.
        molar_mass (float): kg/mol, from the element table (isotope masses for `isotopes`).
        research_molar_mass (float): g/mol as the research part gives it (kept for the audit).
        ir (str or None): `"ir"` an infrared absorber (greenhouse candidate), `"cia"` collision-induced absorption only
            (homonuclear diatomics), `"weak"`, `"none"` (monatomic); `None` where not read.
        volatility (str or None): `UV` below a normal boiling or sublimation point of 30 K, `SV` 30 to 150 K, `V` 150 to 350 K,
            `SmV` 350 to 1000 K, `Ref` above 1000 K.
        photolysis_threshold_nm (float or None): The longest wavelength that splits the molecule, and `photolysis_channel`
            its products.
        absorbs_near_uv (bool or None): Absorbs between the photolysis threshold and the visible.
        isotopes (dict or None): `{isotope name: count}` for a species that is a sum of isotope masses.
        triple_point_k, triple_point_pa, boiling_point_k, critical_temperature_k, critical_pressure_pa, acentric_factor:
            Phase constants, `None` where not reached; `boiling_point_is_sublimation` for CO2-like gases.
        cp_j_mol_k (dict): Ideal-gas heat capacity at 100, 300 and 1000 K (keys 100, 300, 1000).
        gamma_300 (float or None): Cp/(Cp - R) at 300 K.
        latent_heat_vaporisation_kj_mol, latent_heat_sublimation_kj_mol: A number, a `[low, high]` span, or `None`.
        rayleigh (Rayleigh or None): The Rayleigh entry, for the gases that have one.
        source (str), licence (str), state (str), limit (str): Where the row comes from and how far it is checked.
        statuses (dict): The research part's status word of each block (`identity`, `phase points`, `heat capacity`).
        data (dict): The row as the file holds it, for the columns above that are not named fields.
    """

    formula: str
    name: str
    cas: str | None
    molar_mass: float
    research_molar_mass: float
    ir: str | None
    volatility: str | None
    photolysis_threshold_nm: float | None
    photolysis_channel: str | None
    absorbs_near_uv: bool | None
    isotopes: dict | None
    triple_point_k: float | None
    triple_point_pa: float | None
    boiling_point_k: float | None
    boiling_point_is_sublimation: bool
    critical_temperature_k: float | None
    critical_pressure_pa: float | None
    acentric_factor: float | None
    cp_j_mol_k: dict
    gamma_300: float | None
    latent_heat_vaporisation_kj_mol: object
    latent_heat_sublimation_kj_mol: object
    rayleigh: Rayleigh | None
    source: str
    licence: str
    state: str
    limit: str
    statuses: dict
    data: dict

    @property
    def is_infrared_active(self):
        """Whether the gas can be a greenhouse gas by its own bands (`ir == "ir"` or `"weak"`)."""
        return self.ir in ("ir", "weak")


def _molar_mass(formula, isotopes):
    if isotopes:
        return math.fsum(elements.isotope(name).molar_mass * count for name, count in isotopes.items())
    return elements.molar_mass(formula)


def _rayleigh(formula, item):
    if item is None:
        return None
    try:
        sigma = {int(k): float(v) for k, v in item["sigma_cm2"].items()}
        if tuple(sorted(sigma)) != _WAVELENGTHS:
            raise DataError(f"gases.json: {formula} Rayleigh sigma is not given at 450, 550 and 650 nm")
        if not all(math.isfinite(v) and v > 0 for v in sigma.values()):
            raise DataError(f"gases.json: {formula} Rayleigh sigma is not above zero")
        entry = Rayleigh(
            sigma_cm2=sigma, king_factor=float(item["king_factor_550"]), king_factor_assumed=bool(item["king_factor_assumed"]),
            reference_temperature_k=item["reference_temperature_k"], state=item["state"], limit=item["limit"], note=item["note"],
            research_row=item["research_row"])
    except KeyError as exc:
        raise DataError(f"gases.json: {formula} Rayleigh entry lacks {exc}") from None
    _check_state(f"{formula} Rayleigh", entry.state, entry.limit)
    return entry


def _check_state(label, state, limit):
    if state not in STATES:
        raise DataError(f"gases.json: {label} state {state!r} is none of {', '.join(STATES)}")
    if state in (PARTLY_VERIFIED, OPEN_STATE) and not str(limit).strip():
        raise DataError(f"gases.json: {label} is {state} and names no limit")


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    for field in ("sources", "licence", "retrieved", "title"):
        if not data.get(field):
            raise DataError(f"gases.json states no {field}")
    rows = []
    seen = set()
    for item in data.get("gases") or ():
        formula = item.get("formula")
        if not formula or formula in seen:
            raise DataError(f"gases.json: a row has no formula, or {formula!r} appears twice")
        seen.add(formula)
        for field in ("name", "source", "licence", "state", "molar_mass_research_g_mol"):
            if item.get(field) in (None, ""):
                raise DataError(f"gases.json: {formula} has no {field} (write 'open' where no licence is found)")
        _check_state(formula, item["state"], item.get("limit", ""))
        if item.get("ir") not in _IR or item.get("volatility") not in _VOLATILITY:
            raise DataError(f"gases.json: {formula} has an unknown infrared flag or volatility class")
        isotopes = item.get("isotopes")
        try:
            mass = _molar_mass(formula, isotopes)
        except (ValueError, KeyError) as exc:
            raise DataError(f"gases.json: cannot compute the molar mass of {formula}: {exc}") from None
        if abs(mass * 1000.0 - item["molar_mass_research_g_mol"]) > MOLAR_MASS_TOLERANCE_G_MOL:
            raise DataError(f"gases.json: {formula} molar mass {item['molar_mass_research_g_mol']} g/mol is more than "
                            f"{MOLAR_MASS_TOLERANCE_G_MOL} from the element table's {mass * 1000.0:.4f}")
        cp = {int(k): v for k, v in (item.get("cp_j_mol_k") or {}).items()}
        rows.append(Gas(
            formula=formula, name=item["name"], cas=item.get("cas"), molar_mass=mass,
            research_molar_mass=float(item["molar_mass_research_g_mol"]), ir=item.get("ir"), volatility=item.get("volatility"),
            photolysis_threshold_nm=item.get("photolysis_threshold_nm"), photolysis_channel=item.get("photolysis_channel"),
            absorbs_near_uv=item.get("absorbs_near_uv"), isotopes=isotopes, triple_point_k=item.get("triple_point_k"),
            triple_point_pa=item.get("triple_point_pa"), boiling_point_k=item.get("boiling_point_k"),
            boiling_point_is_sublimation=bool(item.get("boiling_point_is_sublimation")),
            critical_temperature_k=item.get("critical_temperature_k"), critical_pressure_pa=item.get("critical_pressure_pa"),
            acentric_factor=item.get("acentric_factor"), cp_j_mol_k={k: cp.get(k) for k in (100, 300, 1000)},
            gamma_300=item.get("gamma_300"), latent_heat_vaporisation_kj_mol=item.get("latent_heat_vaporisation_kj_mol"),
            latent_heat_sublimation_kj_mol=item.get("latent_heat_sublimation_kj_mol"), rayleigh=_rayleigh(formula, item.get("rayleigh")),
            source=item["source"], licence=item["licence"], state=item["state"], limit=item.get("limit", ""),
            statuses={name: item[key] for name, key in (("identity", "identity_status"), ("phase points", "phase_status"),
                                                        ("heat capacity", "thermal_status")) if item.get(key)},
            data=item))
    if not rows:
        raise DataError("gases.json holds no gases")
    return tuple(rows)


GASES = _load()
"""tuple[Gas]: Every gas and compound, in the research table's order."""

_BY_FORMULA = {row.formula: row for row in GASES}
_BY_FORMULA_LOWER = {row.formula.lower(): row for row in GASES}


class UnknownGas(KeyError):
    """A formula that is not a row of the table."""


def gas(formula):
    """
    The row of `formula` as the table writes it (`"CO2"`, `"(CH3)2S"`).

    Raises:
        UnknownGas: The table has no such gas.
    """
    try:
        return _BY_FORMULA[formula]
    except KeyError:
        raise UnknownGas(f"{formula!r} is not in the gas table (physics/data/gases.json); add it with its source") from None


def stored(key):
    """The row of a stored partial (`STORED_PARTIALS`, lower case: `"h2o"`)."""
    if key not in STORED_PARTIALS:
        raise UnknownGas(f"{key!r} is not a stored partial ({', '.join(STORED_PARTIALS)})")
    return _BY_FORMULA_LOWER[key]


def has(formula):
    """Whether the table has the gas, as written."""
    return formula in _BY_FORMULA
