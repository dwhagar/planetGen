# planetgen/physics/gases.py

"""
The gases and compounds that can be in an atmosphere (GEN.410, docs/design/composition/atmosphere-composition.md section 4)
================================================================================================================

The one constants table of the gases and volatile compounds the atmosphere chain reads: formula, name, CAS number, molar
mass, the flags that decide what a gas does (infrared absorber, near-ultraviolet absorber, volatility class), the lowest
photolysis threshold, triple, boiling and critical points, ideal-gas heat capacity at 100, 300 and 1000 K, latent heats,
the Rayleigh cross section with the reference temperature its refractivity fit was made at, the Henry constant of the
gases that dissolve in water and the Lennard-Jones constants that give the viscosity. Every row says where it
comes from (`source`) and how far it is checked (`state`:
`verified`, `partly verified` with its `limit`, `open`). It is a constants table and not a possibility table
(`generation/possibility_table.py`): it has no weights and no class gate; the recipes that weight these gases (GEN.411,
GEN.454 to GEN.456) and the class gate belong to them, and cite these rows by `formula`.

The data are `physics/data/gases.json`, built by `scripts/build_gas_table.py` from the tables of the research parts
(`docs/design/composition/atmosphere-composition-parts/`, parts 01, 02 and 03) and checked in with the values read from them. The loader
refuses a file or a row without the fields above (`DataError`), a duplicate formula, a formula it cannot read, and a molar
mass of the research part that is more than 0.01 g/mol from the sum of the element table (`physics.elements`, GEN.438).

Units: molar mass in kg/mol (computed here from `physics.elements`; the research part's own figure is `research_molar_mass`),
temperatures in K, pressures in Pa, heat capacity in J/(mol K), latent heat in kJ/mol, Rayleigh cross section in cm2 per
molecule (the unit the optics tables are published in), Henry constants in mol/(kg bar), viscosity in uPa s, band centres in cm-1
(`Band.wavelength_um` gives the wavelength), the ultraviolet cross sections in cm2 per molecule over the bins `UV_BINS_NM` (nm). A number the research does not give is `None`.

Use:

    gas("CO2").molar_mass            # kg/mol
    gas("CO2").rayleigh.sigma_cm2[550]
    gas("CO2").henry.kh(288.15)      # mol/(kg bar), 273.15 to 323.15 K only
    gas("N2").viscosity_upa_s(300)   # Chapman-Enskog, the non-polar gases
    gas("O3").uv.at(255.0)           # Hartley band, bin mean in cm2 per molecule
    gas("CO2").infrared.bands        # the band centres the research gives
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

HENRY_REFERENCE_K = 298.15
"""float: The temperature the Henry constants are tabulated at."""
HENRY_RANGE_K = (273.15, 323.15)
"""tuple: The temperatures the research gives the Henry temperature law for (dilute physical solubility)."""
NEUFELD_RANGE = (0.3, 100.0)
"""tuple: The reduced temperatures kT/epsilon the Neufeld collision integral fit is written for."""


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


class NoData(ValueError):
    """A property the table has no law for at the temperature asked, or for the gas."""


@dataclass(frozen=True)
class Henry:
    """
    The Henry constant of one gas in water (part 01 section 8): `kH(T) = kH(298.15 K) exp(d (1/T - 1/298.15))`.

    Attributes:
        kh_298 (float): mol/(kg bar) at 298.15 K, the median of the NIST WebBook compilation [R3909].
        kh_range (tuple): The lowest and highest compiled values, and `kh_n` how many there are.
        d_k (float): -dln kH/d(1/T), K, the compilation's median.
        thermo_kh_298, thermo_d_k (float): The independent fit of the `thermo` package, for the audit.
        dissociates (bool): The gas dissolves more with pH or dissociation (CO2, SO2, NH3, H2S, HCl, HCN); the constant is the
            physical solubility only.
        state (str): One of `STATES`.
        limit (str): For a partly verified entry, the limit.
    """

    kh_298: float
    kh_range: tuple
    kh_n: int
    d_k: float
    thermo_kh_298: float
    thermo_d_k: float
    dissociates: bool
    state: str
    limit: str

    def kh(self, temperature_k):
        """
        The Henry constant at `temperature_k`, mol/(kg bar).

        Raises:
            NoData: The temperature is outside 273.15 to 323.15 K, where the research gives the law.
        """
        if not HENRY_RANGE_K[0] <= temperature_k <= HENRY_RANGE_K[1]:
            raise NoData(f"the Henry law is given for {HENRY_RANGE_K[0]} to {HENRY_RANGE_K[1]} K, not {temperature_k} K")
        return self.kh_298 * math.exp(self.d_k * (1.0 / temperature_k - 1.0 / HENRY_REFERENCE_K))


@dataclass(frozen=True)
class Viscosity:
    """
    The constants that give the dilute-gas viscosity of one gas (part 01 section 9).

    Attributes:
        sigma_angstrom (float): Lennard-Jones collision diameter.
        epsilon_over_k (float): Lennard-Jones well depth over Boltzmann's constant, K.
        research_mu_300_upa_s (float): The research part's own Chapman-Enskog value at 300 K (the test checks the formula against it).
        nist_300_upa_s (float or None): The NIST fluid value at 300 K, where one was read.
        difference_percent (float or None): The research part's Chapman-Enskog against NIST at 300 K.
        polar (bool): The gas is polar and Chapman-Enskog is 10 percent off (water vapour): `Gas.viscosity_upa_s` refuses it.
        state (str): One of `STATES`.
        limit (str): For a partly verified entry, the limit.
    """

    sigma_angstrom: float
    epsilon_over_k: float
    research_mu_300_upa_s: float
    nist_300_upa_s: float | None
    difference_percent: float | None
    polar: bool
    state: str
    limit: str


def _omega22(reduced_temperature):
    """The Neufeld collision integral of a reduced temperature kT/epsilon (part 01 section 9)."""
    return (1.16145 * reduced_temperature ** -0.14874 + 0.52487 * math.exp(-0.7732 * reduced_temperature)
            + 2.16178 * math.exp(-2.43787 * reduced_temperature))


@dataclass(frozen=True)
class Band:
    """
    One infrared band as part 03 table 2.1 writes it.

    Attributes:
        label (str): The research's name of the band (`"nu2 bend"`, `"hot band"`, the molecule's own name for a row of several).
        wavenumber_cm1 (float or None): The band centre; `None` for a band the research gives only as a range.
        range_cm1 (tuple or None): The lowest and highest wavenumber of a band given as a range (the pure rotation band's edge is
            `(0, 600)`, written "below about 600 cm-1").
        approximate (bool): The research writes "about" (the value is recalled, not read from a table).
    """

    label: str
    wavenumber_cm1: float | None
    range_cm1: tuple | None
    approximate: bool

    @property
    def wavelength_um(self):
        """The wavelength of the band centre, um (`None` for a band given as a range)."""
        return None if self.wavenumber_cm1 is None else 1.0e4 / self.wavenumber_cm1


@dataclass(frozen=True)
class Infrared:
    """
    The infrared bands of one gas (part 03 table 2.1).

    Attributes:
        bands (tuple[Band]): The fundamentals and the other bands the research lists; empty for a gas with none (the noble gases,
            N2, O2 and H2 have only collision-induced bands, the metals and oxides electronic ones, and HNO3 has none written).
        saturated_rank (int or None): The gas's place in the [R0620] ranking of saturated blackbody absorption
            (1 is the strongest per molecule), where the research states it.
        strength_text, fundamentals_text (str): The research's own strength and band cells, kept whole because they carry what
            is not a number (forcing statements, overlaps, broadening).
        research_row (str): The row of table 2.1 (one row serves several gases).
        state (str), limit (str): How far the band centres are checked (recalled values are partly verified).
    """

    bands: tuple
    saturated_rank: int | None
    strength_text: str
    fundamentals_text: str
    research_row: str
    state: str
    limit: str


UV_BINS_NM = ()
"""tuple: The 35 wavelength bins `(low, high)`, nm, of the ultraviolet tables (part 02 section 2.3), from 10 to 1100 nm."""


@dataclass(frozen=True)
class UvAbsorption:
    """
    The mean total photoabsorption cross section of one gas over the bins `UV_BINS_NM` (part 02 section 2.3), cm2 per molecule:
    each cell is the mean over the part of the bin the data cover.

    Attributes:
        cross_section_cm2 (tuple): One value per bin; `None` where the bin lies outside the data (the cross section is treated as
            zero there, as the source package does) and `0.0` where the file holds zero.
        partial_bin (tuple): One flag per bin: under half of the bin is covered, so the mean is of less than half of it.
        coverage_nm (tuple): The wavelengths the data cover, `(low, high)`.
        origin (str): Where the data come from (the research names it where it differs from [R1873] or [R2040]).
        state (str), limit (str): How far the bins are checked (the structured O2 and the CH4 visible bins carry a limit).
    """

    cross_section_cm2: tuple
    partial_bin: tuple
    coverage_nm: tuple
    origin: str
    state: str
    limit: str

    def at(self, wavelength_nm):
        """
        The bin mean at `wavelength_nm`, cm2 per molecule (`0.0` outside the data).

        Raises:
            NoData: The wavelength is outside 10 to 1100 nm, the range of the bins.
        """
        for (low, high), value in zip(UV_BINS_NM, self.cross_section_cm2):
            if low <= wavelength_nm < high or (wavelength_nm == high == UV_BINS_NM[-1][1]):
                return 0.0 if value is None else value
        raise NoData(f"the ultraviolet bins cover {UV_BINS_NM[0][0]} to {UV_BINS_NM[-1][1]} nm, not {wavelength_nm} nm")

    def longest_wavelength_above(self, cross_section_cm2):
        """The upper edge of the last bin whose mean is at least `cross_section_cm2`, nm, or `None` when no bin is."""
        edges = [high for (_low, high), value in zip(UV_BINS_NM, self.cross_section_cm2)
                 if value is not None and value >= cross_section_cm2]
        return max(edges) if edges else None


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
        henry (Henry or None): The Henry constant, for the gases that dissolve in water and were compiled.
        viscosity (Viscosity or None): The Lennard-Jones constants, for the gases the research tabulates.
        infrared (Infrared or None): The infrared bands, for the gases part 03 lists.
        uv (UvAbsorption or None): The ultraviolet and visible absorption bins, for the 22 gases part 02 tabulates.
        source (str), state (str), limit (str): Where the row comes from and how far it is checked.
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
    henry: Henry | None
    viscosity: Viscosity | None
    infrared: Infrared | None
    uv: UvAbsorption | None
    source: str
    state: str
    limit: str
    statuses: dict
    data: dict

    @property
    def is_infrared_active(self):
        """Whether the gas can be a greenhouse gas by its own bands (`ir == "ir"` or `"weak"`)."""
        return self.ir in ("ir", "weak")

    def viscosity_upa_s(self, temperature_k):
        """
        The dilute-gas viscosity in uPa s by Chapman-Enskog with the Neufeld collision integral:
        `mu = 2.6693 sqrt(M T)/(sigma^2 Omega22)`, M in g/mol, T in K, sigma in Angstrom.

        Raises:
            NoData: The gas has no Lennard-Jones constants, is polar (water vapour: the research gives only the NIST value at
                300 K), or the reduced temperature kT/epsilon is outside 0.3 to 100.
        """
        entry = self.viscosity
        if entry is None:
            raise NoData(f"the table has no viscosity constants for {self.formula}")
        if entry.polar:
            raise NoData(f"{self.formula} is polar: Chapman-Enskog is 10 percent off and the research has the NIST value at 300 K only "
                         f"({entry.nist_300_upa_s} uPa s)")
        reduced = temperature_k / entry.epsilon_over_k
        if not NEUFELD_RANGE[0] <= reduced <= NEUFELD_RANGE[1]:
            raise NoData(f"{self.formula} at {temperature_k} K has kT/epsilon {reduced:.3g}, outside the Neufeld fit's "
                         f"{NEUFELD_RANGE[0]} to {NEUFELD_RANGE[1]}")
        return 2.6693 * math.sqrt(self.molar_mass * 1000.0 * temperature_k) / (entry.sigma_angstrom ** 2 * _omega22(reduced))


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


def _henry(formula, item):
    if item is None:
        return None
    try:
        low, high = item["kh_range_mol_kg_bar"]
        entry = Henry(
            kh_298=float(item["kh_298_mol_kg_bar"]), kh_range=(float(low), float(high)), kh_n=int(item["kh_n"]), d_k=float(item["d_k"]),
            thermo_kh_298=float(item["thermo_kh_298_mol_kg_bar"]), thermo_d_k=float(item["thermo_d_k"]),
            dissociates=bool(item["dissociates"]), state=item["state"], limit=item["limit"])
    except KeyError as exc:
        raise DataError(f"gases.json: {formula} Henry entry lacks {exc}") from None
    if not (math.isfinite(entry.kh_298) and entry.kh_298 > 0 and 0 < entry.kh_range[0] <= entry.kh_range[1]):
        raise DataError(f"gases.json: {formula} Henry constant or its range is not above zero")
    _check_state(f"{formula} Henry", entry.state, entry.limit)
    return entry


def _viscosity(formula, item):
    if item is None:
        return None
    try:
        entry = Viscosity(
            sigma_angstrom=float(item["sigma_angstrom"]), epsilon_over_k=float(item["epsilon_over_k_k"]),
            research_mu_300_upa_s=float(item["mu_300_research_upa_s"]), nist_300_upa_s=item["nist_300_upa_s"],
            difference_percent=item["difference_percent"], polar=bool(item["polar"]), state=item["state"], limit=item["limit"])
    except KeyError as exc:
        raise DataError(f"gases.json: {formula} viscosity entry lacks {exc}") from None
    if not (entry.sigma_angstrom > 0 and entry.epsilon_over_k > 0):
        raise DataError(f"gases.json: {formula} Lennard-Jones constants are not above zero")
    _check_state(f"{formula} viscosity", entry.state, entry.limit)
    return entry


def _infrared(formula, item):
    if item is None:
        return None
    try:
        bands = []
        for band in item["bands"]:
            centre, span = band["wavenumber_cm1"], band["range_cm1"]
            if (centre is None) == (span is None):
                raise DataError(f"gases.json: {formula} has a band with a centre and a range, or with neither")
            if (centre is not None and not centre > 0) or (span is not None and not 0 <= span[0] < span[1]):
                raise DataError(f"gases.json: {formula} has a band at no positive wavenumber")
            bands.append(Band(band["label"], centre, None if span is None else (float(span[0]), float(span[1])),
                              bool(band["approximate"])))
        entry = Infrared(bands=tuple(bands), saturated_rank=item["saturated_rank"], strength_text=item["strength_text"],
                         fundamentals_text=item["fundamentals_text"], research_row=item["research_row"], state=item["state"],
                         limit=item["limit"])
    except KeyError as exc:
        raise DataError(f"gases.json: {formula} infrared entry lacks {exc}") from None
    if bool(entry.bands) == bool(item["no_band"]):
        raise DataError(f"gases.json: {formula} infrared entry says it has {'no band' if item['no_band'] else 'bands'} and lists the opposite")
    _check_state(f"{formula} infrared", entry.state, entry.limit)
    return entry


def _uv(formula, item, bins):
    if item is None:
        return None
    try:
        values = tuple(None if value is None else float(value) for value in item["cross_section_cm2"])
        entry = UvAbsorption(cross_section_cm2=values, partial_bin=tuple(bool(flag) for flag in item["partial_bin"]),
                             coverage_nm=tuple(float(edge) for edge in item["coverage_nm"]), origin=item["origin"],
                             state=item["state"], limit=item["limit"])
    except KeyError as exc:
        raise DataError(f"gases.json: {formula} ultraviolet entry lacks {exc}") from None
    if len(values) != len(bins) or len(entry.partial_bin) != len(bins):
        raise DataError(f"gases.json: {formula} ultraviolet entry has {len(values)} values for {len(bins)} bins")
    if not all(value is None or (math.isfinite(value) and value >= 0.0) for value in values):
        raise DataError(f"gases.json: {formula} ultraviolet cross section is negative or not finite")
    if not entry.coverage_nm[0] < entry.coverage_nm[1]:
        raise DataError(f"gases.json: {formula} ultraviolet coverage is not an interval")
    _check_state(f"{formula} ultraviolet", entry.state, entry.limit)
    return entry


def _check_state(label, state, limit):
    if state not in STATES:
        raise DataError(f"gases.json: {label} state {state!r} is none of {', '.join(STATES)}")
    if state in (PARTLY_VERIFIED, OPEN_STATE) and not str(limit).strip():
        raise DataError(f"gases.json: {label} is {state} and names no limit")


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    for field in ("sources", "retrieved", "title"):
        if not data.get(field):
            raise DataError(f"gases.json states no {field}")
    global UV_BINS_NM
    bins = tuple((float(low), float(high)) for low, high in data.get("uv_bins_nm") or ())
    if not bins or any(not low < high for low, high in bins) or any(a[1] != b[0] for a, b in zip(bins, bins[1:])):
        raise DataError("gases.json: the ultraviolet bins are missing, empty, or not one after another")
    UV_BINS_NM = bins
    rows = []
    seen = set()
    for item in data.get("gases") or ():
        formula = item.get("formula")
        if not formula or formula in seen:
            raise DataError(f"gases.json: a row has no formula, or {formula!r} appears twice")
        seen.add(formula)
        for field in ("name", "source", "state", "molar_mass_research_g_mol"):
            if item.get(field) in (None, ""):
                raise DataError(f"gases.json: {formula} has no {field}")
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
            henry=_henry(formula, item.get("henry")), viscosity=_viscosity(formula, item.get("viscosity")),
            infrared=_infrared(formula, item.get("infrared")), uv=_uv(formula, item.get("uv"), bins),
            source=item["source"], state=item["state"], limit=item.get("limit", ""),
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
