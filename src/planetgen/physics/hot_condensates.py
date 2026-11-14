# planetgen/physics/hot_condensates.py

"""
Cloud temperatures and trace gases that follow the host star's elements (GEN.1375)
==================================================================================

The condensation temperatures of the cloud species of hot atmospheres (chromium, manganese sulfide, sodium sulfide, zinc sulfide and potassium
chloride) from the elements of the host star, and the elements and gases the atmosphere tables leave out: the trace elements Rb, Cs, Ge, As and Se
and the gases GeH4, AsH3, H2Se, MgH, CaH, NaH, LiCl, LiF, AlCl, PN and HCP. The data are `physics/data/hot_condensates.json`.

The temperature of a cloud species falls out of one relation, [R3165, R4562] (READ):

    10^4 / T_cond = a - b log10 P_T - sum over the limiting elements of c [X/H]

with `P_T` the total pressure in bar and `[X/H]` the host's abundance of the element against the Sun's. The research part prints these with every
`[X/H]` set to `[Fe/H]`; here each takes the host's own, so a star poor in zinc for its iron has a cooler zinc sulfide cloud than one that follows
iron. At the Sun's composition the fits give the tabulated temperatures (Cr 1521 K, MnS 1343 K, Na2S 996 K and ZnS 798 K at 1 bar).
Potassium chloride has the saturation relation `log10 P_sat = 7.611 - 11382/T` for the gas pressure `2.55e-7 P_T 10^[K/H]`; potassium limits it
because the Sun holds less potassium than chlorine (5.07 against 5.31 in log eps).

The abundance of an element against hydrogen by number is `10^(log eps - 12 + [X/H])`. The elements of the star's element set take their
[X/Fe] from the host (`StellarComposition.offsets`); the elements of the table the star does not carry (Rb, Cs, Ge, As, Se) follow iron
(stand-in, `[X/Fe] = 0`; GEN.725 HK1), and chlorine and fluorine take the [Cl/Fe] and [F/Fe] of the rocky-planet inventory
(`mixture_rocky`).

A trace gas can hold no more of an element than the star has: `gas_ceiling` is the most molecules per hydrogen atom, the least of the host's
abundances of its elements divided by how many atoms of each the molecule holds (GEN.725 HK2). The real amount of each gas depends on the
temperature and pressure and the chemistry of the atmosphere (the hot-giant thermochemistry, GEN.994) and is not here.

Use:

    cloud_temperature_k("ZnS", 1.0, host)    # K, host a `StellarComposition`
    abundance_per_h("Zn", host)              # number of zinc atoms per hydrogen atom
    gas_ceiling("HCP", host)                 # most HCP molecules per hydrogen atom
"""

import json
import math
import os
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import elements, mixture_rocky
from planetgen.physics.elements import DataError, positive_number

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "hot_condensates.json")
"""str: The data file."""


@dataclass(frozen=True)
class Fit:
    """One cloud species: `10^4 / T = a - b log10 P_T - sum c [X/H]` over `terms`, `(symbol, c)` pairs."""

    species: str
    a: float
    b: float
    terms: tuple


@dataclass(frozen=True)
class TraceGas:
    """One trace gas: `atoms` maps each element to how many atoms the molecule holds; `where` says what the research gives for it."""

    formula: str
    name: str
    atoms: dict
    where: str


def load(path=DATA_FILE):
    """Reads and checks the file at `path`: `(fits, kcl, trace elements, trace gases)`; a malformed row is refused (`DataError`)."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the hot condensate table {path}: {error}") from error

    number = positive_number      # GEN.1667: the shared check

    fits = {}
    for raw in document.get("fits") or []:
        terms = tuple((symbol, number(coefficient, f"the {raw.get('species')} coefficient")) for symbol, coefficient in raw.get("terms") or [])
        if not terms or raw.get("species") in fits:
            raise DataError(f"fit {raw.get('species')!r} needs terms and a name of its own")
        for symbol, _ in terms:
            elements.element(symbol)
        fits[raw["species"]] = Fit(raw["species"], number(raw.get("a"), "a"), number(raw.get("b"), "b"), terms)
    kcl = document.get("kcl") or {}
    kcl = {"a": number(kcl.get("a"), "the KCl constant"), "b_k": number(kcl.get("b_k"), "the KCl slope"),
           "gas_per_total": number(kcl.get("gas_per_total"), "the KCl gas share"), "element": elements.element(kcl.get("element", "")).symbol}
    trace = {}
    for raw in document.get("trace_elements") or []:
        symbol = elements.element(raw.get("symbol", "")).symbol
        if not isinstance(raw.get("log_eps"), (int, float)) or symbol in trace:
            raise DataError(f"trace element {symbol!r} needs a log eps and a row of its own")
        trace[symbol] = float(raw["log_eps"])
    gases = {}
    for raw in document.get("trace_gases") or []:
        atoms = {elements.element(symbol).symbol: int(count) for symbol, count in (raw.get("atoms") or {}).items()}
        if not atoms or min(atoms.values()) < 1 or raw.get("formula") in gases:
            raise DataError(f"trace gas {raw.get('formula')!r} needs atoms and a row of its own")
        gases[raw["formula"]] = TraceGas(raw["formula"], raw.get("name", ""), atoms, raw.get("where", ""))
    if not fits or not trace or not gases:
        raise DataError("the hot condensate table is missing rows")
    return fits, kcl, trace, gases


FITS, KCL, TRACE_ELEMENTS, HOT_TRACE_GASES = load()
"""The cloud fits by species, the KCl relation, the trace elements' log eps and the trace gases by formula."""

CLOUD_SPECIES = tuple(FITS) + ("KCl",)
"""tuple: The species `cloud_temperature_k` knows."""


def _log_eps(symbol):
    if symbol in tuning.STELLAR_SOLAR_LOG_EPS:
        return tuning.STELLAR_SOLAR_LOG_EPS[symbol]
    if symbol in TRACE_ELEMENTS:
        return TRACE_ELEMENTS[symbol]
    if symbol in HALOGEN_LOG_EPS:
        return HALOGEN_LOG_EPS[symbol]
    raise ValueError(f"no solar abundance for {symbol}")


HALOGEN_LOG_EPS = {"Cl": 5.31, "F": 4.4}
"""dict: [R0168: Table 2] log eps of the halogens the star set does not carry."""


def host_xh(symbol, host):
    """[X/H] (dex) of `symbol` in `host`, a `StellarComposition`: the host's [Fe/H] plus its [X/Fe] (the halogens take `mixture_rocky`'s, the elements
    the host does not carry follow iron)."""
    offset = {"Cl": mixture_rocky.STELLAR_CL_OVER_FE_DEX, "F": mixture_rocky.STELLAR_F_OVER_FE_DEX}.get(symbol)
    if offset is None:
        offset = host.offsets.get(symbol, 0.0)
    return host.fe_h + offset


def abundance_per_h(symbol, host):
    """Atoms of `symbol` per hydrogen atom in `host`: `10^(log eps - 12 + [X/H])`."""
    return 10.0 ** (_log_eps(symbol) - 12.0 + host_xh(symbol, host))


def cloud_temperature_k(species, pressure_bar, host):
    """The temperature (K) at which `species` (one of `CLOUD_SPECIES`) starts to condense from the gas of `host` at `pressure_bar` total pressure.

    Raises:
        ValueError: Not a known species, or the pressure is not positive.
    """
    if not pressure_bar > 0.0:
        raise ValueError("the pressure must be positive")
    log_p = math.log10(pressure_bar)
    if species == "KCl":
        gas_log_p = math.log10(KCL["gas_per_total"]) + log_p + host_xh(KCL["element"], host)
        return KCL["b_k"] / (KCL["a"] - gas_log_p)
    if species not in FITS:
        raise ValueError(f"unknown cloud species {species!r}; known: {', '.join(CLOUD_SPECIES)}")
    fit = FITS[species]
    inverse = fit.a - fit.b * log_p - sum(coefficient * host_xh(symbol, host) for symbol, coefficient in fit.terms)
    return 1e4 / inverse


def gas_ceiling(formula, host):
    """The most molecules of `formula` per hydrogen atom that `host` leaves room for: the least of its elements' abundances over the atoms each
    holds. Raises `ValueError` for a gas that is not a trace gas of the table."""
    if formula not in HOT_TRACE_GASES:
        raise ValueError(f"unknown trace gas {formula!r}; known: {', '.join(HOT_TRACE_GASES)}")
    return min(abundance_per_h(symbol, host) / count for symbol, count in HOT_TRACE_GASES[formula].atoms.items())
