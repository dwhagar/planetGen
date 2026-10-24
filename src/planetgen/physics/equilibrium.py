# planetgen/physics/equilibrium.py

"""
Gas-phase chemical equilibrium by Gibbs-energy minimisation (GEN.518, docs/design/composition/equilibrium-chemistry.md)
=====================================================================================================================

Which molecules a hot atmosphere holds at a temperature, a pressure and an element budget: the free energy of an ideal-gas
mixture is minimised under the conservation of each element, with the element potentials (Lagrange multipliers) as the unknowns
(Gordon and McBride 1994, the method of the NASA CEA program, written here from the equations). The thermodynamics are the NASA
seven-coefficient polynomials of McBride, Gordon and Reno 1993 (NASA TM-4513, `physics/data/thermo_nasa7.json`, built by
`scripts/build_thermo_table.py`). Gas only: no condensates, no ions, ideal gas (valid below about 100 bar), standard state 1 bar.

What it answers for the generator: the equilibrium line between methane and carbon monoxide (`ch4_co_temperature`), the
temperature at which a share of the water is destroyed (`h2o_survival_temperature`) and the temperature at which half the
hydrogen is atomic (`h2_dissociation_temperature`), each at a pressure and a metallicity, which are the edges of the hot-Jupiter
regimes (rows h2 and h3), and the abundances behind the crossover of silicate against haze clouds (row n6), the halogen budget
(rows m5, m6) and the trace gases (row m2).

Inputs are the element amounts (`element_budget` from a metallicity, a carbon-to-oxygen ratio and the solar scale of
`tuning.STELLAR_SOLAR_LOG_EPS`, which is the reference scale and not a preset). The solver reads no random number.

Use:

    eq = solve(1200.0, 1.0, element_budget(metallicity=1.0, c_to_o=0.55))
    eq.fraction("CO"), eq.fraction("CH4")
    ch4_co_temperature(0.1, element_budget(10.0, 0.55))        # K, where CH4 and CO are equally abundant
"""

import json
import math
import os
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import elements
from planetgen.physics.elements import DataError
from planetgen.util import detmath

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "thermo_nasa7.json")
"""str: The data file."""

TOLERANCE = 1e-11
"""float: The relative element balance the iteration stops at."""

MAX_ITERATIONS = 200

HELIUM_PER_HYDROGEN = 0.0851
"""float: Number of helium atoms per hydrogen atom, A21 photosphere (log eps 10.914 - 12); the primordial value differs by under 3 percent."""

_RANGE_CLAMP = 600.0


class EquilibriumError(RuntimeError):
    """The iteration did not converge, or a species or element asked for is not in the table."""


@dataclass(frozen=True)
class Species:
    """
    One gas species and its thermodynamic polynomials.

    Attributes:
        name (str): The NASA name.
        composition (dict): Atoms of each element in the molecule.
        t_low, t_mid, t_high (float): The temperature range of the polynomials and where they join, K.
        low, high (tuple): The seven coefficients below and above `t_mid`.
        note (str): The source note of the fit.
    """

    name: str
    composition: dict
    t_low: float
    t_mid: float
    t_high: float
    low: tuple
    high: tuple
    note: str

    def gibbs_rt(self, temperature):
        """The standard-state (1 bar) Gibbs energy over RT of one mole: H/RT - S/R."""
        a = self.low if temperature < self.t_mid else self.high
        t = temperature
        h_rt = a[0] + a[1] * t / 2 + a[2] * t ** 2 / 3 + a[3] * t ** 3 / 4 + a[4] * t ** 4 / 5 + a[5] / t
        s_r = a[0] * math.log(t) + a[1] * t + a[2] * t ** 2 / 2 + a[3] * t ** 3 / 3 + a[4] * t ** 4 / 4 + a[6]
        return h_rt - s_r


def load(path=DATA_FILE):
    """Reads and checks the thermodynamic table at `path`: a dict of name -> `Species`."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the thermodynamic table {path}: {error}") from error
    table = {}
    for raw in document.get("species") or []:
        owner = f"species {raw.get('name')!r}"
        for key in ("low", "high"):
            values = raw.get(key)
            if not isinstance(values, list) or len(values) != 7 or not all(
                    isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in values):
                raise DataError(f"{owner}: {key} is not seven finite coefficients")
        if not 0 < raw["t_low"] < raw["t_mid"] <= raw["t_high"]:
            raise DataError(f"{owner}: the temperature range is not increasing")
        composition = raw.get("composition")
        if not isinstance(composition, dict) or not composition or not all(
                isinstance(n, int) and n > 0 for n in composition.values()):
            raise DataError(f"{owner}: the composition is missing")
        if raw["name"] in table:
            raise DataError(f"{owner}: listed twice")
        table[raw["name"]] = Species(raw["name"], dict(composition), float(raw["t_low"]), float(raw["t_mid"]), float(raw["t_high"]),
                                     tuple(float(v) for v in raw["low"]), tuple(float(v) for v in raw["high"]), raw.get("note", ""))
    if not table:
        raise DataError("the thermodynamic table is empty")
    return table


SPECIES = load()
"""dict: Species name -> `Species`."""

DEFAULT_SPECIES = tuple(name for name in SPECIES if not set(SPECIES[name].composition) & {"Ti", "V", "Si", "P", "Na", "K", "F", "Cl"})
"""tuple: The species solved for by default: hydrogen, helium, argon, carbon, nitrogen, oxygen and sulphur compounds."""


@dataclass(frozen=True)
class Equilibrium:
    """
    A solved equilibrium.

    Attributes:
        temperature_k, pressure_bar (float): The conditions.
        moles (dict): Moles of each species per the element budget's units (positive ones only).
        total (float): Total moles.
        iterations (int): Newton steps used.
    """

    temperature_k: float
    pressure_bar: float
    moles: dict
    total: float
    iterations: int

    def fraction(self, name):
        """The mole fraction of species `name`."""
        return self.moles.get(name, 0.0) / self.total

    def ratio(self, numerator, denominator):
        """The number ratio of two species."""
        return self.moles.get(numerator, 0.0) / self.moles[denominator]

    def element_in(self, symbol, name):
        """The share of element `symbol` held by species `name`."""
        total = detmath.fsum(n * SPECIES[s].composition.get(symbol, 0) for s, n in self.moles.items())
        return self.moles.get(name, 0.0) * SPECIES[name].composition.get(symbol, 0) / total


def element_budget(metallicity=1.0, c_to_o=None, extra=None):
    """
    The atoms of each element per 1 hydrogen atom for a gas of `metallicity` times the solar metals (A21 log eps,
    `tuning.STELLAR_SOLAR_LOG_EPS`) with helium at the solar ratio. `c_to_o` sets carbon to that multiple of the oxygen
    (None keeps the solar carbon, scaled with the other metals). `extra` adds or overrides elements.
    """
    if not metallicity > 0:
        raise EquilibriumError(f"the metallicity {metallicity!r} is not positive")
    budget = {"H": 1.0, "He": HELIUM_PER_HYDROGEN}
    for symbol in ("C", "N", "O", "S", "Ar"):
        budget[symbol] = metallicity * 10.0 ** (tuning.STELLAR_SOLAR_LOG_EPS[symbol] - 12.0)
    if c_to_o is not None:
        if not c_to_o > 0:
            raise EquilibriumError(f"the C/O {c_to_o!r} is not positive")
        budget["C"] = c_to_o * budget["O"]
    budget.update(extra or {})
    return budget


def _solve_linear(matrix, rhs):
    n = len(rhs)
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-300:
            raise EquilibriumError("the element-potential system is singular")
        a[col], a[pivot] = a[pivot], a[col]
        for r in range(col + 1, n):
            f = a[r][col] / a[col][col]
            if f:
                for c in range(col, n + 1):
                    a[r][c] -= f * a[col][c]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (a[r][n] - detmath.fsum(a[r][c] * x[c] for c in range(r + 1, n))) / a[r][r]
    return x


def _residual(lam, y, atoms, g, ln_p, b, names):
    """(x, E vector): the species fractions and the residuals of the sum-to-one and element equations."""
    ne = len(b)
    x = []
    for i, row in enumerate(atoms):
        e = detmath.fsum(row[j] * lam[j] for j in range(ne)) - g[i] - ln_p
        x.append(math.exp(max(-_RANGE_CLAMP, min(_RANGE_CLAMP, e))))
    res = [detmath.fsum(x) - 1.0]
    for j in range(ne):
        a_j = detmath.fsum(atoms[i][j] * x[i] for i in range(len(x)))
        res.append(y + math.log(a_j) - math.log(b[j]) if a_j > 0 else 1e3)
    return x, res


def solve(temperature_k, pressure_bar, budget, species=None):
    """
    The equilibrium composition of the gas.

    Args:
        temperature_k, pressure_bar (float): The conditions (the polynomials hold from 200 to 6000 K).
        budget (dict): Atoms of each element (only the ratios matter).
        species (iterable of str, optional): The species allowed (default `DEFAULT_SPECIES` limited to the budget's elements).

    Returns:
        Equilibrium: The solution.

    Raises:
        EquilibriumError: A species is unknown, an element has no species, or the iteration does not converge.
    """
    if not temperature_k > 0 or not pressure_bar > 0:
        raise EquilibriumError("the temperature and pressure must be positive")
    names = list(species if species is not None else DEFAULT_SPECIES)
    for name in names:
        if name not in SPECIES:
            raise EquilibriumError(f"unknown species {name}")
    present = sorted(e for e, v in budget.items() if v > 0)
    names = [n for n in names if set(SPECIES[n].composition) <= set(present)]
    for e in present:
        if not any(e in SPECIES[n].composition for n in names):
            raise EquilibriumError(f"no species carries the element {e}")
    atoms = [[SPECIES[n].composition.get(e, 0) for e in present] for n in names]
    b = [budget[e] for e in present]
    g = [SPECIES[n].gibbs_rt(temperature_k) for n in names]
    ln_p = math.log(pressure_bar)
    ne = len(present)
    total_b = detmath.fsum(b)
    # first guess: the element potentials that hold each element in the species that uses the least of it, hydrogen first
    lam = [0.0] * ne
    order = sorted(range(ne), key=lambda j: (present[j] != "H", present[j] not in ("He", "Ar"), -b[j]))
    assigned = set()
    for j in order:
        best = None
        for i, n in enumerate(names):
            composition = SPECIES[n].composition
            if present[j] not in composition or not set(composition) <= {present[k] for k in assigned} | {present[j]}:
                continue
            others = detmath.fsum(composition[present[k]] * lam[k] for k in assigned if present[k] in composition)
            candidate = (g[i] + ln_p + math.log(max(b[j] / total_b, 1e-30)) - others) / composition[present[j]]
            if best is None or candidate < best:
                best = candidate
        lam[j] = best
        assigned.add(j)
    y = math.log(total_b)
    for iteration in range(1, MAX_ITERATIONS + 1):
        x, res = _residual(lam, y, atoms, g, ln_p, b, names)
        norm = max(abs(r) for r in res[1:] + [res[0]])
        if norm < TOLERANCE:
            moles = {n: math.exp(y) * xi for n, xi in zip(names, x) if xi > 0}
            return Equilibrium(temperature_k, pressure_bar, moles, detmath.fsum(moles.values()), iteration)
        a_j = [detmath.fsum(atoms[i][j] * x[i] for i in range(len(x))) for j in range(ne)]
        jac = [[0.0] * (ne + 1) for _ in range(ne + 1)]
        for k in range(ne):
            jac[0][k] = a_j[k]
        for j in range(ne):
            for k in range(ne):
                jac[j + 1][k] = detmath.fsum(atoms[i][j] * atoms[i][k] * x[i] for i in range(len(x))) / a_j[j]
            jac[j + 1][ne] = 1.0
        step = _solve_linear(jac, [-r for r in res])
        scale = min(1.0, 5.0 / max(abs(s) for s in step))
        trial = norm
        for _ in range(40):
            lam_t = [lam[j] + scale * step[j] for j in range(ne)]
            y_t = y + scale * step[ne]
            _x, res_t = _residual(lam_t, y_t, atoms, g, ln_p, b, names)
            trial = max(abs(r) for r in res_t)
            if trial < norm or scale < 1e-6:
                break
            scale /= 2.0
        lam, y = lam_t, y_t
    raise EquilibriumError(f"no convergence at {temperature_k} K and {pressure_bar} bar after {MAX_ITERATIONS} steps")


def _bisect(function, low, high, budget, pressure_bar, tolerance_k=0.5):
    """The temperature in [low, high] where `function(equilibrium)` crosses zero (it must be negative at `low`, positive at `high`)."""
    f_low = function(solve(low, pressure_bar, budget))
    f_high = function(solve(high, pressure_bar, budget))
    if not f_low < 0 < f_high:
        raise EquilibriumError(f"no crossing between {low} and {high} K at {pressure_bar} bar")
    while high - low > tolerance_k:
        mid = 0.5 * (low + high)
        if function(solve(mid, pressure_bar, budget)) < 0:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def ch4_co_temperature(pressure_bar, budget):
    """The temperature, K, where carbon monoxide and methane are equally abundant (above it CO holds the carbon)."""
    return _bisect(lambda e: math.log(e.moles.get("CO", 1e-300) / e.moles.get("CH4", 1e-300)), 300.0, 4000.0, budget, pressure_bar)


def h2o_survival_temperature(pressure_bar, budget, share=0.5):
    """
    The temperature, K, where the water has fallen to `share` of the oxygen left once carbon monoxide has taken its carbon
    (O - C atoms), i.e. where thermal dissociation has destroyed `1 - share` of the water that could exist.
    """
    available = budget["O"] - budget["C"]
    if not available > 0:
        raise EquilibriumError("no oxygen is left for water once carbon monoxide has taken its share (C/O of 1 or more)")
    return _bisect(lambda e: share - e.moles.get("H2O", 0.0) / available, 1500.0, 6000.0, budget, pressure_bar)


def h2_dissociation_temperature(pressure_bar, budget):
    """The temperature, K, where half of the hydrogen nuclei in H2 and H are atomic (`n(H) = 2 n(H2)`)."""
    return _bisect(lambda e: e.moles.get("H", 0.0) - 2.0 * e.moles.get("H2", 0.0), 1000.0, 6000.0, budget, pressure_bar)
