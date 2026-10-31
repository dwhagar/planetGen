# planetgen/physics/gas_equilibrium.py

"""
Ideal-gas chemical equilibrium at a given temperature and pressure (GEN.451, GEN.454, GEN.499)
==============================================================================================

The one equilibrium solver of the atmosphere chain: the mole fractions of every gas of `physics/species_thermo.py` that can
be made of the given elements, at the temperature and total pressure of the layer. It minimises the Gibbs energy of an ideal gas
at constant temperature and pressure under the element balance by the method of Gordon and McBride (1994, NASA RP-1311,
section 2): Newton's method on the element potentials, with the step damped so that no species changes by more than a factor
of e^2 per iteration and a negligible species stays negligible. Two things beyond their problem are added: an element whose
potential is held (the oxygen of a mantle whose oxygen fugacity is given, `fixed_o2_log_fugacity`), and a result that carries
the mole fractions of every species, the ones below 1e-300 included as zero.

    result = solve({"H": 1.0, "C": 0.1, "N": 0.02, "S": 0.02}, 1200.0, 1.0, fixed_o2_log_fugacity=-12.2)
    result.mole_fractions["H2O"]        # 0.30
    result.moles["H2"]                  # moles per the abundance scale given

Abundances are moles of atoms in any common scale (only their ratios and the pressure decide the mole fractions). With the oxygen
potential held the oxygen abundance is not asked for: the mantle supplies or takes the oxygen that the buffer requires, so the
mole fractions are those of a gas in contact with it, and the moles of the gas are per the other abundances. Condensed
phases are not in the solution: the caller caps what the saturation curves allow (`physics/saturation.py`), because the vapour
pressure of an ice is a pure-phase curve and the freeze-out belongs to the pressure bookkeeping of the generator (GEN.422).

Species whose fit does not reach the temperature are left out and named in `skipped`. The standard state is 1 bar, as the
fits' (`species_thermo`); pressures are in bar.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from planetgen.physics import species_thermo
from planetgen.util import detmath

NEGLIGIBLE_LOG_FRACTION = -18.420680743952367
"""float: ln 1e-8, the mole fraction below which a species is treated as negligible when the step is damped (Gordon and McBride)."""

BALANCE_TOLERANCE = 1e-8
"""float: The relative element-balance error the converged solution may have."""

CONVERGED_LOG_STEP = 1e-6
"""float: The largest change of ln n_j, over the species above `CONVERGED_FRACTION`, that counts as converged. The step has been
applied when it is tested, and Newton's error after a step is the square of the step, so the answer is good to 1e-12; a smaller
test is out of reach for a trace element, whose residual carries the rounding of the major ones."""

STALL_ITERATIONS = 60
"""int: A solve whose smallest step has not improved for this many iterations, and is under `STALL_LOG_STEP`, has converged to the
precision the rounding allows."""

STALL_LOG_STEP = 1e-4

CONVERGED_FRACTION = 1e-12
"""float: Species rarer than this are not waited for: the trace species of a stiff system jitter at 1e-8 in the Newton step while
the rest has converged to 1e-12, and they are recomputed exactly from the converged element potentials at the end."""

TRACE_ELEMENT_FRACTION = 1e-12
"""float: An element rarer than this against the commonest one is left out of a second try when the first does not settle: the
Newton step of a gas that is all one element (a nitrogen air with 1e-15 of everything else) can stall with the rarest elements
in, and what they would add is below the precision of the rest."""

MOLES_FLOOR = 1e-290
"""float: The least moles a species is allowed, so that no logarithm is minus infinity."""

MAX_ITERATIONS = 400


def _times_vector(matrix, vector):
    """`matrix @ vector` without BLAS (a result that goes to the database must not depend on the library)."""
    return detmath.matmul(matrix, np.asarray(vector)[:, None])[:, 0]


class EquilibriumError(RuntimeError):
    """The solver did not converge, or the elements cannot make a gas (an element no species holds)."""


@dataclass(frozen=True)
class Equilibrium:
    """
    The solution.

    Attributes:
        temperature_k (float): The temperature, K.
        pressure_bar (float): The total pressure, bar.
        moles (dict): Species name -> moles, in the scale of the abundances given (species below 1e-300 of the total are 0).
        mole_fractions (dict): Species name -> mole fraction, summing to 1.
        total_moles (float): The moles of gas.
        iterations (int): How many Newton iterations it took.
        skipped (tuple): The species left out because their fit does not reach `temperature_k`.
    """

    temperature_k: float
    pressure_bar: float
    moles: dict
    mole_fractions: dict
    total_moles: float
    iterations: int
    skipped: tuple = field(default=())


def solve(abundances, temperature_k, pressure_bar, *, fixed_o2_log_fugacity=None, species=None):
    """
    The equilibrium gas of `abundances`, tried again without the trace elements (`TRACE_ELEMENT_FRACTION`) when it does not settle
    with them in.  See `_solve` for the arguments.
    """
    try:
        return _solve(abundances, temperature_k, pressure_bar, fixed_o2_log_fugacity=fixed_o2_log_fugacity, species=species)
    except EquilibriumError:
        given = {symbol: float(value) for symbol, value in abundances.items() if value and value > 0}
        if not given:
            raise
        kept = {symbol: value for symbol, value in given.items() if value >= TRACE_ELEMENT_FRACTION * max(given.values())}
        if len(kept) == len(given):
            raise
        return _solve(kept, temperature_k, pressure_bar, fixed_o2_log_fugacity=fixed_o2_log_fugacity, species=species)


def _solve(abundances, temperature_k, pressure_bar, *, fixed_o2_log_fugacity=None, species=None):
    """
    The equilibrium gas of `abundances` at `temperature_k` and `pressure_bar`.

    Args:
        abundances (dict): Element symbol -> moles of atoms (any common scale, positive; an element of 0 is left out).
        temperature_k (float): Kelvin.
        pressure_bar (float): Total pressure, bar.
        fixed_o2_log_fugacity (float, optional): log10 of the oxygen fugacity (bar) the gas is held at, for gas over a buffered
            mantle. The oxygen is then not an abundance: `abundances` must not carry it.
        species (iterable, optional): The gases to allow, by name in `species_thermo`; default every gas made of the elements.

    Returns:
        Equilibrium: The solution.

    Raises:
        EquilibriumError: No convergence, or an element that no allowed species holds.
        ValueError: A non-positive temperature, pressure or abundance, or oxygen given as an abundance with its fugacity held.
    """
    if not temperature_k > 0 or not pressure_bar > 0:
        raise ValueError("the temperature and the pressure must be positive")
    elements_given = {symbol: float(value) for symbol, value in abundances.items() if value}
    if any(value < 0 for value in elements_given.values()):
        raise ValueError("an abundance must not be negative")
    if not elements_given:
        raise ValueError("no element has an abundance")
    held = fixed_o2_log_fugacity is not None
    if held and "O" in elements_given:
        raise ValueError("oxygen is not an abundance when its fugacity is held")

    free = sorted(elements_given)
    allowed = set(free) | ({"O"} if held else set())
    candidates = species if species is not None else species_thermo.GASES
    chosen, skipped = [], []
    for name in candidates:
        row = species_thermo.SPECIES[name]
        if row.condensed or not set(row.elements) <= allowed:
            continue
        low, high = row.range_k
        if not low <= temperature_k <= high:
            skipped.append(name)
            continue
        chosen.append(name)
    if not chosen:
        raise EquilibriumError("no species can be made of these elements at this temperature")
    row_of = {symbol: index for index, symbol in enumerate(free)}
    matrix = np.zeros((len(free), len(chosen)))
    oxygen_atoms = np.zeros(len(chosen))
    for column, name in enumerate(chosen):
        for symbol, count in species_thermo.SPECIES[name].elements.items():
            if symbol in row_of:
                matrix[row_of[symbol], column] = count
            else:
                oxygen_atoms[column] = count
    for symbol, row in row_of.items():
        if not matrix[row].any():
            raise EquilibriumError(f"no species holds {symbol} at {temperature_k} K")

    scale = sum(elements_given.values())
    totals = np.array([elements_given[symbol] / scale for symbol in free])
    chemical_constant = np.array([species_thermo.SPECIES[name].g_over_rt(temperature_k) for name in chosen]) + math.log(pressure_bar)
    if held:
        # The potential of an oxygen atom is half that of O2 at the held fugacity: (g_O2 + ln fO2) / 2; a species with
        # n oxygen atoms is shifted by n times it.
        potential_of_oxygen = 0.5 * (species_thermo.SPECIES["O2"].g_over_rt(temperature_k) + fixed_o2_log_fugacity * math.log(10.0))
        chemical_constant = chemical_constant - oxygen_atoms * potential_of_oxygen

    count = len(chosen)
    for anchored in (False, True):
        moles = np.full(count, 0.1 / count)
        total = 0.1
        best_step, best_at = math.inf, 0
        for iterations in range(1, MAX_ITERATIONS + 1):
            log_fraction = np.log(moles / total)
            potentials = chemical_constant + log_fraction
            weighted = matrix * moles
            row_sums = weighted.sum(axis=1)
            lhs = np.empty((len(free) + 1, len(free) + 1))
            lhs[:-1, :-1] = detmath.matmul(weighted, matrix.T)
            lhs[:-1, -1] = row_sums
            lhs[-1, :-1] = row_sums
            lhs[-1, -1] = moles.sum() - total
            rhs = np.empty(len(free) + 1)
            rhs[:-1] = totals - row_sums + _times_vector(weighted, potentials)
            rhs[-1] = total - moles.sum() + float(_times_vector(moles[None, :], potentials)[0])
            # The rows of a trace element weigh 1e-9 or less against the main ones: scale the symmetric system by its diagonal so that
            # the solve keeps its digits (the unscaled Newton step jitters at 1e-8 and never meets `CONVERGED_LOG_STEP`).
            diagonal = np.abs(np.diag(lhs))
            if anchored:
                # The last row holds `sum(n) - total`, which is rounding noise (1e-16) once the solve is near its end; scaled by its own
                # diagonal it blows the coupling to a trace element up and the step jitters for ever (GEN.1003: near-stoichiometric water with
                # 1e-12 of argon at 726 K). The retry scales it by the total, its natural size. A solve that converges first time never gets here.
                diagonal[-1] = total
            scale_of = 1.0 / np.sqrt(np.where(diagonal > 0.0, diagonal, 1.0))
            scaled = lhs * scale_of[:, None] * scale_of[None, :]
            try:
                solution = np.linalg.solve(scaled, rhs * scale_of) * scale_of
            except np.linalg.LinAlgError:
                solution = np.linalg.lstsq(scaled, rhs * scale_of, rcond=None)[0] * scale_of
            element_potentials, d_ln_total = solution[:-1], solution[-1]
            d_ln_moles = _times_vector(matrix.T, element_potentials) - potentials + d_ln_total
            significant = log_fraction > NEGLIGIBLE_LOG_FRACTION
            biggest = max(5.0 * abs(d_ln_total), float(np.max(np.abs(d_ln_moles[significant]))) if significant.any() else 0.0)
            factor = 1.0 if biggest <= 2.0 else 2.0 / biggest
            rising = ~significant & (d_ln_moles >= 0.0)
            if rising.any():
                gap = d_ln_moles[rising] - d_ln_total
                gap = gap[gap != 0.0]
                if gap.size:
                    room = np.abs((-log_fraction[rising][d_ln_moles[rising] - d_ln_total != 0.0] - 9.2103) / gap)
                    factor = min(factor, float(room.min()))
            moles = np.maximum(moles * np.exp(factor * d_ln_moles), MOLES_FLOOR)
            total = total * math.exp(factor * d_ln_total)
            live = moles > CONVERGED_FRACTION * total
            step = float(np.max(np.abs(d_ln_moles[live]))) if live.any() else 0.0
            if factor == 1.0 and step < CONVERGED_LOG_STEP and abs(d_ln_total) < CONVERGED_LOG_STEP:
                break
            if factor == 1.0 and step < best_step:
                best_step, best_at = step, iterations
            if iterations - best_at > STALL_ITERATIONS and best_step < STALL_LOG_STEP:
                break
        else:
            if anchored:
                raise EquilibriumError(f"no convergence at {temperature_k} K and {pressure_bar} bar")
            continue
        break
    # one undamped pass on every species from the converged potentials, so the trace species are exact too
    moles = np.maximum(np.exp(math.log(total) - chemical_constant + _times_vector(matrix.T, element_potentials)), 0.0)
    if np.abs(_times_vector(matrix, moles) - totals).max() / totals.max() > BALANCE_TOLERANCE:
        raise EquilibriumError(f"the element balance is off at {temperature_k} K and {pressure_bar} bar")

    fractions = moles / moles.sum()
    return Equilibrium(
        temperature_k=temperature_k,
        pressure_bar=pressure_bar,
        moles={name: float(value) * scale for name, value in zip(chosen, moles)},
        mole_fractions={name: float(value) for name, value in zip(chosen, fractions)},
        total_moles=float(moles.sum()) * scale,
        iterations=iterations,
        skipped=tuple(skipped),
    )
