# planetgen/physics/proton_range.py

"""
Proton range in a column of gas or rock (GEN.453, docs/design/radiation-model.md section 5)
============================================================================================

The kinetic energy a vertical proton needs to cross `X` g/cm2 of a material, from the CSDA range and total stopping power
of NIST PSTAR (ICRU 49). The data are `physics/data/proton_range.json`, built by `scripts/build_proton_range.py`: the
range and the stopping power at 26 energies from 1 MeV to 10 GeV for 11 elements and 6 compounds (the table's `materials`).
This replaces the power law `R = 0.0022 E^1.77`, which holds below 300 MeV and is 5 times too long at 10 GeV.

A CSDA range leaves out nuclear interactions (an interaction length of about 90 g/cm2 in air), so the energy is the
threshold of the primary proton, not of the dose behind it: above about 100 g/cm2 the dose is cascade secondaries.

Use:

    energy_mev(22.5)                              # air: 171 MeV
    energy_mev(21.0, "CO2")                       # Mars
    energy_mev(22.5, {"N": 0.7, "O": 0.3})        # a mixture: Bragg sum of the mass stopping powers
    rigidity_gv(energy_mev(100.0))                # the proton's rigidity, GV
    csda_range_g_cm2(300.0, "air")                # the range of a 300 MeV proton, g/cm2

Between the table's energies the range and the stopping power follow a power law (a straight line in log-log); beyond
10 GeV the last segment's power law continues (the column of Venus or Titan, 10 000 g/cm2, reaches 10 GeV). A material
that is not a row is an error (`DataError`): name its elements and give mass fractions.
"""

import json
import math
import os

from planetgen.physics import elements
from planetgen.physics.elements import DataError

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "proton_range.json")
"""str: The data file."""
PROTON_REST_MEV = 938.272
"""float: The proton's rest energy, MeV."""


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    energies = tuple(data["energies_mev"])
    if any(b <= a for a, b in zip(energies, energies[1:])):
        raise DataError("proton_range.json: the energies must rise")
    tables = {}
    for key, row in data["materials"].items():
        ranges, stopping = tuple(row["csda_range_g_cm2"]), tuple(row["stopping_mev_cm2_g"])
        if len(ranges) != len(energies) or len(stopping) != len(energies):
            raise DataError("proton_range.json: %s does not have a value at every energy" % key)
        if any(b <= a for a, b in zip(ranges, ranges[1:])) or min(ranges) <= 0 or min(stopping) <= 0:
            raise DataError("proton_range.json: %s needs rising ranges and positive stopping powers" % key)
        tables[key] = (ranges, stopping)
    return data, energies, tables


DATA, ENERGIES_MEV, TABLES = _load()
"""The file's contents, its energies (MeV) and `{material: (ranges g/cm2, stopping MeV cm2/g)}`."""
MATERIALS = tuple(TABLES)
"""tuple: The materials of the table."""


def _table(material):
    try:
        return TABLES[material]
    except KeyError:
        raise DataError("no PSTAR table for %r; the materials are %s" % (material, ", ".join(MATERIALS))) from None


def _interpolate(x, xs, ys):
    """`y` at `x` along the straight lines of log y against log x through `(xs, ys)`; the first and last segments
    continue past the ends."""
    if x <= xs[0]:
        low = 0
    elif x >= xs[-1]:
        low = len(xs) - 2
    else:
        low = 0
        high = len(xs) - 1
        while high - low > 1:
            middle = (low + high) // 2
            if xs[middle] <= x:
                low = middle
            else:
                high = middle
    slope = math.log(ys[low + 1] / ys[low]) / math.log(xs[low + 1] / xs[low])
    return ys[low] * (x / xs[low]) ** slope


def csda_range_g_cm2(energy_mev, material="air"):
    """The CSDA range, g/cm2, of a proton of `energy_mev` in `material` (a row of the table, or a mixture)."""
    if energy_mev <= 0:
        return 0.0
    if isinstance(material, dict):
        return _mixture_range(energy_mev, material)
    ranges, _ = _table(material)
    return _interpolate(energy_mev, ENERGIES_MEV, ranges)


def energy_mev(column_g_cm2, material="air"):
    """The kinetic energy, MeV, a vertical proton needs to cross `column_g_cm2` of `material`: a row of the table
    (`"air"`, `"CO2"`, `"H"`, ...) or a mixture `{row: mass fraction}`, whose fractions are normalised."""
    if column_g_cm2 <= 0:
        return 0.0
    if isinstance(material, dict):
        nodes = _mixture_ranges(material)
        return _interpolate(column_g_cm2, nodes, ENERGIES_MEV)
    ranges, _ = _table(material)
    return _interpolate(column_g_cm2, ranges, ENERGIES_MEV)


def rigidity_gv(kinetic_mev):
    """The rigidity, GV, of a proton of `kinetic_mev` kinetic energy."""
    if kinetic_mev <= 0:
        return 0.0
    return math.sqrt(kinetic_mev ** 2 + 2.0 * kinetic_mev * PROTON_REST_MEV) / 1000.0


def _mixture_ranges(fractions):
    """The CSDA range at each table energy of a mixture of mass fractions: Bragg's rule adds the mass stopping powers
    (S = sum w_i S_i), and the range is the integral of dT / S from the first energy on, a power law between nodes, with
    the range at the first energy taken as the mass-weighted harmonic mean of the rows' (the same rule at 1 MeV)."""
    total = sum(fractions.values())
    if total <= 0 or any(w < 0 for w in fractions.values()):
        raise DataError("a mixture needs positive mass fractions")
    weights = {key: w / total for key, w in fractions.items() if w > 0}
    stopping = [sum(w * _table(key)[1][i] for key, w in weights.items()) for i in range(len(ENERGIES_MEV))]
    first = 1.0 / sum(w / _table(key)[0][0] for key, w in weights.items())
    nodes = [first]
    for i in range(len(ENERGIES_MEV) - 1):
        t0, t1 = ENERGIES_MEV[i], ENERGIES_MEV[i + 1]
        s0, s1 = stopping[i], stopping[i + 1]
        p = math.log(s1 / s0) / math.log(t1 / t0)
        # S = s0 (T/t0)^p, so the integral of dT / S is t0 / s0 times the integral of u^-p from 1 to t1/t0.
        u = t1 / t0
        step = t0 / s0 * (math.log(u) if abs(p - 1.0) < 1e-9 else (u ** (1.0 - p) - 1.0) / (1.0 - p))
        nodes.append(nodes[-1] + step)
    return nodes


def _mixture_range(energy, fractions):
    return _interpolate(energy, ENERGIES_MEV, _mixture_ranges(fractions))
