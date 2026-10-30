# planetgen/physics/shoreline_caps.py

"""
The pressure an air can keep, by how far below the cosmic shoreline its body lies (GEN.483, row d6 of the atmosphere recipes).

The shoreline of row d1 is the one of Berta-Thompson, Wachiraphan and Murray 2025 (arXiv:2507.02136): a body keeps an atmosphere
with the probability `1 / (1 + 10^(delta / w))`, where `delta = log10 S - (log10 f0 + p log10 (v_esc / 11.186) + q log10 L)` is
the decimal distance of its bolometric flux `S` (Earth = 1) above the line. They call a body "atmosphere-bearing" when its
surface pressure exceeds 1e-6 bar (all the Solar System planets but Mercury, Titan) or when it holds volatile reservoirs (Triton,
Pluto, Makemake, Eris). The line says whether an air exists; how thick it may be, deep below the line, is a second question the
paper does not answer, and the bracketing caps of row d6 (1e-5, 0.01 and 0.1 bar) were a design.

The measured bodies of the Solar System that the paper labels as having an atmosphere and that have a measured surface pressure
(`BODIES`) give an upper envelope of pressure against the depth below the line: Pluto and Triton hold about 1e-5 bar within 0.4 dex of
it, Titan 1.5 bar at 0.9 dex, Venus 92 bar at 2.1 dex (Mars, with 6 mbar at 0.9 dex, and the Earth are under it). `pressure_cap_bar`
reads the envelope between those record points, linearly in the logarithm of the pressure, and `window_caps` takes it at the deep edge of
each depth window of the recipe. A cap is what the shoreline lets an air reach and not what a body has: the recipe still takes the
smaller of the cap, the supply and the vapour pressure.

Limits. Seven bodies, three of them icy, give four record points: the envelope is a stand-in of tier "educated guess from physics"
and level L2 (registry row SH2), and a rocky body near the line could hold more than the nearby icy ones do. The surface pressures of
Titan and Triton are second-hand textbook values.
"""

import bisect
import math
from dataclasses import dataclass

LOG_F0, P_V, Q_L = 2.62, 5.9, 1.17
"""float: The shoreline of Berta-Thompson, Wachiraphan and Murray 2025: log10 f0 at Earth's escape velocity and a solar luminosity, the
escape-velocity index and the luminosity index (the parameters of row d1, written in `scripts/atmosphere_recipes_data.py`). [Stand-in
SH1, tier T1, tweak L0.]"""
EARTH_ESCAPE_KMS = 11.186


@dataclass(frozen=True)
class Body:
    """A body of the Solar System with a measured surface pressure: its distance from the Sun (AU), escape velocity (km/s) and pressure (bar)."""

    name: str
    distance_au: float
    escape_kms: float
    pressure_bar: float
    source: str


BODIES = (
    Body("Pluto", 39.48, 1.21, 1.3e-5, "NASA Planetary Fact Sheet (about 13 microbar)"),
    Body("Triton", 30.07, 1.455, 1.4e-5, "Voyager 2 occultation, 14 microbar (2ND-HAND, textbook value)"),
    Body("Titan", 9.58, 2.639, 1.467, "Huygens HASI, Fulchignoni et al. 2005 (2ND-HAND, textbook value)"),
    Body("Mars", 1.5237, 5.03, 6.36e-3, "NASA Planetary Fact Sheet (4.0 to 8.7 mbar over the year)"),
    Body("Venus", 0.7233, 10.36, 92.0, "NASA Planetary Fact Sheet"),
    Body("Earth", 1.0, 11.186, 1.014, "NASA Planetary Fact Sheet"),
)
"""tuple: The Solar System bodies of Berta-Thompson et al.'s atmosphere sample with a measured surface pressure. The shoreline is for the
present flux of the Sun (luminosity 1). [Stand-in SH2 reads their pressures; tier T2, tweak L2.]"""


def delta(distance_au, escape_kms, luminosity_solar=1.0):
    """The decimal distance (dex) of a body's flux above the shoreline: positive means above it, with an atmosphere unlikely."""
    flux = luminosity_solar / distance_au ** 2
    return math.log10(flux) - (LOG_F0 + P_V * math.log10(escape_kms / EARTH_ESCAPE_KMS) + Q_L * math.log10(luminosity_solar))


def depth(body):
    """How far (dex) below the shoreline `body` lies (`-delta`)."""
    return -delta(body.distance_au, body.escape_kms)


def envelope(bodies=BODIES):
    """
    The record points of pressure against depth: scanning from the shoreline downward, each body whose pressure is above all those nearer
    the line.

    Returns:
        list: `(depth, log10 pressure_bar)`, shallowest first, the log pressure rising.
    """
    points, highest = [], -math.inf
    for body in sorted(bodies, key=depth):
        level = math.log10(body.pressure_bar)
        if level > highest:
            points.append((depth(body), level))
            highest = level
    return points


def pressure_cap_bar(depth_dex, points=None):
    """
    The most pressure an air may keep at `depth_dex` below the shoreline: the envelope interpolated linearly in log pressure, and the
    first point's pressure from the line down to it. None below the last record point (no cap).
    """
    points = envelope() if points is None else points
    if depth_dex > points[-1][0]:
        return None
    if depth_dex <= points[0][0]:
        return 10.0 ** points[0][1]
    depths = [d for d, _ in points]
    upper = bisect.bisect_left(depths, depth_dex)
    (d0, l0), (d1, l1) = points[upper - 1], points[upper]
    return 10.0 ** (l0 + (l1 - l0) * (depth_dex - d0) / (d1 - d0))


def window_caps(edges, points=None):
    """
    The cap of each depth window of the recipe, taken at its deep edge: `edges` are the depths (dex, rising) at which the windows end,
    the last window being open below its edge. Returns one cap (bar, two significant figures) per edge, then None for the open window.
    """
    caps = []
    for edge in edges:
        cap = pressure_cap_bar(edge, points)
        caps.append(None if cap is None else float(f"{cap:.2g}"))
    caps.append(None)
    return caps
