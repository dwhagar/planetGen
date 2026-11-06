# planetgen/web/lib/cluster_profile.py

"""
MAP.381: how bright a globular cluster looks at each distance from its middle, for the Sector Map's ball.

The ball is a picture of a sphere seen from outside, so what it shows is the *column* density: the King density
(`physics.king.density_per_pc3`) added up along the line of sight at each projected distance from the center. The profile is
returned relative to the center (1 there, 0 at the tidal radius), in equal steps of the tidal radius, so the picture can use it
as it is.
"""

import math

from planetgen.physics import king

PROFILE_STEPS = 8
"""int: Intervals between the center and the tidal radius; the profile has one more value than this."""

LINE_OF_SIGHT_STEPS = 32
"""int: Steps of the sum along one line of sight (midpoint rule)."""


def surface_profile(component, steps=PROFILE_STEPS):
    """
    The projected density of a King cluster at `steps + 1` equal steps from its center to its tidal radius, over the center's.

    Args:
        component (dict): The cluster's `component` as the container API gives it: `mass_solar`, `half_mass_radius_pc`,
            `tidal_radius_pc` and `king_w0`.
        steps (int): Intervals between the center and the tidal radius.

    Returns:
        list[float] | None: `steps + 1` values from 1.0 falling to 0.0, or `None` when the component lacks what the King model needs.
    """
    try:
        w0 = int(component["king_w0"])
        mass = float(component["mass_solar"])
        half_mass = float(component["half_mass_radius_pc"])
        tidal = float(component["tidal_radius_pc"])
    except (KeyError, TypeError, ValueError):
        return None
    if not (mass > 0.0 and half_mass > 0.0 and tidal > 0.0) or steps < 1:
        return None
    try:
        king.king_table(w0)
    except ValueError:
        return None
    edge = min(tidal, king.tidal_radius_pc(w0, half_mass))     # the density is zero beyond the King model's own edge
    columns = []
    for i in range(steps + 1):
        projected = tidal * i / steps
        if projected >= edge:
            columns.append(0.0)
            continue
        reach = math.sqrt(edge * edge - projected * projected)
        dz = reach / LINE_OF_SIGHT_STEPS
        columns.append(sum(king.density_per_pc3(w0, mass, half_mass, math.hypot(projected, (k + 0.5) * dz))
                           for k in range(LINE_OF_SIGHT_STEPS)) * dz)
    if columns[0] <= 0.0:
        return None
    return [round(value / columns[0], 4) for value in columns]
