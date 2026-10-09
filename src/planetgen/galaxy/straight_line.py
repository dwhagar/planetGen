# planetgen/galaxy/straight_line.py

"""
Where an object that is not on a galactic orbit is at a given time (GEN.137).

The rule, for course routing (NAV.47) and the orbit update alike: most
objects are bound to the galaxy and are carried round it by
`db.store.advance_galactic_positions` (their position turns about the
galactic axis, their velocity with it). A hypervelocity star is not bound.
It was ejected from the center and keeps going, so its position at time `t`
is

    p(t) = p0 + v * (t - t0)

with `p0` its stored position, `v` its stored galactic velocity (constant)
and `t0` the time `p0` holds at: its `epoch_unix` (a star system's, set by
every orbit update), or, for an object still waiting in `phenomenon_scatter`,
that row's `epoch_unix`, the orbit epoch the scatter was drawn at. A row with
no epoch holds at the database's own orbit epoch.
"""

from planetgen.physics import constants as physical_constants

KM_PER_PC = physical_constants.PARSEC_M / 1000.0
"""float: Kilometers in a parsec."""


def straight_line_position_pc(position_pc, velocity_kms, elapsed_seconds):
    """
    `position_pc` moved by `velocity_kms` for `elapsed_seconds`.

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        velocity_kms (tuple): `(vx, vy, vz)`, km/s in the same frame.
        elapsed_seconds (float): Seconds since the position held (negative
            looks back).

    Returns:
        tuple: The `(x, y, z)` position, parsecs.
    """
    factor = elapsed_seconds / KM_PER_PC
    return tuple(p + v * factor for p, v in zip(position_pc, velocity_kms))


def scatter_position_mpc(row, at_unix):
    """
    A `phenomenon_scatter` row's position (milliparsecs) at `at_unix`: the
    stored point for anything but a hypervelocity star, else the straight
    line from its `epoch_unix`; the stored point when either time is
    unknown.
    """
    point = (row["position_x_mpc"], row["position_y_mpc"], row["position_z_mpc"])
    velocity = (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])
    if row["kind"] != "hypervelocity-star" or None in velocity or row["epoch_unix"] is None or at_unix is None:
        return point
    moved = straight_line_position_pc(tuple(c / 1000.0 for c in point), velocity, at_unix - row["epoch_unix"])
    return tuple(c * 1000.0 for c in moved)
