# planetgen/physics/lagrange_geometry.py

"""
Where the Lagrange points of a pair are, for the System Map (MAP.354, and the gravity layer MAP.169; docs/design/systems/lagrange-trojans.md)
============================================================================================================================================

Pure geometry in the system frame, no drawing and no database. For a host at `r1` (mass `m1`) and a secondary at `r2` (mass `m2`), with the secondary's velocity
relative to the host `v` (the direction of motion; its cross product with the separation gives the orbit's normal), the five points of the circular restricted
problem are (units of the arguments; `a` is the separation, `mu = m2 / (m1 + m2)`, `B` the barycenter, `u` the unit vector host to secondary, `w` the unit vector along the motion):

    L1  between the two, L2 beyond the secondary, L3 beyond the host (the roots of the collinear force balance, found by bisection)
    L4  B + a ((1/2 - mu) u + (sqrt(3)/2) w)     60 degrees ahead of the secondary
    L5  B + a ((1/2 - mu) u - (sqrt(3)/2) w)     60 degrees behind

`points` returns all five. `occupied` keeps those that hold objects and gives each its center and its label, `Jupiter-L4` (the same as `trojans.point_label`, which names only L4 and L5; here L1 to L3 are labelled too for the gravity layer), the icon the
map draws at the center. Only L4 and L5 hold natural bodies (`trojans`), but L1 to L3 are given for the gravity layer.
"""

import math

from planetgen.physics import trojans, vec3

_ALL = ("L1", "L2", "L3", "L4", "L5")


def _unit(a):
    if vec3.norm(a) == 0.0:
        raise trojans.TrojanError("a zero vector has no direction")
    return vec3.unit(a)


def collinear_offsets(mu):
    """`(x1, x2, x3)`: the x-coordinates of L1, L2, L3 in the rotating frame in units of the separation, the barycenter at 0, the host at `-mu`, the secondary at `1 - mu`."""
    if not 0.0 < mu < 0.5:
        raise trojans.TrojanError(f"the mass parameter must be between 0 and 0.5, got {mu!r}")

    def force(x):
        return x - (1.0 - mu) * (x + mu) / abs(x + mu) ** 3 - mu * (x - 1.0 + mu) / abs(x - 1.0 + mu) ** 3

    def root(low, high):
        f_low = force(low)
        for _ in range(200):
            mid = 0.5 * (low + high)
            f_mid = force(mid)
            if (f_mid > 0.0) == (f_low > 0.0):
                low, f_low = mid, f_mid
            else:
                high = mid
        return 0.5 * (low + high)

    eps = 1e-12
    x1 = root(-mu + eps, 1.0 - mu - eps)
    x2 = root(1.0 - mu + eps, 3.0)
    x3 = root(-3.0, -mu - eps)
    return x1, x2, x3


def points(host_pos, host_mass, secondary_pos, secondary_mass, secondary_velocity, host_velocity=(0.0, 0.0, 0.0)):
    """`{"L1": [x, y, z], ..., "L5": [...]}`: the five points in the frame of the positions, for a secondary moving with `secondary_velocity` (relative to the host's `host_velocity`)."""
    separation_vector = vec3.sub(secondary_pos, host_pos)
    a = vec3.norm(separation_vector)
    u = _unit(separation_vector)
    normal = _unit(vec3.cross(separation_vector, vec3.sub(secondary_velocity, host_velocity)))
    w = _unit(vec3.cross(normal, u))
    mu = trojans.mass_parameter(host_mass, secondary_mass)
    total = host_mass + secondary_mass
    barycenter = [(host_mass * host_pos[i] + secondary_mass * secondary_pos[i]) / total for i in range(3)]
    x1, x2, x3 = collinear_offsets(mu)
    out = {}
    for name, x in (("L1", x1), ("L2", x2), ("L3", x3)):
        out[name] = [barycenter[i] + a * x * u[i] for i in range(3)]
    for name, sign in (("L4", 1.0), ("L5", -1.0)):
        out[name] = [barycenter[i] + a * ((0.5 - mu) * u[i] + sign * (math.sqrt(3.0) / 2.0) * w[i]) for i in range(3)]
    return out


def occupied(pair_points, secondary_name, held_points):
    """The points of `held_points` (names such as `"L4"`) that have objects, each as `{"point", "position", "label"}` with the label `Jupiter-L4`."""
    result = []
    for name in _ALL:
        if name in held_points:
            result.append({"point": name, "position": pair_points[name], "label": f"{secondary_name}-{name}"})
    return result
