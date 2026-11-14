# planetgen/physics/vec3.py

"""
The few 3-vector operations on plain tuples the physics modules share (GEN.1667): dot, cross, norm, unit and scale.
They were written four times over (`state_vectors`, `spin`, `sky_suns`, `daylight`); numpy is not used here because
these are called on single vectors, where a tuple is cheaper than an array.
"""

import math


def dot(a, b):
    """The dot product of `a` and `b`."""
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    """The cross product `a x b`."""
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    """The length of `a`."""
    return math.sqrt(dot(a, a))


def unit(a):
    """`a` scaled to length 1."""
    length = norm(a)
    return (a[0] / length, a[1] / length, a[2] / length)


def scale(a, k):
    """`a` times the scalar `k`."""
    return (a[0] * k, a[1] * k, a[2] * k)


def add(a, b):
    """The sum `a + b`."""
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
