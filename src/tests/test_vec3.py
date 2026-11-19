"""GEN.1671: the one set of 3-vector operations (`physics.vec3`) the physics, galaxy and database modules share."""

import math

from planetgen.physics import vec3


def test_the_operations_agree_with_the_plain_formulas():
    a, b = (1.0, 2.0, 3.0), (-4.0, 0.5, 2.0)
    assert vec3.add(a, b) == (-3.0, 2.5, 5.0)
    assert vec3.sub(a, b) == (5.0, 1.5, 1.0)
    assert vec3.dot(a, b) == -4.0 + 1.0 + 6.0
    assert vec3.cross(a, b) == (2.0 * 2.0 - 3.0 * 0.5, 3.0 * -4.0 - 1.0 * 2.0, 1.0 * 0.5 - 2.0 * -4.0)
    assert vec3.scale(a, 2.0) == (2.0, 4.0, 6.0)
    assert vec3.norm((3.0, 4.0, 12.0)) == 13.0
    assert math.isclose(vec3.norm(vec3.unit(a)), 1.0, rel_tol=1e-15)
