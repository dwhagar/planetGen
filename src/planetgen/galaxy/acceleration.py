# planetgen/galaxy/acceleration.py

"""
The acceleration kernel of the orbital update (GEN.342; docs/design/orbits/orbital-updates.md section 4)
=========================================================================================================

One function, `pull`, returns the acceleration at a point: the pull of the object's influence set (the pairwise point masses of
`galaxy.influence`, a folded sector as one point at its centre of mass) plus the galactic gradient of `galaxy.potential`, and one,
`kick`, turns it into a new velocity. The object itself is left out of its own set (the caller builds the set with the object's key in
`exclude_keys`, as `influence.influence_set` does), so nothing here needs to know which member is the object.

- **The set's term** is Newtonian, `a = G m d / (|d|^2 + eps^2)^(3/2)`, summed over the exact members and the monopoles. It is unsoftened
  to the tiny `MEMBER_SOFTENING_PC` (0.001 pc, about 200 AU) that keeps the pull finite for two masses at one point; the central black
  hole's own 1 pc softening is the potential's (`potential.BLACK_HOLE_SOFTENING_PC`), where it lives, not in the set (the black hole is left out of
  the sets). There is no minimum-acceleration cut: every member inside the radius counts.
- **A star** takes `pull` as the extra term in the stepper's kick (GEN.340). **A planet or moon** works in its host's frame: `host_frame_pull`
  is the pull at the object less the pull at the host, both with the host left out of the set, so the host's own gravity (the orbit solver's)
  is not counted twice and only the perturbation is left.
- Nothing is moved and nothing is stored: the kernel is a pure function of positions, masses and the potential, and the same inputs give the
  same numbers. A caller's order of objects does not change it (the chunk rule of GEN.105 is the caller's).

Units: positions in parsecs, masses in solar masses (the set's kilograms are converted here), the pull in (km/s)^2 per pc, velocities in km/s
and the step in Myr (`potential.KMS_TO_PC_PER_MYR`).
"""

import numpy as np

from planetgen.galaxy import potential as galactic
from planetgen.physics import constants

MEMBER_SOFTENING_PC = 1.0e-3
"""float: The softening of a member's pull, parsecs: the stand-in that keeps the pull of two masses at one point finite (registry row AK1)."""

_KG_PER_MSUN = constants.SOLAR_MASS_TO_KG


def set_pull(position_pc, influence, softening_pc=MEMBER_SOFTENING_PC):
    """
    The pull of an `influence.InfluenceSet` at `position_pc` (x, y, z, parsecs), in (km/s)^2 per pc: every exact member and every monopole as a point
    mass with the softening `softening_pc`. An empty set pulls with `(0.0, 0.0, 0.0)`.
    """
    here = np.asarray(position_pc, dtype=float)
    positions = np.asarray(influence.positions_pc, dtype=float).reshape(-1, 3)
    masses = np.asarray(influence.masses_kg, dtype=float).reshape(-1)
    if influence.monopoles:
        positions = np.vstack([positions, [m.position_pc for m in influence.monopoles]])
        masses = np.concatenate([masses, [m.mass_kg for m in influence.monopoles]])
    if not len(masses):
        return 0.0, 0.0, 0.0
    offset = positions - here
    distance2 = np.einsum("ij,ij->i", offset, offset) + softening_pc ** 2
    weight = galactic.G_PC_KMS2_PER_MSUN * (masses / _KG_PER_MSUN) / distance2 ** 1.5
    pull = weight @ offset
    return float(pull[0]), float(pull[1]), float(pull[2])


def pull(position_pc, pot, influence=None, softening_pc=MEMBER_SOFTENING_PC):
    """
    The acceleration at `position_pc`, (km/s)^2 per pc: the galactic gradient of `pot` (`potential.Potential`) plus the pull of `influence` (the
    object's `InfluenceSet`, `None` for none). With no set it is the potential's pull alone.
    """
    gx, gy, gz = galactic.acceleration(pot, *position_pc)
    if influence is None:
        return gx, gy, gz
    sx, sy, sz = set_pull(position_pc, influence, softening_pc)
    return gx + sx, gy + sy, gz + sz


def host_frame_pull(position_pc, host_position_pc, pot, influence=None, softening_pc=MEMBER_SOFTENING_PC):
    """
    The perturbing acceleration of a planet or moon at `position_pc` in the frame of its host at `host_position_pc`: `pull` at the object less `pull` at the
    host. The caller leaves the host out of `influence` (its pull is the orbit's own), so what is left is the tide of the galaxy and of the other masses.
    """
    return tuple(a - b for a, b in zip(pull(position_pc, pot, influence, softening_pc), pull(host_position_pc, pot, influence, softening_pc)))


def kick(velocity_kms, acceleration, step_myr):
    """The velocity, km/s, after `step_myr` Myr of the acceleration `acceleration` ((km/s)^2 per pc): `v + a dt`, with `dt` in pc per (km/s) (a Myr is 1.0227 of them)."""
    step = step_myr * galactic.KMS_TO_PC_PER_MYR
    return tuple(v + a * step for v, a in zip(velocity_kms, acceleration))


def velocity_update(position_pc, velocity_kms, pot, influence=None, step_myr=0.0, host_position_pc=None, softening_pc=MEMBER_SOFTENING_PC):
    """
    The velocity after a kick of `step_myr` Myr at `position_pc`: for a star `host_position_pc` is `None` and the whole pull kicks; for a planet or moon
    it is the host's position and only the host-frame perturbation does.
    """
    if host_position_pc is None:
        accel = pull(position_pc, pot, influence, softening_pc)
    else:
        accel = host_frame_pull(position_pc, host_position_pc, pot, influence, softening_pc)
    return kick(velocity_kms, accel, step_myr)
