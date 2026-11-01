# planetgen/galaxy/potential.py

"""
The galactic potential as one pure module (GEN.338; docs/design/orbits/orbital-solvers-and-integrators.md)
==============================================================================================================

One function of position returns the potential and the pull of a galaxy, and one the circular speed in its plane, from three kinds of mass:

- **Hernquist spheres**, `Phi = -G M / (r + a)`: the bulge (the bar's three axes are not kept: a sphere of the bar's long
  axis), the stellar halo, the nuclear star cluster;
- **Miyamoto-Nagai disks** [R3103], `Phi = -G M / sqrt(R^2 + (a + sqrt(z^2 + b^2))^2)`: the thin and thick disks and the nuclear
  disc, with `a` the disk's scale length and `b` its scale height;
- **an NFW halo** [R3248] of the dark matter, `Phi = -G M_s ln(1 + r / r_s) / r`, of the mass chain's halo less the galaxy's
  stars, its virial radius and its concentration (`enclosed_mass.dark_halo`).

and the central black hole, a Plummer point of 1 pc softening, which is the only softened part. Every pull is finite at the origin.

**No Milky Way value is a default.** `from_shape` takes every mass and scale from the galaxy: the six components' masses and scales from the density model
(`density.density_components`, GEN.379) and the halo from `enclosed_mass.dark_halo`; the circular speed follows from them and is not a constant. A component
of mass 0 (a galaxy without a disk, an elliptical) adds nothing and is legal. The matching of a Miyamoto-Nagai disk to the exponential disk of the density
model (`a` the scale length, `b` the scale height) and of a Hernquist sphere to the bulge, the stellar halo and the nuclear cluster (`a` their scale radius,
break radius and cluster scale) are stand-ins (registry rows GP1 to GP3 of 6.157 in `gap-filling-plan.md`); the three-term fit of the Milky Way's
rotation curve in the literature ([R0496]'s MWPotential2014: a Hernquist-like bulge, a Miyamoto-Nagai disk, an NFW halo) is the same kind of model with
its parameters fitted, not derived.

Units: parsecs, solar masses and km/s. The pull is in (km/s)^2 per pc (`ACCELERATION_TO_SI` turns it into m/s^2); a time in Myr is `KMS_TO_PC_PER_MYR` (1.0227)
pc per (km/s), the unit a pull times a time gives a speed in. The arrays are numpy arrays, so one call evaluates many points, and the same inputs give the same
numbers.
"""

import collections
import dataclasses
import math

import numpy as np

from planetgen.galaxy import density, enclosed_mass
from planetgen.physics import constants

G_PC_KMS2_PER_MSUN = constants.G * constants.SOLAR_MASS_TO_KG / (constants.PARSEC_M * 1.0e6)
"""float: The gravitational constant in pc (km/s)^2 per solar mass (4.3009e-3)."""
ACCELERATION_TO_SI = 1.0e6 / constants.PARSEC_M
"""float: One (km/s)^2 per pc in m/s^2."""
KMS_TO_PC_PER_MYR = 1.0e3 * 365.25 * 86400.0 * 1.0e6 / constants.PARSEC_M
"""float: One km/s in pc per Myr (1.0227)."""
BLACK_HOLE_SOFTENING_PC = 1.0
"""float: The Plummer softening of the central black hole, parsecs; the only softened part of the potential (GEN.338)."""

HERNQUIST_COMPONENTS = {"bulge": "scale_x_pc", "halo": "break_pc", "nuclear_cluster": "scale_pc"}
"""dict: The components of the density model a Hernquist sphere stands for, and the parameter of each that is its scale radius (stand-ins GP2 and GP3)."""
DISK_COMPONENTS = ("thin_disk", "thick_disk", "nuclear_disc")
"""tuple: The components a Miyamoto-Nagai disk stands for, with `a` the component's `scale_length_pc` and `b` its `scale_height_pc` (stand-in GP1)."""

Hernquist = collections.namedtuple("Hernquist", ("name", "mass_solar", "scale_pc"))
Disk = collections.namedtuple("Disk", ("name", "mass_solar", "scale_length_pc", "scale_height_pc"))
Halo = collections.namedtuple("Halo", ("mass_solar", "concentration", "virial_radius_pc"))


@dataclasses.dataclass(frozen=True)
class Potential:
    """A galaxy's potential: Hernquist spheres, Miyamoto-Nagai disks, one NFW halo and the central black hole (`black_hole_solar`)."""

    spheres: tuple = ()
    disks: tuple = ()
    halo: Halo = Halo(0.0, 1.0, 1.0)
    black_hole_solar: float = 0.0
    softening_pc: float = BLACK_HOLE_SOFTENING_PC

    def __post_init__(self):
        for part in (*self.spheres, *self.disks):
            if part.mass_solar < 0.0 or not math.isfinite(part.mass_solar):
                raise ValueError(f"{part.name}: the mass must be zero or positive and finite")
        if any(s.scale_pc <= 0.0 for s in self.spheres) or any(d.scale_length_pc < 0.0 or d.scale_height_pc <= 0.0 for d in self.disks):
            raise ValueError("a scale radius must be positive (a disk's scale length may be 0)")
        if self.halo.mass_solar < 0.0 or self.halo.concentration <= 0.0 or self.halo.virial_radius_pc <= 0.0:
            raise ValueError("the halo needs a mass of zero or more, a concentration and a virial radius above zero")

    @property
    def nfw_scale_pc(self):
        """The NFW scale radius `r_s = r_virial / c`, parsecs."""
        return self.halo.virial_radius_pc / self.halo.concentration

    @property
    def nfw_characteristic_solar(self):
        """`M_s = M_halo / (ln(1 + c) - c / (1 + c))`, solar masses: the halo's mass inside the virial radius is `M_s m(c)`."""
        c = self.halo.concentration
        return self.halo.mass_solar / (math.log1p(c) - c / (1.0 + c))


def _arrays(x, y, z):
    x, y, z = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (x, y, z)))
    return x, y, z


def _nfw_mass_function(u):
    """`m(u) = ln(1 + u) - u / (1 + u)`, with its series `u^2 / 2 - 2 u^3 / 3` near zero where the difference loses its digits."""
    u = np.asarray(u, dtype=float)
    small = u < 1.0e-4
    safe = np.where(small, 1.0, u)
    return np.where(small, u * u / 2.0 - 2.0 * u ** 3 / 3.0, np.log1p(safe) - safe / (1.0 + safe))


def potential(pot, x, y, z):
    """The potential at `(x, y, z)` (parsecs, the galaxy's axes with the centre at the origin), in (km/s)^2, zero at infinity. A float for float input."""
    x, y, z = _arrays(x, y, z)
    r = np.sqrt(x * x + y * y + z * z)
    big_g = G_PC_KMS2_PER_MSUN
    total = np.zeros_like(r)
    for sphere in pot.spheres:
        total -= big_g * sphere.mass_solar / (r + sphere.scale_pc)
    cylinder2 = x * x + y * y
    for disk in pot.disks:
        total -= big_g * disk.mass_solar / np.sqrt(cylinder2 + (disk.scale_length_pc + np.sqrt(z * z + disk.scale_height_pc ** 2)) ** 2)
    if pot.black_hole_solar:
        total -= big_g * pot.black_hole_solar / np.sqrt(r * r + pot.softening_pc ** 2)
    if pot.halo.mass_solar:
        scale = pot.nfw_scale_pc
        u = r / scale
        small = u < 1.0e-8
        # ln(1 + u) / u -> 1 at the centre, so the potential there is -G M_s / r_s.
        ratio = np.where(small, 1.0, np.log1p(np.where(small, 1.0, u)) / np.where(small, 1.0, u))
        total -= big_g * pot.nfw_characteristic_solar / scale * ratio
    return float(total) if total.ndim == 0 else total


def acceleration(pot, x, y, z):
    """
    The pull `(ax, ay, az)` at `(x, y, z)`, parsecs, in (km/s)^2 per pc: minus the gradient of the potential, finite everywhere including the origin (it is 0
    there). Arrays in, arrays out; floats in, floats out.
    """
    x, y, z = _arrays(x, y, z)
    r2 = x * x + y * y + z * z
    r = np.sqrt(r2)
    big_g = G_PC_KMS2_PER_MSUN
    ax, ay, az = np.zeros_like(r), np.zeros_like(r), np.zeros_like(r)
    safe_r = np.where(r > 0.0, r, 1.0)

    def radial(magnitude):
        """Add a central pull `-magnitude * r_hat`, magnitude >= 0, nothing at r = 0."""
        nonlocal ax, ay, az
        factor = np.where(r > 0.0, magnitude / safe_r, 0.0)
        ax, ay, az = ax - factor * x, ay - factor * y, az - factor * z

    for sphere in pot.spheres:
        radial(big_g * sphere.mass_solar / (r + sphere.scale_pc) ** 2)
    if pot.black_hole_solar:
        radial(big_g * pot.black_hole_solar * r / (r2 + pot.softening_pc ** 2) ** 1.5)
    if pot.halo.mass_solar:
        scale = pot.nfw_scale_pc
        u = r / scale
        radial(big_g * pot.nfw_characteristic_solar * _nfw_mass_function(u) / np.where(r > 0.0, r2, 1.0))
    cylinder2 = x * x + y * y
    for disk in pot.disks:
        root = np.sqrt(z * z + disk.scale_height_pc ** 2)
        denominator = (cylinder2 + (disk.scale_length_pc + root) ** 2) ** 1.5
        pull = big_g * disk.mass_solar / denominator
        ax, ay = ax - pull * x, ay - pull * y
        az = az - pull * z * (disk.scale_length_pc + root) / root
    if ax.ndim == 0:
        return float(ax), float(ay), float(az)
    return ax, ay, az


def circular_speed_kms(pot, radius_pc):
    """The speed, km/s, of a circular orbit of the given cylindrical radius in the plane `z = 0`: `sqrt(R |a_R|)`, so `v_c^2 = R x pull` holds by construction."""
    radius = np.asarray(radius_pc, dtype=float)
    pull_x, _, _ = acceleration(pot, radius, 0.0, 0.0)
    speed = np.sqrt(np.maximum(-radius * np.asarray(pull_x), 0.0))
    return float(speed) if speed.ndim == 0 else speed


def vertical_force(pot, radius_pc, height_pc):
    """`K_z`, the vertical pull on a star at `(R, z)` in (km/s)^2 per pc (positive towards the plane), for the vertical force of the disk near the Sun."""
    _, _, pull_z = acceleration(pot, radius_pc, 0.0, height_pc)
    return -pull_z


def from_shape(shape, galaxy_seed=None):
    """
    The `Potential` of the galaxy of `shape` (see the module doc): the density model's six components (`density.density_components`; a spheroid-bearing
    shape the density model has no component for is refused there) and the dark halo of `enclosed_mass.dark_halo` (`galaxy_seed` `None`: the median
    concentration). Masses are in solar masses, so the potential and `enclosed_mass`'s total are the same galaxy's.
    """
    spheres, disks = [], []
    for component in density.density_components(shape):
        if component.name in HERNQUIST_COMPONENTS:
            scale = component.parameters[HERNQUIST_COMPONENTS[component.name]]
            spheres.append(Hernquist(component.name, component.mass_solar, scale))
        elif component.name in DISK_COMPONENTS:
            parameters = component.parameters
            disks.append(Disk(component.name, component.mass_solar, parameters["scale_length_pc"], parameters["scale_height_pc"]))
    dark = enclosed_mass.dark_halo(shape, galaxy_seed)
    return Potential(tuple(spheres), tuple(disks), Halo(dark.dark_solar, dark.concentration, dark.virial_radius_pc),
                     density.central_black_hole_mass_solar(shape))
