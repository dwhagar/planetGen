# planetgen/physics/dust_model.py

"""
The dust of a galaxy from its own gas (VIEW.43, goal #4907; reports/view20-realistic-dust-2026-11-05.md).

The first sky used the dust of the Milky Way written in fixed parsecs (`sky_photometry.dust_av_per_pc` without a profile): 8,000 pc for the
Sun, a scale length of 2,260 pc, a layer 134 pc thick. In a galaxy of another size it was wrong by the ratio of the sizes, 9.6 mag per kpc in a
galaxy of 481 pc scale length and 2.3 mag from the Sun to the Galactic Center where the center is seen through 30 to 45. Dust is gas that carries
metals, and the generator has the gas of every ring of the disk (`galaxy.disk_history`, GEN.469), the metallicity of that gas and the
dust-to-gas law of its metals (`physics.nebula_composition`, GEN.1269), and the gravity that sets how thick a gas layer is. So the dust of a place is:

1. **The column.** The gas surface density `Sigma_gas(R)` of the ring, its hydrogen share `X`, and `N_H / A_V = 2.21e21` hydrogen atoms per cm2 per
   magnitude ([R5130: Eq. 1], 22 supernova remnants, X-ray N_H against optical A_V) make the V extinction through the whole layer, both sides:

       A_V(column) = f_vis * Sigma_gas * X / (m_H * 2.21e21) * [dust-to-gas(Z(R)) / dust-to-gas(solar)]

   `Z(R)` is the gas metallicity at the ring's radius (`nebula_composition.GAS_GRADIENT_DEX_PER_KPC`, zero at the galaxy's reference radius) and the
   dust-to-gas ratio is the broken power law of [R5102] (`nebula_composition.gas_to_dust_law`), linear in the metals above 12 + log(O/H) 7.96.
   `f_vis` is the share of the gas column that dims the stars a catalogue or a person can use (stand-in DU2: the model column at the Sun is twice the
   measured high-latitude one, 0.24 mag through the layer, [R5132]; stars behind dense clouds are not in the maps and part of the gas is not in the thin layer).

2. **The layer.** Gas in the disk's own gravity is a sech^2 sheet of scale height `h(R) = sigma_g^2 / (pi G (Sigma_gas + Sigma_stars))`, with the
   velocity dispersion of the neutral gas `sigma_g` (stand-in DU3). It is 113 pc at the Sun and flares to 470 pc at 12.5 kpc, as the dust does
   (thin and thick disks of [R5131] at 12 to 13 kpc: 340 and 1,000 pc), because the surface density falls outward. A floor `LAYER_FLOOR_PC` (DU4) stands in for
   the bulge's pull and the central molecular zone's thickness, which the surface density of the disk alone leaves out (VIEW.46 replaces it with the vertical force of the
   potential, GEN.115).

3. **The extinction per parsec** is `A_V(column) / (2 h) sech^2(z / h)`: the sech^2 layer holds `2 h` of height, so the column is its integral. Two table lookups
   and a hyperbolic cosine a point, the cost of the fixed law it replaces.

Against what is measured (the Milky Way, scale length 2,600 pc, 5e10 solar masses; none but the first was a tuning target): the column through the layer at the Sun
0.24 mag; the plane at the Sun 1.06 mag per kpc (0.85 +- 0.20, [R5132]; 0.8 to 1.0, [R5131]); the plane at 9.5, 10.5 and 12.5 kpc 0.51, 0.24 and 0.08 (0.75, 0.60 and 0.14 from
the thin and thick disks of [R5131]); the Sun to the Galactic Center 39 mag with the 40 pc floor (A_V 31 to 59 in the literature, 38 to 42 for the model of [R5134]). The center moves
by a factor of 6 with the floor, so it is consistent, not confirmed. A galaxy with less gas, fewer metals and a thicker layer is clearer: the small test galaxy gives 0.14 mag per kpc.

What is left out: the spiral arms' contrast (the gas rings are azimuthally smooth), the warp, a thin and a thick disk of dust (a third of the dust is thick, [R5131]), the
central molecular zone ring (GEN.1292), clumpy clouds (the mean is used; dense clouds and nebulae are VIEW.8's) and reddening (VIEW.11). Nothing here is stored: a profile is worked out from the
galaxy's shape, which is stored, and is kept in memory.
"""

import functools
import hashlib
import math

import numpy as np

from planetgen import tuning
from planetgen.physics import constants, nebula_composition, stellar_abundances, unit_registry

HYDROGEN_ATOMS_PER_AV_CM2 = 2.21e21
"""float: Hydrogen atoms per cm2 for one magnitude of V extinction in the Galaxy: N_H = (2.21 +- 0.09)e21 A_V, from 22 supernova remnants [R5130: Eq. 1]. Stand-in DU1,
tier T1, tweak L0."""

VISIBLE_SHARE = 0.5
"""float: The share of the gas column that dims the stars a catalogue or a person can use: the model column at the Sun (0.47 mag for 11.8 solar masses per pc2) over the
measured one through a half layer at high latitude doubled (0.24, [R5132]); the maps leave out the stars behind dense clouds. Stand-in DU2, tier T6, tweak L2."""

GAS_DISPERSION_KMS = 9.0
"""float: The velocity dispersion of the neutral gas, km/s, that holds a layer of scale height sigma^2 / (pi G Sigma) up (the HI of the Galaxy is 7 to 11, recalled, not read). Stand-in
DU3, tier T2, tweak L1."""

LAYER_FLOOR_PC = 40.0
"""float: The thinnest layer, parsecs: the thickness of the central molecular zone (a few tens of pc), which stands in for the bulge's pull the disk's own gravity leaves out. Stand-in
DU4, tier T6, tweak L2; VIEW.46 replaces it with the vertical force of the potential."""

SOLAR_12LOG = nebula_composition.SOLAR_OXYGEN_12LOG
"""float: The solar 12 + log(O/H) the dust-to-gas law is written against."""

_PC_CM = constants.PARSEC_M * 100.0
_MSUN_KG = unit_registry.converter("Msun", "kg")(1.0)
"""float: One solar mass in kg, through the unit registry (GEN.568); the same number as `constants.SOLAR_MASS_TO_KG`."""
_MSUN_PC2_IN_HYDROGEN_PER_CM2 = _MSUN_KG / (_PC_CM / 100.0) ** 2 / constants.HYDROGEN_ATOM_MASS_KG / 1.0e4
"""float: Hydrogen atoms per cm2 in one solar mass of hydrogen per pc2 (1.25e20)."""

_G_PC_KMS2_PER_MSUN = constants.G * _MSUN_KG / constants.PARSEC_M / 1.0e6
"""float: The gravitational constant in pc (km/s)^2 per solar mass (4.30e-3)."""


@functools.lru_cache(maxsize=1)
def hydrogen_mass_fraction():
    """X, the hydrogen share of the gas by mass: that of the solar thin-disk birth composition the nebulae and the stars share."""
    offsets = stellar_abundances.abundance_offsets(0.0, stellar_abundances.BRANCH_THIN, (0.0,) * len(stellar_abundances.RESIDUAL_ELEMENTS))
    return stellar_abundances.composition(0.0, offsets).x


def gas_oh_dex(radius_pc, reference_radius_pc):
    """[O/H] of the gas at `radius_pc` (array or number), dex: zero at `reference_radius_pc`, falling by `nebula_composition.GAS_GRADIENT_DEX_PER_KPC` a kpc."""
    return -nebula_composition.GAS_GRADIENT_DEX_PER_KPC * (np.asarray(radius_pc, dtype=float) - reference_radius_pc) / 1000.0


def dust_to_gas_ratio(oh_dex):
    """The dust-to-gas ratio at gas [O/H] `oh_dex` (array or number) over the solar one: the law of [R5102] (`nebula_composition.gas_to_dust_law`)."""
    solar = nebula_composition.gas_to_dust_law(SOLAR_12LOG)
    return np.array([solar / nebula_composition.gas_to_dust_law(SOLAR_12LOG + float(v)) for v in np.atleast_1d(oh_dex)]).reshape(np.shape(oh_dex))


def layer_height_pc(surface_density_msun_pc2, dispersion_kms=GAS_DISPERSION_KMS, floor_pc=LAYER_FLOOR_PC):
    """
    The scale height, pc, of a gas sheet of velocity dispersion `dispersion_kms` in the gravity of `surface_density_msun_pc2` (gas and stars together):
    `sigma^2 / (pi G Sigma)`, no thinner than `floor_pc`. Arrays or numbers.
    """
    sigma = np.asarray(surface_density_msun_pc2, dtype=float)
    return np.maximum(dispersion_kms ** 2 / (math.pi * _G_PC_KMS2_PER_MSUN * np.maximum(sigma, 1.0e-12)), floor_pc)


def column_av(gas_msun_pc2, oh_dex, visible_share=VISIBLE_SHARE):
    """The V extinction through the whole layer, both sides, of gas of `gas_msun_pc2` and gas [O/H] `oh_dex` (the column of the module text)."""
    hydrogen = np.asarray(gas_msun_pc2, dtype=float) * hydrogen_mass_fraction() * _MSUN_PC2_IN_HYDROGEN_PER_CM2
    return visible_share * hydrogen / HYDROGEN_ATOMS_PER_AV_CM2 * dust_to_gas_ratio(oh_dex)


class DustProfile:
    """
    The dust of one galaxy as two tables against galactocentric radius: the V extinction through the layer (`column_av`) and the layer's scale height (pc). Built by `profile_from_history`
    and `profile_for_shape`; hashable by its `key`.

    Attributes:
        radius_pc (tuple): The radii of the table, pc, increasing.
        column_av (tuple): The extinction through the whole layer, both sides, at each radius.
        height_pc (tuple): The scale height of the sech^2 layer at each radius.
    """

    __slots__ = ("radius_pc", "column_av", "height_pc", "_r", "_c", "_h")

    def __init__(self, radius_pc, column_av, height_pc):
        self.radius_pc = tuple(float(v) for v in radius_pc)
        self.column_av = tuple(float(v) for v in column_av)
        self.height_pc = tuple(float(v) for v in height_pc)
        if not (len(self.radius_pc) == len(self.column_av) == len(self.height_pc) >= 2):
            raise ValueError("the tables must be of one length, two or more")
        if any(b <= a for a, b in zip(self.radius_pc, self.radius_pc[1:])):
            raise ValueError("the radii must increase")
        if any(v < 0 for v in self.column_av) or any(v <= 0 for v in self.height_pc):
            raise ValueError("a column must be zero or more and a height positive")
        self._r, self._c, self._h = (np.array(self.radius_pc), np.array(self.column_av), np.array(self.height_pc))

    def __eq__(self, other):
        return isinstance(other, DustProfile) and (self.radius_pc, self.column_av, self.height_pc) == (other.radius_pc, other.column_av, other.height_pc)

    def __hash__(self):
        return hash((self.radius_pc, self.column_av, self.height_pc))

    def column_at(self, radius_pc):
        """The extinction through the layer at `radius_pc`; the first ring's inside it, none beyond the last."""
        return np.interp(np.asarray(radius_pc, dtype=float), self._r, self._c, right=0.0)

    def height_at(self, radius_pc):
        """The layer's scale height at `radius_pc`, held at the ends of the table."""
        return np.interp(np.asarray(radius_pc, dtype=float), self._r, self._h)

    def av_per_pc(self, x_pc, y_pc, z_pc):
        """The extinction in magnitudes per parsec at galactocentric `(x, y, z)` in parsecs (arrays or numbers, z up)."""
        radius = np.hypot(np.asarray(x_pc, dtype=float), np.asarray(y_pc, dtype=float))
        height = self.height_at(radius)
        return self.column_at(radius) / (2.0 * height) / np.cosh(np.asarray(z_pc, dtype=float) / height) ** 2

    def plane_av_per_kpc(self, radius_pc):
        """The extinction in the plane at `radius_pc`, magnitudes per kpc."""
        return 1000.0 * self.column_at(radius_pc) / (2.0 * self.height_at(radius_pc))

    def key(self):
        """A short text that names these tables (for a cache key): a digest of the values rounded to six figures."""
        text = repr((tuple(round(v, 6) for v in self.radius_pc), tuple(float("%.6g" % v) for v in self.column_av), tuple(float("%.6g" % v) for v in self.height_pc)))
        return hashlib.sha1(text.encode(), usedforsecurity=False).hexdigest()[:16]


def profile_from_history(history, disk_scale_length_pc, reference_radius_pc=None, *, visible_share=VISIBLE_SHARE, dispersion_kms=GAS_DISPERSION_KMS,
                         floor_pc=LAYER_FLOOR_PC):
    """
    The `DustProfile` of the disk whose rings are `history` (a `galaxy.disk_history.DiskHistory`) and whose stellar disk has scale length `disk_scale_length_pc`.

    The gas and the metallicity are the ring's own; the layer is held up by the gravity of the ring's gas and stars. `reference_radius_pc` is where the gas is solar
    (`tuning.REFERENCE_RADIUS_SCALE_LENGTHS` scale lengths by default).
    """
    if reference_radius_pc is None:
        reference_radius_pc = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * disk_scale_length_pc
    radius = np.asarray(history.radius_scale_lengths, dtype=float) * disk_scale_length_pc
    gas = np.asarray(history.gas_now, dtype=float)
    stars = np.asarray(history.stars_now, dtype=float)
    oh = gas_oh_dex(radius, reference_radius_pc)
    return DustProfile(radius, column_av(gas, oh, visible_share), layer_height_pc(gas + stars, dispersion_kms, floor_pc))


def profile_for_shape(shape):
    """The `DustProfile` of the galaxy of this shape (`galaxy.disk_history.history_of`): the gas it ends its history with, at its own scale length."""
    return _profile_for_shape(shape._replace(k_norm=1.0))


@functools.lru_cache(maxsize=64)
def _profile_for_shape(shape):
    from planetgen.galaxy import disk_history
    return profile_from_history(disk_history.history_of(shape), shape.disk_scale_length_pc)


def dust_for_shape(shape):
    """
    The `dust=` dict of the sky functions for the galaxy of this shape (VIEW.44): its own `DustProfile` under `"profile"`. The sky, the glow band and the tier table take it
    as they take any override (`sky_photometry.dust_av_per_pc`), so the default dust of a galaxy is the dust of its own gas.
    """
    return {"profile": profile_for_shape(shape)}
