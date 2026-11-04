"""
The air a world's field keeps from the stellar wind, and the air an open world loses to it (GEN.804, the lock-free half of GEN.800).

A pure module: a field's magnetopause, a planet's size and the XUV in, the mass of ions that leave out. It reads nothing from a body and draws
nothing. Thermal escape (`thermosphere`) lifts neutrals out of the well; this is the other channel, the one a field changes: the wind sweeps up
ions the XUV made at the top of the air (pickup, sputtering, the polar wind) and carries them off. Venus and Mars, with no global field, lose
them at about 1e25 and 2e24 ions per second (Venus Express and MAVEN; second-hand); a world with a dipole is closed to the wind except at its
polar caps.

1. **Open share.** The field lines that reach out past the magnetopause at `r_mp` planet radii leave the dipole at the colatitude where
   `sin^2(theta) = 1 / r_mp`, so the two caps hold `1 - sqrt(1 - 1 / r_mp)` of the surface: 0.05 at ten radii (Earth), 0.18 at three, all of it
   as `r_mp` falls to 1 (no field, or one the wind presses to the surface). A world with no field has every ion exposed.
2. **Rate.** The ions the wind takes are a share `ION_ESCAPE_EFFICIENCY` of the energy-limited rate the XUV could lift (`thermosphere.energy_limited_mass_loss_kg_s`)
   over the open share. STAND-IN ION1: the efficiency is tuned so Venus and Mars lose what the spacecraft saw (both within a factor 2).

Every number that is not a published value is a stand-in in the gap registry (docs/design/research/gap-filling-plan.md section 6.67).
"""

import math

ION_ESCAPE_EFFICIENCY = 1.2e-5
"""float: STAND-IN ION1: the share of the energy-limited XUV mass loss that leaves as ions through the open field. Venus (about 0.2 kg/s at 1.9 times Earth's
XUV, 21.7 t/s energy-limited) gives 9e-6 and Mars (about 0.13 kg/s at 0.43 times, 6.5 t/s) 2e-5; the geometric mean."""


def open_share(magnetopause_rp):
    """The share of the surface the stellar wind reaches: 1 for no field (0 or `None`) or a magnetopause at or inside the surface, the polar caps of
    `1 - sqrt(1 - 1 / r_mp)` otherwise."""
    if not magnetopause_rp or magnetopause_rp <= 1.0:
        return 1.0
    return 1.0 - math.sqrt(1.0 - 1.0 / magnetopause_rp)


def ion_loss_kg_s(energy_limited_kg_s, magnetopause_rp=None):
    """The mass the stellar wind carries off as ions, kg/s, from the XUV's energy-limited rate (`thermosphere.energy_limited_mass_loss_kg_s` at an
    efficiency of 1) and the magnetopause in planet radii (0 or `None` for no field)."""
    return max(energy_limited_kg_s, 0.0) * ION_ESCAPE_EFFICIENCY * open_share(magnetopause_rp)
