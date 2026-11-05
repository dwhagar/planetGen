# planetgen/physics/nitrogen_ice.py

"""
Nitrogen ice from the disk model (GEN.1218a, GEN.1225)
=======================================================

The comet table has no N2, so no body held nitrogen ice and no Triton-like nitrogen sea could form. The research of GEN.1225 puts nitrogen ice in two steps.

1. **Where it forms.** N2 condenses at the temperature where its saturation pressure equals its partial pressure in the disk gas: the nitrogen over the hydrogen
   (`DISK_N2_PER_H2`, all of the Sun's nitrogen as N2: log eps 7.83 against 12 for hydrogen [R0169]) times the disk's total pressure (`DISK_PRESSURE_PA`, 1e-4 bar, the
   pressure of the condensation table `condensate_shares`). The vapor pressure is beta-N2 [R1445] (`saturation`). That is 28.7 K; a body whose ice formed in a zone warmer
   than that (the blackbody disk temperature at its formation distance, as the water line is found in `bulk_draw`) gets no N2 ice, and its nitrogen is the ammonia
   of its ice, as before.
2. **How much.** The solar ceiling is 0.107 of the water's mass as N2 (all nitrogen as N2, all oxygen in water; `SOLAR_N2_PER_WATER`). What cold bodies keep is 1e-4 to 1.4e-3
   of the water mass (67P bulk 1.4e-3 [R3835: Table 2]; the cometary model of Pluto 1.2e-4 to 7.8e-4 [R1599]), a trapping efficiency of 0.1 to 1.3 percent of the ceiling
   (`TRAPPING_EFFICIENCY_RANGE`), drawn log-uniform per body.

Limits. Retention is not modelled (a small body at 20 to 50 K loses N2; the sources give the regime and no radius cutoff), so every body that forms cold keeps its N2
(stand-in GN3). The two N2 to CO ratios of 67P (the coma value against the bulk) disagree five-fold; the range spans both.
"""

import math

from planetgen.physics import saturation

SOLAR_LOG_EPSILON_N = 7.83
"""float: log10 of the Sun's N atoms per 1e12 H atoms, with the hydrogen at 12 ([R0169], Asplund and others 2009)."""

DISK_N2_PER_H2 = 10.0 ** (SOLAR_LOG_EPSILON_N - 12.0)
"""float: N2 molecules per H2 molecule of the disk gas with all the nitrogen as N2: N2 / H2 = (N / 2) / (H / 2) = N / H, 6.76e-5."""

DISK_PRESSURE_PA = 10.0
"""float: The disk's total pressure, 1e-4 bar (the pressure of the equilibrium condensation of `condensate_shares`). Stand-in GN2."""

SOLAR_N2_PER_WATER = 0.107
"""float: The mass of N2 over the mass of water if all the Sun's nitrogen is N2 and the oxygen the rock leaves is water (GEN.1225 section 3.4). Source tier T3, tweak level L0."""

TRAPPING_EFFICIENCY_RANGE = (0.001, 0.013)
"""tuple: The share of the solar ceiling that cold ice keeps, from 67P and the Pluto cometary model: 1e-4 to 1.4e-3 of the water mass. Stand-in GN1."""


def condensation_temperature_k():
    """The temperature (K) below which the disk's N2 is ice: the beta-N2 saturation curve at the nitrogen's partial pressure in the disk gas."""
    return saturation.condensation_temperature_k("N2", DISK_N2_PER_H2 * DISK_PRESSURE_PA)


def zone_is_cold(temperature_k):
    """Whether a disk zone at `temperature_k` is colder than N2 condenses at."""
    return temperature_k < condensation_temperature_k()


def percent_of_water(trapping):
    """N2 as percent of the water's mass for a trapping efficiency `trapping` (a share of the ceiling), the unit of the comet volatiles."""
    if not 0.0 <= trapping <= 1.0:
        raise ValueError(f"a trapping efficiency is a share, not {trapping!r}")
    return 100.0 * SOLAR_N2_PER_WATER * trapping
