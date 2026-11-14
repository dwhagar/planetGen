# planetgen/physics/surface_dose.py

"""
Surface dose kernels (GEN.444, docs/design/composition/radiation-model.md sections 5 and 6)
==========================================================================================

Pure pieces of the ground-level dose that `physics/radiation.py` reads: the Band form of the stellar-particle spectrum
(`solar_particles` builds its measured events from it) and the nuclear-interaction shielding of a mixed column. The share of
an event that gets past a column and a magnetic cutoff is `solar_particles.dose_share` with `proton_range`
(NIST PSTAR, mixtures by Bragg's rule); the first slice's own copy of it, its rigidity function and its ICRP neutron weight
were retired (GEN.444): the share duplicated the production code, and the cosmic-ray table is already a dose weighted for
its neutrons, so a second weight would count them twice.

Shielding (section 5). A column of several gases is as thick for a proton as `sum X_k (90.1 / lambda_k)` g/cm2 of air,
with the nuclear interaction lengths of the Particle Data Group (`INTERACTION_LENGTH_G_CM2`). A gas not in the table is
an error, not a guess.
"""

import math

from planetgen.physics.elements import DataError

# --- Stellar particles ------------------------------------------------------------------------------------------


class BandSpectrum:
    """
    The integral rigidity spectrum of a solar particle event, `J(>R) = J0 R^-g1 exp(-R/R0)` up to `R1 = (g2 - g1) R0`
    and `J0 R1^-g1 exp(-R1/R0) (R/R1)^-g2` above, continuous at `R1` (the published form, not the
    discontinuous one of the research note). Rigidity in GV; `J0` sets the unit of the flux and cancels in a share.
    """

    def __init__(self, j0, gamma1, gamma2, r0_gv):
        if j0 <= 0 or r0_gv <= 0 or gamma2 <= gamma1:
            raise DataError("a Band spectrum needs J0 > 0, R0 > 0 and gamma2 > gamma1")
        self.j0, self.gamma1, self.gamma2, self.r0_gv = j0, gamma1, gamma2, r0_gv
        self.r1_gv = (gamma2 - gamma1) * r0_gv

    def flux(self, rigidity_gv):
        """`J(>R)`."""
        r = max(rigidity_gv, 1e-9)
        if r <= self.r1_gv:
            return self.j0 * r ** -self.gamma1 * math.exp(-r / self.r0_gv)
        return (self.j0 * self.r1_gv ** -self.gamma1 * math.exp(-self.r1_gv / self.r0_gv)
                * (r / self.r1_gv) ** -self.gamma2)


# --- Shielding by a mixed column --------------------------------------------------------------------------------

AIR_INTERACTION_LENGTH_G_CM2 = 90.1
INTERACTION_LENGTH_G_CM2 = {"air": AIR_INTERACTION_LENGTH_G_CM2, "N2": AIR_INTERACTION_LENGTH_G_CM2, "O2": AIR_INTERACTION_LENGTH_G_CM2, "CO2": 88.9, "H2": 52.0}
"""dict: The nuclear interaction length, g/cm2, of the gases the research lanes read from the Particle Data Group
(air 90.1, CO2 88.9, H2 52.0; N2 and O2 take air's). A gas not here is an error."""


def effective_air_column_g_cm2(columns_g_cm2):
    """The column, g/cm2 of air, that stops as many primaries as `columns_g_cm2`, a mapping of a gas formula (or
    "air") to its own column in g/cm2: `sum X_k (90.1 / lambda_k)`."""
    total = 0.0
    for gas, column in columns_g_cm2.items():
        try:
            length = INTERACTION_LENGTH_G_CM2[gas]
        except KeyError:
            raise DataError("no nuclear interaction length for %r; the gases are %s"
                            % (gas, ", ".join(INTERACTION_LENGTH_G_CM2))) from None
        total += column * AIR_INTERACTION_LENGTH_G_CM2 / length
    return total
