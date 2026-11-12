# planetgen/physics/surface_dose.py

"""
The Band spectrum of a stellar particle event (GEN.444, docs/design/composition/radiation-model.md section 6)
============================================================================================================

The integral rigidity spectrum of Band form (R14) that `solar_particles` sums over the 47 measured events. The share of it past a column and a magnetic cutoff is
`solar_particles.dose_share`, the column's stopping of a proton is `proton_range` (a mixed air through `proton_range.gas_material`), and the neutron and other
radiation weights are inside the CORSIKA total biological dose of `radiation` (GEN.453); the first slice of GEN.444 kept its own copies of those three, and they were
retired (GEN.1559) when the production code took their place. Nothing here draws, stores or reads a body.
"""

import math

from planetgen.physics.elements import DataError


class BandSpectrum:
    """
    The integral rigidity spectrum of a solar particle event, `J(>R) = J0 R^-g1 exp(-R/R0)` up to `R1 = (g2 - g1) R0`
    and `J0 R1^-g1 exp(-R1/R0) (R/R1)^-g2` above, continuous at `R1` (the published form, not the
    discontinuous one of the research note). Rigidity in GV; `J0` sets the unit of the flux and cancels in `share`.
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


EVENT_2006_12_13 = BandSpectrum(1.33e8, 1.045, 5.80, 0.1765)
EVENT_1976 = BandSpectrum(7.62e6, 1.710, 6.33, 0.2054)
"""The two events of the SEP table that research lanes read (10B 1.2a): the 2006-12-13 event and the 1976 event."""
