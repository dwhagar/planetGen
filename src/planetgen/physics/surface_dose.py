# planetgen/physics/surface_dose.py

"""
Surface dose kernels (GEN.444, docs/design/composition/radiation-model.md sections 5 and 6)
==========================================================================================

Pure functions for the parts of the ground-level dose that the proxies of `physics/radiation.py` stand in for: the
stellar-particle spectrum and the share of it that gets past a column and a magnetic cutoff, the nuclear-interaction
shielding of a mixed column, and the radiation weights that turn absorbed dose into equivalent dose. Nothing here
draws, stores or reads a body; `radiation.generate_dose` calls these once GEN.443, GEN.441 and GEN.459 give it the
inputs (the second slice of GEN.444). The cosmic-ray table and the muon floor stay in `radiation.py` (GEN.453); the
ground's gamma dose and radon come from the radiogenic module of GEN.459.

Stellar particles (R14). The integral rigidity spectrum of Band form with the parameters, and the
share of the fluence above 30 MeV that has a rigidity above both the atmospheric cutoff and the geomagnetic cutoff
(area-weighted over latitude with the dipole's `cos^4` law).

Shielding (section 5). A column of several gases is as thick for a proton as `sum X_k (90.1 / lambda_k)` g/cm2 of air,
with the nuclear interaction lengths of the Particle Data Group (`INTERACTION_LENGTH_G_CM2`). A gas not in the table is
an error, not a guess.
"""

import math

from planetgen.physics.elements import DataError

# --- Stellar particles ------------------------------------------------------------------------------------------

PROTON_REST_MEV = 938.272
SEP_REFERENCE_MEV = 30.0
"""float: The energy the event's fluence is quoted above."""


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


def rigidity_gv(kinetic_mev):
    """The rigidity, GV, of a proton of this kinetic energy."""
    return math.sqrt(kinetic_mev * kinetic_mev + 2.0 * kinetic_mev * PROTON_REST_MEV) / 1000.0


def share_past_cutoffs(spectrum, atmospheric_cutoff_gv, geomagnetic_cutoff_gv=0.0, steps=360):
    """The share of the event's fluence above 30 MeV that reaches the ground: the fluence above
    `max(atmospheric cutoff, Rc cos^4(latitude))`, averaged over the latitudes with the weight `cos(latitude)` (equal
    areas of a sphere), over the fluence above 30 MeV; the floor is never below the 30 MeV rigidity, so no cutoff gives 1."""
    if steps < 2:
        raise ValueError("share_past_cutoffs needs at least 2 latitude steps")
    lowest = rigidity_gv(SEP_REFERENCE_MEV)
    reference = spectrum.flux(lowest)
    total = 0.0
    weights = 0.0
    for i in range(steps + 1):
        latitude = (math.pi / 2.0) * i / steps
        weight = math.cos(latitude) * (0.5 if i in (0, steps) else 1.0)
        floor = max(lowest, atmospheric_cutoff_gv, geomagnetic_cutoff_gv * math.cos(latitude) ** 4)
        total += weight * spectrum.flux(floor)
        weights += weight
    return min(total / weights / reference, 1.0)


# --- Shielding by a mixed column --------------------------------------------------------------------------------

AIR_INTERACTION_LENGTH_G_CM2 = 90.1
INTERACTION_LENGTH_G_CM2 = {"air": 90.1, "N2": 90.1, "O2": 90.1, "CO2": 88.9, "H2": 52.0}
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


# --- Radiation weights ------------------------------------------------------------------------------------------

WEIGHT_PHOTON_ELECTRON_MUON = 1.0
WEIGHT_PROTON_PION = 2.0
WEIGHT_ALPHA_HEAVY_ION = 20.0
"""float: ICRP 103 radiation weights wR (secondary source, R3; checked against the IAEA GSR Part 3 values)."""


def neutron_weight(energy_mev):
    """The ICRP 103 neutron radiation weight at `energy_mev`: `2.5 + 18.2 exp(-[ln E]^2 / 6)` below 1 MeV,
    `5.0 + 17.0 exp(-[ln 2E]^2 / 6)` from 1 to 50 MeV and `2.5 + 3.25 exp(-[ln 0.04E]^2 / 6)` above (peak 20.7 near 1 MeV,
    floor 2.5)."""
    if energy_mev <= 0:
        raise ValueError("neutron_weight needs a positive energy")
    if energy_mev < 1.0:
        return 2.5 + 18.2 * math.exp(-math.log(energy_mev) ** 2 / 6.0)
    if energy_mev <= 50.0:
        return 5.0 + 17.0 * math.exp(-math.log(2.0 * energy_mev) ** 2 / 6.0)
    return 2.5 + 3.25 * math.exp(-math.log(0.04 * energy_mev) ** 2 / 6.0)
