# src/planetgen/physics/cloud_albedo.py

"""
The cloud albedo of a world with an ocean, on the surface-air fixed point (GEN.974; design in `docs/design/composition/cloud-albedo.md`).

A world that warms thins its low cloud and reflects less; a world locked to its star builds a deck of cloud under the sun that rises with the light. The Bond albedo of
a body with surface water is the ice law of GEN.972 (`climate_branches.IceAlbedoLaw`, the drawn albedo is its ice-free value) plus the two cloud terms here:

    A(T, S) = A_ice(T) + dA_amount(T) + dA_deck(S)

* **Amount:** `dA_amount = -BETA (T - T_REF)`, clipped to plus or minus `DELTA_MAX`. `BETA` makes the cloud feedback of the Earth at 288 K the assessment's central
  +0.42 W/m2 per K (IPCC AR6 2021, table 7.10): `BETA = 0.42 x 4 / S_EARTH`. A body at `T_REF` keeps exactly the albedo it drew.
* **Deck:** `dA_deck = +GAMMA (S - S_REF)` for a body locked 1:1 to its star in a flux above `S_REF`, clipped at `DELTA_DECK_MAX`: an ocean planet of an M-star spectrum
  has albedos of 0.42 at 1,000 W/m2 and 0.53 at 2,200 ([R0914] table 1). A body that is not locked has no deck.

Both terms are nil where there is no liquid ocean at the step. The law is pure numbers; `surface_air` reads it each step at the temperature of the step before, as it does the
ice. Stand-ins CA1 to CA5 of `gap-filling-plan.md`.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants

CLOUD_FEEDBACK_W_M2_PER_K = 0.42
"""float: The cloud feedback of the Earth at 288 K, W/m2 per K (IPCC AR6 2021, table 7.10, central value; very likely range -0.10 to +0.94). Row CA1."""
REFERENCE_SOLAR_CONSTANT_W_M2 = constants.SOLAR_CONSTANT_W_M2
"""float: The starlight the albedo feedback is sized against (the Earth's, 1,361 W/m2: the flux that gives a feedback in W/m2 per K its albedo per K). Row CA1."""
BETA_PER_K = CLOUD_FEEDBACK_W_M2_PER_K * 4.0 / REFERENCE_SOLAR_CONSTANT_W_M2
"""float: The albedo lost per kelvin of warming, 1.2e-3 per K: the feedback of `CLOUD_FEEDBACK_W_M2_PER_K` is `-BETA x S / 4`. Row CA1."""
REFERENCE_TEMPERATURE_K = 288.0
"""float: The mean surface temperature at which a wet body keeps the albedo it drew (the Earth's; the wet albedo range is centered on the Earth's 0.30). Row CA2."""
DELTA_MAX = 0.12
"""float: The most the amount term moves the albedo either way (the Earth's shortwave cloud effect, about 47 W/m2 over 340, less the share that cannot thin). Row CA2."""
DECK_GAMMA_PER_W_M2 = (0.53 - 0.42) / (2200.0 - 1000.0)
"""float: The albedo a locked ocean body gains per W/m2 of starlight above `DECK_REFERENCE_W_M2`, 9.2e-5 ([R0914] table 1: 0.42 at 1,000 and 0.53 at 2,200, an M-star spectrum). Row CA3."""
DECK_REFERENCE_W_M2 = 1000.0
"""float: The starlight below which a locked ocean body has no substellar deck. Row CA3."""
DECK_DELTA_MAX = 0.12
"""float: The most the deck term adds. Row CA3."""
MAX_ALBEDO = 0.95
"""float: The albedo the sum is held under, so the absorbed flux stays positive for any ice law."""

LOCK_Q_OVER_K2 = 100.0 / 0.3
"""float: The tidal dissipation over the Love number of a rocky planet in the despin time (Q 100, k2 0.3). A planetGen stand-in, row CA5: the spin history is drawn after the air, so the
air cannot read the state it will end in."""
LOCK_BIRTH_SPIN_HOURS = 10.0
"""float: The birth spin the despin time is worked from (a rocky planet's draw is a few hours to a day). Row CA5."""
MOMENT_FACTOR = 0.4
"""float: The moment of inertia of a uniform sphere over `m R^2`. Row CA5."""


def amount_shift(temperature_k):
    """The albedo the cloud amount adds at a mean surface temperature: negative above `REFERENCE_TEMPERATURE_K` (the low cloud thins), positive below, within `DELTA_MAX`."""
    return min(max(-BETA_PER_K * (temperature_k - REFERENCE_TEMPERATURE_K), -DELTA_MAX), DELTA_MAX)


def deck_shift(starlight_w_m2):
    """The albedo the substellar deck of a body locked to its star adds in `starlight_w_m2`: nil up to `DECK_REFERENCE_W_M2`, rising with the light to `DECK_DELTA_MAX`."""
    return min(max(DECK_GAMMA_PER_W_M2 * (starlight_w_m2 - DECK_REFERENCE_W_M2), 0.0), DECK_DELTA_MAX)


@dataclass(frozen=True)
class CloudLaw:
    """The cloud terms of a body: the starlight it is under and whether it is locked 1:1 to its star (the deck needs both)."""

    starlight_w_m2: float
    locked: bool = False

    def shift(self, temperature_k, ocean=True):
        """The albedo the clouds add at `temperature_k`; nil with no liquid ocean."""
        if not ocean:
            return 0.0
        return amount_shift(temperature_k) + (deck_shift(self.starlight_w_m2) if self.locked else 0.0)


def despin_time_s(planet_mass_kg, planet_radius_m, star_mass_kg, semi_major_axis_m):
    """The time the tide of the star takes to stop a planet born spinning at `LOCK_BIRTH_SPIN_HOURS` (Gladman and others 1996, the standard simplified form): `omega a^6 I Q /
    (3 G M^2 k2 R^5)` with `I = MOMENT_FACTOR m R^2`. Seconds."""
    omega = 2.0 * math.pi / (LOCK_BIRTH_SPIN_HOURS * 3600.0)
    return (omega * semi_major_axis_m ** 6 * MOMENT_FACTOR * planet_mass_kg * LOCK_Q_OVER_K2
            / (3.0 * constants.G * star_mass_kg ** 2 * planet_radius_m ** 3))


def tidally_locked(planet_mass_kg, planet_radius_m, star_mass_kg, semi_major_axis_m, age_s):
    """Whether the tide of the star has stopped the spin of a planet by `age_s` (the despin time is under the age): the estimate the deck is gated by."""
    if not (planet_mass_kg > 0.0 and planet_radius_m > 0.0 and star_mass_kg > 0.0 and semi_major_axis_m > 0.0):
        return False
    return despin_time_s(planet_mass_kg, planet_radius_m, star_mass_kg, semi_major_axis_m) < age_s
