# planetgen/physics/giant_colour.py

"""
The atmosphere colour of gas and ice giants (GEN.419; docs/design/composition/atmosphere-colour.md section 11).

A giant has no ground, so its colour is the colour of the light seen from a reference level: the air above that level scatters and
absorbs (Rayleigh from hydrogen and helium, the gas table's bands for methane and the rest), and below it lies a cloud deck or the deep
gas as the lower boundary. It is `atmosphere_colour.atmosphere_colour` with the column cut at the reference pressure and the cloud
deck's reflectance as the ground.

* **The reference level** (decision 23) is `P_ref = min(1 bar, cloud top)`. The cloud top follows from the cloud's depth: the condensate
  decays upward as `exp(-f_sed z / L)` while the air thins as `exp(-z / H)`, so the depth measured from the top reaches 1 at
  `P_top = P_base tau_c^(-1 / (1 + f_sed H / L))` (not built: no cloud deck is stored on a giant yet, so `reference_pressure_pa` is given none and the level stays at 1 bar). A cloud-free giant keeps 1 bar.
* **The cloud deck** is an `Aerosol` (the depth, single-scattering albedo and asymmetry of GEN.414); its diffuse reflectance over a black
  base, per wavelength, is the Lambertian ground of the colour function (`deck_albedo`). A clear giant has the dark `DEEP_ALBEDO`
  below the level (a stand-in: the deep gas absorbs; the alkali wings, TiO and VO of the hot giants are not in the gas table).
* **Band strengths** come from the gas table's coarse bins. The methane, ammonia and water band strengths of Karkoschka and of the
  Oxford k-tables (GEN.422) are not built, so a giant with any of them above `PARTIAL_ABOVE` in mole fraction carries
  `colour_quality = "partial"`.

The air above the level is the stored mixture rescaled to `P_ref` (well mixed; the depletion of a condensable above its own cloud is the
caller's, through the fractions it passes).
"""

from dataclasses import dataclass

import numpy as np

from planetgen.physics import air_column, atmosphere_colour as colour_model, constants, photometry, two_stream

REFERENCE_BAR_PA = 1.0e5
"""float: The radius level of a giant (decision 9), Pa; the reference level is never deeper."""

DEEP_ALBEDO = 0.03
"""float: The reflectance of the deep clear gas below the reference level (stand-in)."""

PARTIAL_ABOVE = 1e-4
"""float: A band-bearing gas above this mole fraction makes the colour `partial` until the band strengths are read."""

BAND_GASES = ("CH4", "NH3", "H2O")
"""tuple: The gases whose visible bands are only the gas table's coarse bins."""

SUDARSKY_EDGES_K = ((150.0, "I"), (350.0, "II"), (900.0, "III"), (1500.0, "IV"))
"""tuple: Upper equilibrium temperatures of the Sudarsky classes (approximate); hotter is V."""


@dataclass(frozen=True)
class GiantColour:
    """
    The colour of a giant's air.

    Attributes:
        atmosphere (AtmosphereColour): The colour at the reference level.
        p_ref_pa (float): The reference pressure, Pa.
        colour_quality (str): `full`, or `partial` while a band-bearing gas is present and its strengths are the coarse bins.
        has_cloud_deck (bool): The lower boundary is a cloud deck, not the deep gas.
    """

    atmosphere: colour_model.AtmosphereColour
    p_ref_pa: float
    colour_quality: str
    has_cloud_deck: bool


def reference_pressure_pa(cloud_top_pa=None):
    """`min(1 bar, cloud top)`, Pa; 1 bar for a giant with no cloud."""
    return REFERENCE_BAR_PA if cloud_top_pa is None else float(min(REFERENCE_BAR_PA, cloud_top_pa))


def deck_albedo(deck, mu0=colour_model.DEFAULT_MU0):
    """The diffuse reflectance of a cloud deck over a black base on the colour grid: a Lambertian ground for the air above it."""
    tau, omega, g = deck.optics(colour_model.WAVELENGTHS_NM)
    return two_stream.scattering_layer(tau, omega, g, mu0).reflectance_above


def giant_atmosphere_colour(partials_pa, gravity_m_s2, *, temperature_k, p_ref_pa=REFERENCE_BAR_PA, cloud_deck=None,
                            mu0=colour_model.DEFAULT_MU0, haze=None, spread_x=0.0, spread_y=0.0,
                            light=None):
    """
    The colour of a giant's air at its reference level.

    Args:
        partials_pa (dict): The mixture, gas-table formula to partial pressure (any scale; only the fractions are used).
        gravity_m_s2 (float): The gravity at the reference level.
        temperature_k (float): The star's effective temperature, K.
        p_ref_pa (float, optional): The reference pressure (`reference_pressure_pa`).
        cloud_deck (Aerosol, optional): The cloud at the level, the lower boundary; the dark deep gas when `None`.
        mu0 (float, optional): The cosine of the star's zenith angle.
        haze (Aerosol, optional): Haze above the level.
        spread_x, spread_y (float, optional): A cool star's seeded chromaticity draws.
        light (ndarray, optional): The star's relative spectral power (`atmosphere_colour.body_light`); the blackbody when `None`.

    Returns:
        GiantColour

    Raises:
        air_column.EmptyColumn: There is no mixture or no gravity.
        ValueError: `p_ref_pa` is not above zero.
    """
    if not p_ref_pa > 0.0:
        raise ValueError("the reference pressure must be above zero")
    total = sum(value for value in partials_pa.values() if value and value > 0.0)
    if total <= 0.0:
        raise air_column.EmptyColumn("the giant has no mixture to take a colour from")
    column = air_column.from_partials({gas: value / total * p_ref_pa for gas, value in partials_pa.items()}, gravity_m_s2)
    albedo = deck_albedo(cloud_deck, mu0) if cloud_deck is not None and cloud_deck.tau_550 > 0.0 else DEEP_ALBEDO
    result = colour_model.atmosphere_colour(column, temperature_k=temperature_k, ground_albedo=albedo, mu0=mu0, haze=haze,
                                            spread_x=spread_x, spread_y=spread_y, light=light)
    bands = any(column.fraction(gas) > PARTIAL_ABOVE for gas in BAND_GASES)
    return GiantColour(result, float(p_ref_pa), "partial" if bands else "full", isinstance(albedo, np.ndarray))


def body_giant_colour(body):
    """
    The `GiantColour` of a gas or ice giant row from what it stores (its mixture, gravity and star's temperature), read at 1 bar over the
    dark deep gas: the stored cloud decks do not exist yet (DB.68), so this is the colour of the air above the deck, which is also what a
    limb of the planet shows. `None` for a body with no mixture, no gravity or no star temperature.
    """
    temperature = getattr(getattr(body, "star", None), "temperature", None)
    gravity = (getattr(body, "gravity", None) or 0.0) * constants.EARTH_GRAVITY
    if not temperature or not photometry.gives_visible_light(temperature) or gravity <= 0.0:
        return None
    mixture = air_column.partials_pa(body)
    if sum(mixture.values()) <= 0.0:
        return None
    return giant_atmosphere_colour(mixture, gravity, temperature_k=temperature, light=colour_model.body_light(body.star))
