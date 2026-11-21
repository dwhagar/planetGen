# planetgen/physics/hot_cloud_deck.py

"""
Where the clouds of a hot atmosphere sit, and the most each trace gas can hold (GEN.994 and GEN.1577)
=====================================================================================================

`hot_condensates` gives the temperature at which a cloud species starts to condense at a pressure (`cloud_temperature_k`) and the most molecules of a trace gas the host
leaves room for (`gas_ceiling`). This module places them on a giant: the *deck*, where each cloud species meets the giant's temperature-pressure profile, and the *ceiling* of each
trace gas, packed into the one text column `cloud_decks` of its `planets` row (`pack`, `unpack`; JSON, as `sky_suns` is).

**The profile** (`giant_profile`) is the grey radiative-equilibrium profile of `hot_giant_profile` on the giant's equilibrium temperature, its own heat flux and its gravity, floored from
below by a dry adiabat down from the giant's 1 bar temperature (`ADIABAT_EXPONENT`): the grey profile alone is nearly isothermal below 1 bar and never meets a cloud curve inside the grid
(GEN.1577 measurement), where a giant's interior is convective and hotter with depth. The 1 bar temperature is the one every giant's air is solved at (`giant_draw.ONE_BAR_OVER_EFFECTIVE`
times its effective temperature, as `mixture_subject` reads it). Stand-in 6.317 DA1.

**The deck.** A profile is a function `T(P)` in K of the pressure in bar that rises with depth. A species condenses where the profile is cooler than its condensation
curve, so the cloud base is the deepest pressure where the two meet and the cloud lies above it. `cloud_base` finds that meeting on a log-spaced grid of pressures
(`GRID_BAR`) and refines it by bisection. It returns None for a species the profile never cools enough to condense, and flags `reaches_bottom` when the profile is cooler than
the curve down to the deepest grid pressure (the cloud has no base inside the grid: its base is at least `GRID_BAR[-1]`).

**What is stored.** `generate` sets `cloud_decks` on a giant planet of a star (not a moon, not a rogue planet): the base pressure (bar) and temperature (K) of each cloud species that
forms, deepest first, and the ceiling of each trace gas in molecules per hydrogen atom. NULL for a body that is not a giant and for a giant with no star. No draw is made.
"""

import json
import math
from dataclasses import dataclass

from planetgen.physics import bulk_draw, climate_draws, constants, giant_draw, hot_condensates, hot_giant_profile

FIELDS = ("cloud_decks",)
"""tuple: The stored column of a giant planet (`generate`; `pack` and `unpack` make and read its text)."""

ADIABAT_EXPONENT = (constants.ADIABATIC_INDEX_H2_HE - 1.0) / constants.ADIABATIC_INDEX_H2_HE
"""float: `d ln T / d ln P` of the dry adiabat below a giant's 1 bar level, (gamma - 1) / gamma of the H2-He gas (2/7 at `constants.ADIABATIC_INDEX_H2_HE` 1.4, the value `mixture_giant`
takes down to its quench level). Stand-in 6.317 DA1."""

GRID_BAR = tuple(10.0 ** (exponent / 4.0) for exponent in range(-24, 13))
"""tuple: The pressures, bar, the profile is searched on: 1e-6 to 1e3, four to the decade."""

TOLERANCE = 1.0e-3
"""float: The bisection stops when the bracket is this narrow in log10 of the pressure. [Stand-in 6.317 DA3, tier T9, tweak L0.]"""


@dataclass(frozen=True)
class CloudBase:
    """One species' cloud: the deepest pressure where the profile meets the condensation curve, and the temperature there (`reaches_bottom`: no meeting inside the grid)."""

    species: str
    base_bar: float
    base_k: float
    reaches_bottom: bool = False


def _gap(species, host, profile, pressure_bar):
    """Profile minus condensation temperature, K: positive where the gas is too hot to condense."""
    return profile(pressure_bar) - hot_condensates.cloud_temperature_k(species, pressure_bar, host)


def cloud_base(species, host, profile):
    """The `CloudBase` of `species` in `host`'s gas along `profile` (a function of pressure in bar giving K), or None when the profile is nowhere cool enough."""
    gaps = [_gap(species, host, profile, pressure) for pressure in GRID_BAR]
    if gaps[-1] < 0.0:
        return CloudBase(species, GRID_BAR[-1], profile(GRID_BAR[-1]), reaches_bottom=True)
    for index in range(len(GRID_BAR) - 1, 0, -1):
        if gaps[index] >= 0.0 > gaps[index - 1]:
            low, high = math.log10(GRID_BAR[index - 1]), math.log10(GRID_BAR[index])
            while high - low > TOLERANCE:
                mid = 0.5 * (low + high)
                if _gap(species, host, profile, 10.0 ** mid) < 0.0:
                    low = mid
                else:
                    high = mid
            base = 10.0 ** (0.5 * (low + high))
            return CloudBase(species, base, profile(base))
    return None


def deck(host, profile, species=None):
    """The clouds that form along `profile`, deepest base first, as a tuple of `CloudBase`; `species` defaults to every species of `hot_condensates.CLOUD_SPECIES`."""
    found = (cloud_base(name, host, profile) for name in (species or hot_condensates.CLOUD_SPECIES))
    return tuple(sorted((base for base in found if base is not None), key=lambda base: -base.base_bar))


def trace_gas_ceilings(host):
    """The most molecules per hydrogen atom each trace gas of `hot_condensates.HOT_TRACE_GASES` can hold in `host`'s gas, by formula."""
    return {formula: hot_condensates.gas_ceiling(formula, host) for formula in hot_condensates.HOT_TRACE_GASES}


def one_bar_temperature_k(planet, starlight_w_m2):
    """The temperature at a giant's 1 bar level, K: `giant_draw.ONE_BAR_OVER_EFFECTIVE` times the temperature it radiates at (the starlight `starlight_w_m2` and its own heat), as its air is solved."""
    return giant_draw.ONE_BAR_OVER_EFFECTIVE * giant_draw.effective_temperature_k(starlight_w_m2, giant_draw.own_flux_w_m2(planet))


def giant_profile(planet, starlight_w_m2):
    """The profile (a function of pressure in bar giving K) of a giant lit by `starlight_w_m2`: the grey profile floored by the dry adiabat from its 1 bar temperature, or None without a Bond albedo."""
    albedo = getattr(planet, "bond_albedo", None)
    if albedo is None:
        return None
    grey = hot_giant_profile.profile_function(climate_draws.equilibrium_temperature_k(albedo, starlight_w_m2), giant_draw.own_flux_w_m2(planet),
                                              planet.gravity * constants.EARTH_GRAVITY)
    one_bar_k = one_bar_temperature_k(planet, starlight_w_m2)
    return lambda pressure_bar: max(grey(pressure_bar), one_bar_k * pressure_bar ** ADIABAT_EXPONENT)


def _rounded(value):
    return None if value is None else float(f"{value:.4g}")


def pack(clouds, ceilings):
    """The text of the column: JSON `{"clouds": [[species, base bar, base K, 1 if no base inside the grid], ...], "ceilings": {formula: molecules per H atom}}`, deepest cloud first."""
    return json.dumps({"clouds": [[base.species, _rounded(base.base_bar), _rounded(base.base_k), int(base.reaches_bottom)] for base in clouds],
                       "ceilings": {formula: _rounded(value) for formula, value in ceilings.items()}}, separators=(",", ":"))


@dataclass(frozen=True)
class Record:
    """What a giant's column holds: its `clouds` (`CloudBase`, deepest first) and the `ceilings` of its trace gases, formula to molecules per hydrogen atom."""

    clouds: tuple
    ceilings: dict


def unpack(text):
    """The `Record` of a column's text, or None for NULL or an empty text."""
    if not text:
        return None
    data = json.loads(text)
    return Record(tuple(CloudBase(species, bar, kelvin, bool(bottom)) for species, bar, kelvin, bottom in data["clouds"]), dict(data["ceilings"]))


def generate(planet, starlight_w_m2):
    """
    Sets `cloud_decks` on `planet` (call after its air, which sets its Bond albedo and gravity; no draw): the deck along `giant_profile` and the trace-gas ceilings of its star's gas,
    packed; None for a moon, a body that is not a giant, a giant with no star or no [Fe/H], and one with no starlight.
    """
    planet.cloud_decks = None
    star = getattr(planet, "star", None)
    if planet.is_moon or planet.body_type != "g" or star is None or getattr(star, "starless", False) or not starlight_w_m2 or starlight_w_m2 <= 0.0:
        return
    host = bulk_draw.composition_of(star)
    profile = giant_profile(planet, starlight_w_m2)
    if host is None or profile is None:
        return
    planet.cloud_decks = pack(deck(host, profile), trace_gas_ceilings(host))
