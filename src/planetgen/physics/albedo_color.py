# planetgen/physics/albedo_color.py

"""
How the colour of the star changes the Bond albedo of a bright surface (GEN.912, registry 6.122)
================================================================================================

The Bond albedo draw of `climate_draws.albedo_range` was the same for every star. Ice, snow and the ocean under its clouds reflect the blue
and visible light of a hot star well and the near-infrared light of a cool star badly (snow falls from about 0.8 in the visible to 0.3 to 0.5
beyond 1 micron), and Rayleigh scattering is weaker in the red; so a planet of a K or M dwarf reflects less of the light it gets than the same
planet of a G star and is warmer at the same flux. [R2202] find the broadband albedo of snow 0.8 in the Sun's
light and 0.6 in the light of a 3300 K black body, and of ice 0.5 and 0.3. [R4065] give three measures for the
Sun-like G star, the F2V HD128167, the K2V HD22049 and the M3V AD Leo:

* Table 2 (a 1-D column with an Earth-like air, clouds and Rayleigh scattering, read from a text extraction): the broadband planetary albedo of
  a snow-covered planet 0.668 (F), 0.646 (G), 0.604 (K), 0.398 (M) and of an ocean planet 0.329, 0.319, 0.302, 0.234;
* Sec. 4.2.1 (the GCM): the mean albedo of the mostly ice-covered planet at the edge of global ice is 0.58 (F, at 0.99 of the modern solar
  constant), 0.55 (G, 0.93) and 0.51 (M, 0.74);
* the freezing line itself: a planet of the M star stays free of global ice with 19 percent less flux than one of the G star, and one of the F
  star needs 6 percent more (0.99 against 0.93), which is the same absorbed power at an albedo shifted by a factor 0.79 (M) and 1.05 (F) from 0.55.

The three disagree for M by a factor of 0.62, 0.79 and 0.93; the module takes the middle for the **frozen** range (0.79 M, 0.935 K from the table, 1.05
F) and, for the **wet** range (ocean and cloud), the ocean row of Table 2 at `WET_SHARE` (half) of its shift, since clouds are grey and snow
is the surface that loses most. The factors are ratios to the G star, as a function of the star's effective temperature (the four stars'
temperatures are not in the paper; the module takes the main-sequence values of the spectral types). The dry range is left alone: the table's
land entry is rock under a thick air of clouds and Rayleigh scattering, and a dry body has none of that to lose (the colour effect of bare
regolith is not in the papers). The freezing line of `climate_draws.albedo_range` (`tuning.ALBEDO_FROZEN_BELOW_K`) is left where it was.

Use:

    factor = spectral_factor("frozen", star_temperature_k)         # 1.0 for the Sun
    low, high = shifted_range(climate_draws.ALBEDO_FROZEN_RANGE, "frozen", star_temperature_k)
"""

SPECTRAL_TEMPERATURES_K = (3410.0, 5040.0, 5770.0, 6900.0)
"""tuple: The effective temperatures of the M3V, K2V, G2V and F2V stars of [R4065], K (main-sequence values of the types; the paper gives none, row SAB1)."""

SPECTRAL_FACTOR = {
    "frozen": (0.79, 0.935, 1.0, 1.05),
    "wet": (0.234 / 0.319, 0.302 / 0.319, 1.0, 0.329 / 0.319),
}
"""dict: The Bond albedo of a bright world of each star of `SPECTRAL_TEMPERATURES_K` over the Sun-like one. `frozen`: M 0.79 and F 1.05 from the freezing
lines of [R4065] (19 percent less flux for M, 6 percent more for F; the GCM mean albedos 0.51 and 0.58 over 0.55 and Table 2's snow row 0.398 over
0.646 bracket the M one), K 0.935 from Table 2's snow row (0.604 over 0.646). `wet`: the ocean row of Table 2 over the G entry (row SAB1)."""

WET_SHARE = 0.5
"""float: The share of the ocean row's shift that the wet range (ocean under cloud) takes: clouds are spectrally grey and snow is the surface that loses most
(row SAB2)."""

def spectral_factor(regime, star_temperature_k):
    """
    The ratio of the Bond albedo of a `regime` ("frozen" or "wet") planet of a star of `star_temperature_k` to the Sun-like one: linear in temperature between
    the four stars of `SPECTRAL_TEMPERATURES_K`, held at the end values beyond them. 1.0 for any other regime (dry rock) and for a star with no temperature.
    """
    factors = SPECTRAL_FACTOR.get(regime)
    if factors is None or not star_temperature_k or star_temperature_k <= 0.0:
        return 1.0
    temps = SPECTRAL_TEMPERATURES_K
    t = min(max(star_temperature_k, temps[0]), temps[-1])
    factor = factors[-1]
    for i in range(len(temps) - 1):
        if t <= temps[i + 1]:
            share = (t - temps[i]) / (temps[i + 1] - temps[i])
            factor = factors[i] + share * (factors[i + 1] - factors[i])
            break
    if regime == "wet":
        factor = 1.0 + WET_SHARE * (factor - 1.0)
    return factor


def shifted_range(albedo_range, regime, star_temperature_k):
    """`albedo_range` (low, high) times `spectral_factor`, held under 0.95 so a bright world of a hot star stays a possible albedo."""
    factor = spectral_factor(regime, star_temperature_k)
    return min(albedo_range[0] * factor, 0.95), min(albedo_range[1] * factor, 0.95)
