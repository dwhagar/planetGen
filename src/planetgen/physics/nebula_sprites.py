# planetgen/physics/nebula_sprites.py

"""
The sprites of nebulae and supernova remnants in the sky of a system (VIEW.28, S10 of reports/system-sky-from-the-system-2026-11-04.md)
=========================================================================================================================================

A nebula or a remnant is not a point: the sky draws it as an ellipse sprite with an angular size, a surface brightness and the extinction it puts in
front of what lies behind it. All of it is worked out on arrays from the geometry the light-travel code gives (`light_travel.retarded_offsets`: the
centre where it was when its light left it) and the stored contents of the object (`nebulae`, `supernova_remnants`):

* **Angular radius** `asin(R / d)` for an observer outside the ball of radius `R` at distance `d` (the silhouette of a sphere; the plan's `atan(R / d)`
  is the same to 1 percent beyond three radii); a ball that holds the observer covers the sky (`pi`).
* **Extinction** `extinction_av` is the dimming through the middle of the ball (a diameter); a line of sight at impact parameter `b` crosses a chord
  `2 sqrt(R^2 - b^2)`, so it dims by `extinction_av * sqrt(1 - (b / R)^2)` (`chord_extinction_mag`, the first tier of the dust of the sky plan section 2.5).
  An observer inside the ball sees half the diameter (`extinction_av / 2`, the line from the centre to the edge).
* **Surface brightness** (V magnitudes per square arcsecond, the mean over the disc) by what makes the light:
  - *emission* (H II regions and the emission-reflection nebulae, the planetary nebulae and the Wolf-Rayet rings, family `emission` and `planetary`):
    the hydrogen-alpha light of an ionised gas. Its emission measure over the mean chord of a uniform ball (`4R/3`) is `EM = n^2 * 4R/3` (cm^-6 pc), the
    line is `EM / 2.75` Rayleighs at 10,000 K, and a Rayleigh of H-alpha is worked out into the V band through the filter (`rayleigh_to_v_mag`);
  - *reflection* (family `reflection`): the light of its central star scattered by the grains, `albedo * (1 - exp(-tau)) / (4 pi)` of the star's
    flux at half the radius;
  - *supernova remnants*: the swept-up shell of a strong shock, four times the ambient density over a twelfth of the radius (the Sedov-Taylor
    shell), seen as hydrogen-alpha at a small share of a photoionised gas of the same emission measure; a plerion has no model (`None`);
  - *dark* and *diffuse* clouds give no light of their own (`None`): they dim, they do not shine.
* **Visibility**: an object is a naked-eye sprite only when its mean surface brightness is under about 23 magnitudes per square arcsecond
  (docs/design/orbits/sky-view.md 2.10).
* **Remnants** that exploded less than their light time ago have not been seen (`light_travel.remnant_seen`): the sprite is `seen = False`.

Every constant that is a choice is a registry row (GEN.725, section 6.173, rows NS1 to NS7).

Use:

    sprites = nebula_sprites.sprite_data(family, radius_ly, distance_ly, density_cm3, extinction_av, foreground_av=0.0, age_years=None, light_time_s=None)
    # {"angular_radius_rad", "inside", "extinction_mag", "surface_brightness", "total_magnitude", "visible", "seen"}  (arrays)
"""

import math

import numpy as np

from planetgen.physics import constants, light_travel

PC_PER_LY = constants.LY_TO_M / constants.PARSEC_M
"""float: Parsecs in a light year."""

ARCSEC2_PER_SR = (180.0 * 3600.0 / math.pi) ** 2
"""float: Square arcseconds in a steradian."""

EM_PER_RAYLEIGH = 2.75
"""float: STAND-IN NS1 (Confident Unverified, L0): emission measure (cm^-6 pc) per Rayleigh of hydrogen-alpha for gas at 10,000 K (the standard relation, recalled)."""

RAYLEIGH_PHOTONS_CM2_S_SR = 1.0e6 / (4.0 * math.pi)
"""float: Photons per square centimetre per second per steradian in one Rayleigh (its definition)."""

HALPHA_WAVELENGTH_A = 6562.8
"""float: The hydrogen-alpha wavelength (angstroms)."""

V_ZERO_FLUX_ERG_CM2_S_A = 3.63e-9
"""float: STAND-IN NS2 (Confident Unverified, L0): the flux density of a star of V = 0 (erg per square centimetre per second per angstrom), the Johnson V zero point (recalled)."""

V_EFFECTIVE_WIDTH_A = 880.0
"""float: STAND-IN NS2: the effective width of the V band (angstroms, recalled)."""

HALPHA_V_TRANSMISSION = 0.4
"""float: STAND-IN NS3 (Marginal Unverified, L3): the share of a hydrogen-alpha line that the V filter passes (the filter's red edge, recalled)."""

SHELL_COMPRESSION = 4.0
"""float: STAND-IN NS4 (Confident educated guess from physics, L0): the density of a strong adiabatic shock over the ambient density, (gamma + 1) / (gamma - 1) for gamma = 5/3."""

SHELL_THICKNESS_OVER_RADIUS = 1.0 / 12.0
"""float: STAND-IN NS4: the thickness of a Sedov-Taylor shell over the radius (recalled)."""

SHELL_HALPHA_SHARE = 0.5
"""float: STAND-IN NS5 (Marginal Unverified, L3): the share of the hydrogen-alpha of a photoionised gas of the same emission measure that a shocked shell gives, tuned so that an old shell remnant (n = 1 cm^-3, 50 ly) is about 1 to 10 Rayleighs, the order of the dim optical remnants (recalled, not looked up)."""

REFLECTION_ALBEDO = 0.5
"""float: STAND-IN NS6 (Marginal Unverified, L1): the albedo of interstellar grains in the V band (recalled)."""

REFLECTION_CENTRAL_ABSOLUTE_V = -3.0
"""float: STAND-IN NS7 (Marginal educated guess from physics, L3): the absolute V magnitude of the B star that lights a reflection nebula (class F, about 1,000 solar luminosities), set so a 10 ly reflection nebula comes out near 22 to 23, the order of the bright ones (recalled)."""

NAKED_EYE_LIMIT_MAG_ARCSEC2 = 23.0
"""float: A sprite is visible when its mean surface brightness is under this (sky-view.md 2.10; a dark site's sky is 21.8)."""

EMISSION_FAMILIES = ("emission", "planetary")
"""tuple: The nebula families whose light is the hydrogen-alpha of an ionised gas."""

REMNANT_FAMILIES = ("supernova-remnant", "plerion")
"""tuple: The families of a supernova remnant: its shell (`"supernova-remnant"`) and a plerion, whose pulsar wind nebula has no model of its light."""


def angular_radius_rad(radius_ly, distance_ly):
    """The angular radius of balls of `radius_ly` seen from `distance_ly` from their centres: `asin(R / d)`; `pi` where the observer is inside (d <= R)."""
    radius = np.asarray(radius_ly, dtype=float)
    distance = np.asarray(distance_ly, dtype=float)
    inside = distance <= radius
    ratio = np.where(inside, 1.0, radius / np.where(distance > 0.0, distance, 1.0))
    return np.where(inside, math.pi, np.arcsin(np.minimum(ratio, 1.0)))


def chord_extinction_mag(extinction_av, radius_ly, impact_ly):
    """The dimming of a line of sight at `impact_ly` from the centre of a ball: `extinction_av * sqrt(1 - (b / R)^2)`; 0 outside the ball's silhouette."""
    radius = np.asarray(radius_ly, dtype=float)
    ratio = np.asarray(impact_ly, dtype=float) / np.where(radius > 0.0, radius, 1.0)
    inside_silhouette = ratio < 1.0
    return np.where(inside_silhouette, np.asarray(extinction_av, dtype=float) * np.sqrt(np.where(inside_silhouette, 1.0 - ratio * ratio, 0.0)), 0.0)


def sprite_extinction_mag(extinction_av, radius_ly, distance_ly):
    """The dimming a sprite puts in front of what lies behind its centre: the diameter's `extinction_av` for an observer outside, half of it for one inside."""
    extinction = np.asarray(extinction_av, dtype=float)
    return np.where(np.asarray(distance_ly, dtype=float) <= np.asarray(radius_ly, dtype=float), 0.5 * extinction, extinction)


def emission_measure_cm6_pc(density_cm3, radius_ly):
    """`n^2 * 4R/3` in cm^-6 pc: the emission measure over the mean chord of a uniform ball."""
    return np.asarray(density_cm3, dtype=float) ** 2 * (4.0 / 3.0) * np.asarray(radius_ly, dtype=float) * PC_PER_LY


def rayleigh_to_v_mag(rayleigh):
    """The V magnitude per square arcsecond of a hydrogen-alpha line of `rayleigh` Rayleighs, through the V filter (`HALPHA_V_TRANSMISSION`); 1 R is about 30.4."""
    photon_erg = 2.0 * math.pi * constants.REDUCED_PLANCK * constants.SPEED_OF_LIGHT_M_S / (HALPHA_WAVELENGTH_A * 1.0e-10) * 1.0e7
    line = np.asarray(rayleigh, dtype=float) * RAYLEIGH_PHOTONS_CM2_S_SR * photon_erg / ARCSEC2_PER_SR      # erg / cm2 / s / arcsec2
    reference = V_ZERO_FLUX_ERG_CM2_S_A * V_EFFECTIVE_WIDTH_A
    with np.errstate(divide="ignore"):
        return -2.5 * np.log10(line * HALPHA_V_TRANSMISSION / reference)


def emission_surface_brightness(density_cm3, radius_ly):
    """The mean V surface brightness (mag per square arcsecond) of an ionised gas of `density_cm3` in a ball of `radius_ly`, from its hydrogen-alpha."""
    return rayleigh_to_v_mag(emission_measure_cm6_pc(density_cm3, radius_ly) / EM_PER_RAYLEIGH)


def remnant_surface_brightness(density_cm3, radius_ly):
    """The mean V surface brightness of a supernova remnant's shocked shell: `SHELL_COMPRESSION` times the ambient density over `SHELL_THICKNESS_OVER_RADIUS` of the radius, at `SHELL_HALPHA_SHARE`."""
    density = np.asarray(density_cm3, dtype=float) * SHELL_COMPRESSION
    measure = density ** 2 * SHELL_THICKNESS_OVER_RADIUS * np.asarray(radius_ly, dtype=float) * PC_PER_LY
    return rayleigh_to_v_mag(SHELL_HALPHA_SHARE * measure / EM_PER_RAYLEIGH)


def reflection_surface_brightness(radius_ly, extinction_av, central_absolute_v=REFLECTION_CENTRAL_ABSOLUTE_V):
    """The mean V surface brightness of a reflection nebula: the central star's V flux at half the radius, scattered by `REFLECTION_ALBEDO` of the grains' optical depth, per steradian."""
    radius_pc = np.asarray(radius_ly, dtype=float) * PC_PER_LY
    star_mag = central_absolute_v + 5.0 * np.log10(np.maximum(radius_pc / 2.0, 1e-6) / 10.0)
    tau = np.asarray(extinction_av, dtype=float) / 1.0857
    scattered = REFLECTION_ALBEDO * (1.0 - np.exp(-tau)) / (4.0 * math.pi)
    with np.errstate(divide="ignore"):
        return star_mag - 2.5 * np.log10(scattered) + 2.5 * math.log10(ARCSEC2_PER_SR)


def surface_brightness(family, density_cm3, radius_ly, extinction_av):
    """
    The mean V surface brightness (mag per square arcsecond) of each object by its family, or `nan` where it has no light of its own (dark and diffuse
    clouds, and the plerions, which have no model). `family` is an array of the family names of `tuning.NEBULA_FAMILIES`, with `"supernova-remnant"`
    for a remnant and `"plerion"` for one whose morphology is a plerion.
    """
    family = np.asarray(family)
    out = np.full(family.shape, np.nan)
    emission = np.isin(family, EMISSION_FAMILIES)
    out = np.where(emission, emission_surface_brightness(density_cm3, radius_ly), out)
    out = np.where(family == "reflection", reflection_surface_brightness(radius_ly, extinction_av), out)
    out = np.where(family == "supernova-remnant", remnant_surface_brightness(density_cm3, radius_ly), out)
    return out


def total_magnitude(surface_brightness_mag, angular_radius_rad_):
    """The integrated V magnitude of a disc of mean surface brightness `surface_brightness_mag` and angular radius `angular_radius_rad_`: `S - 2.5 log10(pi theta^2)` with theta in arcseconds."""
    theta_arcsec = np.asarray(angular_radius_rad_, dtype=float) * (180.0 * 3600.0 / math.pi)
    with np.errstate(divide="ignore"):
        return np.asarray(surface_brightness_mag, dtype=float) - 2.5 * np.log10(math.pi * theta_arcsec * theta_arcsec)


def sprite_data(family, radius_ly, distance_ly, density_cm3, extinction_av, foreground_av=0.0, age_years=None, light_time_s=None):
    """
    The sprite of each object at `distance_ly` from the observer (its retarded centre, `light_travel.retarded_offsets`).

    Args:
        family (array): The family of each object (see `surface_brightness`).
        radius_ly, distance_ly, density_cm3, extinction_av (array): The stored radius, the distance to the retarded centre, the stored density and the stored extinction.
        foreground_av (float or array): The dust of the sky between the observer and the object (magnitudes), which dims its surface brightness.
        age_years (array, optional): The age of each remnant (`supernova_remnants.age_years`; `nan` for a nebula). With `light_time_s` it decides `seen`.
        light_time_s (array, optional): The light time of each object in seconds.

    Returns:
        dict: Arrays `angular_radius_rad`, `inside` (the observer is within the ball), `extinction_mag` (what the sprite puts in front of what lies behind it),
            `surface_brightness` (observed, with the foreground; `nan` for no light of its own), `total_magnitude`, `visible` (under
            `NAKED_EYE_LIMIT_MAG_ARCSEC2`) and `seen` (a remnant's explosion has been seen; True for every nebula).
    """
    radius = np.asarray(radius_ly, dtype=float)
    distance = np.asarray(distance_ly, dtype=float)
    theta = angular_radius_rad(radius, distance)
    observed = surface_brightness(family, density_cm3, radius, extinction_av) + np.asarray(foreground_av, dtype=float)
    seen = np.ones(radius.shape, dtype=bool)
    if age_years is not None and light_time_s is not None:
        age_s = np.asarray(age_years, dtype=float) * constants.SECONDS_PER_YEAR
        remnant = np.isin(np.asarray(family), REMNANT_FAMILIES)
        seen = np.where(remnant, light_travel.remnant_seen(age_s, light_time_s), True)
    return {
        "angular_radius_rad": theta,
        "inside": distance <= radius,
        "extinction_mag": sprite_extinction_mag(extinction_av, radius, distance),
        "surface_brightness": observed,
        "total_magnitude": total_magnitude(observed, theta),
        "visible": np.isfinite(observed) & (observed < NAKED_EYE_LIMIT_MAG_ARCSEC2) & seen,
        "seen": seen,
    }
