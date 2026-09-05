"""
Calibration tests for the unified evolved-star wind model.

`Star._calculate_heliosphere_radius_static`'s giant/subgiant/bright-giant/
supergiant/hypergiant tiers were three independently hand-fit pieces (a
separate hypergiant power-law, plus Reimers' Law with two different eta
values) stitched together at Yerkes-class boundaries -- exactly the
structure that let the hypergiant constant drift ~9 orders of magnitude out
of calibration unnoticed earlier this session. They're now one formula:
Nieuwenhuijzen & de Jager (1990), Mdot = 9.63e-15 * L^1.42 * M^0.16 * R^0.81
(L, M, R in solar units), with a literature-cited 10x correction above
log(L/Lsun) > 5 (Mauron & Josselin 2011).

These tests use fixed, hand-computed inputs (not the randomized star
generator) so the underlying mass-loss rate can be inverted from the
returned heliosphere radius and checked directly against real measured
mass-loss rates, independent of any other part of the generation pipeline.

Run with: pytest tests/test_heliosphere_model.py
"""
import math

from stellarObjects.starData import Star
from stellarObjects import physical_constants as pc


def _implied_mass_loss_rate_msun_per_year(heliosphere_au, mass_kg, radius_km, wind_velocity_factor):
    """
    Inverts Star._calculate_heliosphere_radius_static's momentum-balance
    equation (R = sqrt(Mdot * v_inf / (4 * pi * P_ism))) to recover the
    Mdot (in Msun/yr) that produced a given heliosphere radius, given the
    star's mass/radius (for escape velocity) and which wind-velocity tier
    applied (GIANT_WIND_VELOCITY_FACTOR or HYPERGIANT_WIND_VELOCITY_FACTOR,
    both multiples of escape velocity).
    """
    radius_m = radius_km * pc.KM_TO_M_FACTOR
    escape_velocity = math.sqrt(pc.ESCAPE_VELOCITY_CONSTANT * pc.G * mass_kg / radius_m)
    wind_velocity = wind_velocity_factor * escape_velocity

    heliopause_radius_m = heliosphere_au * pc.AU_TO_M
    momentum_flux = (heliopause_radius_m ** 2) * pc.FOUR_PI * pc.ISM_PRESSURE
    mass_loss_rate_kgs = momentum_flux / wind_velocity
    return mass_loss_rate_kgs * pc.SECONDS_PER_YEAR / pc.SOLAR_MASS_TO_KG


def _solar_radius_from_teff(lum_sol, teff_k, sun_teff_k=5778):
    """R/Rsun from the Stefan-Boltzmann law, given L/Lsun and Teff (K)."""
    return math.sqrt(lum_sol) / (teff_k / sun_teff_k) ** 2


def test_betelgeuse_like_supergiant_matches_measured_mass_loss_rate():
    """
    Betelgeuse (M-type red supergiant): L~1.26e5 Lsun, M~17 Msun, Teff~3600K.
    Measured mass-loss rate estimates cluster around 2e-6 to 1e-5 Msun/yr.
    This luminosity is just above the log(L/Lsun) > 5 correction threshold,
    so this also exercises that correction.
    """
    lum_sol, mass_sol, teff_k = 1.26e5, 17.0, 3600
    radius_sol = _solar_radius_from_teff(lum_sol, teff_k)
    mass_kg = mass_sol * pc.SOLAR_MASS_TO_KG
    luminosity_w = lum_sol * pc.SOLAR_LUMINOSITY
    radius_km = radius_sol * pc.SOLAR_RADIUS_M / pc.KM_TO_M_FACTOR

    helio_au = Star._calculate_heliosphere_radius_static(mass_kg, luminosity_w, radius_km, "M2IA", "IA")
    implied_mdot = _implied_mass_loss_rate_msun_per_year(helio_au, mass_kg, radius_km, pc.GIANT_WIND_VELOCITY_FACTOR)

    assert 5e-7 <= implied_mdot <= 5e-5, f"implied Mdot {implied_mdot:.3g} Msun/yr is far from Betelgeuse's measured rate"


def test_hypergiant_matches_lbv_mass_loss_rate_range():
    """
    A hot (Teff~20,000K) hypergiant at L=1e6 Lsun, M=100 Msun: real
    LBVs/hypergiants (e.g. Eta Carinae's quiescent wind) are typically
    ~1e-5 to 1e-4 Msun/yr.
    """
    lum_sol, mass_sol, teff_k = 1e6, 100.0, 20000
    radius_sol = _solar_radius_from_teff(lum_sol, teff_k)
    mass_kg = mass_sol * pc.SOLAR_MASS_TO_KG
    luminosity_w = lum_sol * pc.SOLAR_LUMINOSITY
    radius_km = radius_sol * pc.SOLAR_RADIUS_M / pc.KM_TO_M_FACTOR

    helio_au = Star._calculate_heliosphere_radius_static(mass_kg, luminosity_w, radius_km, "B0", "0")
    implied_mdot = _implied_mass_loss_rate_msun_per_year(helio_au, mass_kg, radius_km, pc.HYPERGIANT_WIND_VELOCITY_FACTOR)

    assert 1e-6 <= implied_mdot <= 5e-4, f"implied Mdot {implied_mdot:.3g} Msun/yr is outside the LBV/hypergiant range"


def test_ordinary_red_giant_mass_loss_is_much_weaker_than_supergiant():
    """
    A modest red giant (L=100 Lsun, well below the high-luminosity
    correction threshold) should have a mass-loss rate many orders of
    magnitude below a supergiant's -- real RGB rates are ~1e-11 to 1e-8
    Msun/yr, versus supergiants' ~1e-7 to 1e-4 Msun/yr.
    """
    lum_sol, mass_sol, teff_k = 100.0, 1.5, 4000
    radius_sol = _solar_radius_from_teff(lum_sol, teff_k)
    mass_kg = mass_sol * pc.SOLAR_MASS_TO_KG
    luminosity_w = lum_sol * pc.SOLAR_LUMINOSITY
    radius_km = radius_sol * pc.SOLAR_RADIUS_M / pc.KM_TO_M_FACTOR

    helio_au = Star._calculate_heliosphere_radius_static(mass_kg, luminosity_w, radius_km, "K0III", "III")
    implied_mdot = _implied_mass_loss_rate_msun_per_year(helio_au, mass_kg, radius_km, pc.GIANT_WIND_VELOCITY_FACTOR)

    assert 1e-13 <= implied_mdot <= 1e-8, f"implied Mdot {implied_mdot:.3g} Msun/yr is not giant-like"


def test_high_luminosity_correction_only_applies_above_threshold():
    """
    The Mauron & Josselin (2011) 10x correction should produce a visible
    downward kink in Mdot (and therefore heliosphere) right at the
    log(L/Lsun) > 5 threshold, holding mass/radius/Teff fixed.
    """
    mass_sol, teff_k = 15.0, 4500
    mass_kg = mass_sol * pc.SOLAR_MASS_TO_KG

    def helio_at(lum_sol):
        radius_sol = _solar_radius_from_teff(lum_sol, teff_k)
        radius_km = radius_sol * pc.SOLAR_RADIUS_M / pc.KM_TO_M_FACTOR
        luminosity_w = lum_sol * pc.SOLAR_LUMINOSITY
        return Star._calculate_heliosphere_radius_static(mass_kg, luminosity_w, radius_km, "K0IB", "IB")

    just_below = helio_at(pc.NDJ_HIGH_LUMINOSITY_THRESHOLD_LSUN * 0.999)
    just_above = helio_at(pc.NDJ_HIGH_LUMINOSITY_THRESHOLD_LSUN * 1.001)

    # Without the correction, heliosphere would keep rising smoothly across
    # this boundary (Mdot ~ L^1.42, radius ~ sqrt(Mdot) ~ L^0.71, so a 0.2%
    # increase in L should only raise helio by ~0.14%). The 10x Mdot cut
    # applied just above the threshold instead makes it visibly *drop*.
    assert just_above < just_below * 0.5


def test_evolved_tier_never_produces_galactic_scale_heliosphere():
    """
    Regression guard for the original bug: even at the absolute extreme of
    the hypergiant mass/luminosity/radius ranges, heliosphere must stay far
    below galactic scale (~1e11 AU was the pre-fix value; the Milky Way is
    ~100,000 ly across, i.e. ~6.3e9 AU).
    """
    mass_kg = 150.0 * pc.SOLAR_MASS_TO_KG
    lum_sol = 2_000_000.0
    luminosity_w = lum_sol * pc.SOLAR_LUMINOSITY
    # A large, cool hypergiant radius (extreme end -- red hypergiants like
    # VY Canis Majoris are ~1000-2000 Rsun).
    radius_km = 2000.0 * pc.SOLAR_RADIUS_M / pc.KM_TO_M_FACTOR

    helio_au = Star._calculate_heliosphere_radius_static(mass_kg, luminosity_w, radius_km, "M5", "0")
    assert helio_au < 1e9, f"heliosphere {helio_au:.3g} AU approaches galactic scale"
