# planetgen/physics/body_radiation.py

"""
The light at a planet's or moon's ground, from what the body stores (GEN.443; docs/design/composition/radiation-model.md sections 3 and 7).

`surface_spectrum(body)` gathers the stored values the kernel of `atmospheric_radiation` needs and calls it: the air column of
`air_column.column_of`, the star's temperature and its flux at the body's distance, the star's UV-C share (`activity`), the ozone column
(`ozone_column_du`, GEN.413) and the haze, dust and cloud (`aerosols.of_body`, GEN.414). Nothing is drawn and nothing is stored: the spectrum is a
pure function of the body's record, so the dose of GEN.444, the habitability of GEN.445 and the sky colour of GEN.418 read one place.

The ground's reflectance is a stand-in until GEN.462 gives it by band: `DEFAULT_GROUND_ALBEDO`, the Earth-like land mean of 0.15 (labelled).
A star with no stored UV-C share (a remnant, or a row saved before schema v137) gives the blackbody (1).

`uv_index` and `dna_dose_w_m2` are the two numbers the habitability code wants; both are for the star at the zenith angle given (`mu0`, the default
is overhead) with the ground mixture, and 0 for a body with no air or no star.

`dna_dose_index(body)` is the number a body stores as `uv_surface_index` (GEN.445): the Setlow DNA-damage dose at its ground as a multiple of the
dose at Earth's (`earth_dna_dose_w_m2`, so Earth is 1 by definition and the limits of `protection_limits` stay in "x Earth"). A body with air reads the
spectrum below its air; a body with none reads the star's unfiltered light (`bare_dna_dose_w_m2`), the bare-ground case the proxy gave its ceiling
for; a body with no star or no distance has no index (`None`).
"""

import functools

import numpy as np

from planetgen.physics import activity, aerosols, air_column, atmospheric_radiation, constants, stellar_spectrum, uv_cross_sections

DEFAULT_GROUND_ALBEDO = 0.15
"""float: The ground reflectance used until GEN.462 gives it by band (stand-in: the Earth-like land mean)."""
EARTH_OZONE_DU = 300.0
"""float: The ozone column of the reference Earth, DU (measured: 279 to 313 DU in the models and 300 DU in the record)."""
EARTH_BOND_ALBEDO = 0.306
"""float: The Earth's Bond albedo (measured), for the temperature of its ozone layer in the reference dose."""
SUN_EFFECTIVE_K = 5772.0
"""float: The Sun's effective temperature, K (IAU 2015)."""
EARTH_AIR_PA = {"N2": 79100.0, "O2": 21200.0, "Ar": 950.0}
"""dict: The reference Earth's dry air, Pa by gas formula (101.25 kPa)."""


def _star_distance_au(body):
    """The body's distance from its star: the one `magnetism.update_exposure` set, else a planet's own orbit (a moon's is its planet's, so a
    moon without it has none)."""
    distance = getattr(body, "_star_distance_au", None)
    if distance is None and not getattr(body, "is_moon", False):
        distance = getattr(body, "distance", None)
    return distance


def _star_flux_w_m2(body, distance_au):
    star = getattr(body, "star", None)
    luminosity = getattr(star, "luminosity", None)
    temperature = getattr(star, "temperature", None)
    if not luminosity or not temperature or not distance_au:
        return None
    return temperature, stellar_spectrum.bolometric_flux_w_m2(luminosity, distance_au)


def surface_spectrum(body, mu0=1.0, ground_albedo=DEFAULT_GROUND_ALBEDO, distance_au=None):
    """
    The `SurfaceSpectrum` at the ground of `body` for the star at cosine `mu0` of the zenith angle (0.2 to 1), or `None` for a body with no air,
    no star or no gravity.
    `distance_au` is the body's distance from its star when the body does not carry it.
    """
    column = air_column.column_of(body)
    lit = _star_flux_w_m2(body, distance_au or _star_distance_au(body))
    if column is None or lit is None:
        return None
    temperature, flux = lit
    star = body.star
    share = getattr(star, "uv_c_share", None)
    stored = aerosols.of_body(body)
    centres = uv_cross_sections.centres(uv_cross_sections.standard_edges())
    haze, cloud = aerosols.layers(stored, centres)
    return atmospheric_radiation.surface_radiation(
        column, temperature_k=temperature, flux_w_m2=flux, uv_c_share=1.0 if share is None else share, uv_shares=activity.uv_shares_of(star),
        ozone_du=getattr(body, "ozone_column_du", None) or 0.0,
        ozone_temperature_k=atmospheric_radiation.ozone_layer_temperature_k(flux, getattr(body, "bond_albedo", None)), mu0=mu0,
        ground_albedo=ground_albedo,
        haze=haze, low_aerosol=cloud, cloud_cover=stored.cloud_cover)


def uv_index(body, mu0=1.0, ground_albedo=DEFAULT_GROUND_ALBEDO, distance_au=None):
    """The CIE UV index at the ground of `body`, 0 where there is no air to compute it in."""
    spectrum = surface_spectrum(body, mu0, ground_albedo, distance_au)
    return 0.0 if spectrum is None else spectrum.uv_index()


def dna_dose_w_m2(body, mu0=1.0, ground_albedo=DEFAULT_GROUND_ALBEDO, distance_au=None):
    """The Setlow DNA-damage weighted irradiance at the ground of `body`, W/m2, 0 where there is no air to compute it in."""
    spectrum = surface_spectrum(body, mu0, ground_albedo, distance_au)
    return 0.0 if spectrum is None else spectrum.dna_dose_w_m2()


def bare_dna_dose_w_m2(body, mu0=1.0, distance_au=None):
    """The Setlow DNA-damage weighted irradiance of the star's unfiltered light on a horizontal surface, W/m2, or `None` for a body with no star or no
    distance. The ground of a body with no air."""
    lit = _star_flux_w_m2(body, distance_au or _star_distance_au(body))
    if lit is None:
        return None
    temperature, flux = lit
    share = getattr(body.star, "uv_c_share", None)
    edges = uv_cross_sections.standard_edges()
    top = stellar_spectrum.irradiance_w_m2(edges, temperature, flux, 1.0 if share is None else share, activity.uv_shares_of(body.star))
    return float(mu0 * np.sum(top * atmospheric_radiation.setlow_weight(uv_cross_sections.centres(edges))))


@functools.lru_cache(maxsize=1)
def earth_dna_dose_w_m2():
    """The DNA dose at the ground of the reference Earth, W/m2: the Sun at 1 au (blackbody UV-C), `EARTH_AIR_PA` at 1 g, `EARTH_OZONE_DU`, a clear sky,
    the default ground and the Sun overhead. The unit of `dna_dose_index`."""
    column = air_column.from_partials(EARTH_AIR_PA, constants.EARTH_GRAVITY)
    flux = stellar_spectrum.bolometric_flux_w_m2(constants.SOLAR_LUMINOSITY, 1.0)
    spectrum = atmospheric_radiation.surface_radiation(
        column, temperature_k=SUN_EFFECTIVE_K, flux_w_m2=flux, ozone_du=EARTH_OZONE_DU,
        ozone_temperature_k=atmospheric_radiation.ozone_layer_temperature_k(flux, EARTH_BOND_ALBEDO))
    return spectrum.dna_dose_w_m2()


def dna_dose_index(body, mu0=1.0, ground_albedo=DEFAULT_GROUND_ALBEDO, distance_au=None):
    """The DNA dose at the ground of `body` as a multiple of Earth's (`earth_dna_dose_w_m2`), or `None` for a body with no star or no distance."""
    if air_column.column_of(body) is None:
        dose = bare_dna_dose_w_m2(body, mu0, distance_au)
    else:
        spectrum = surface_spectrum(body, mu0, ground_albedo, distance_au)
        dose = None if spectrum is None else spectrum.dna_dose_w_m2()
    return None if dose is None else dose / earth_dna_dose_w_m2()
