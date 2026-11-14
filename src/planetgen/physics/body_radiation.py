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

from planetgen.physics import activity, aerosols, air_column, atmospheric_radiation, constants, radiation_field, stellar_spectrum, uv_cross_sections

DEFAULT_GROUND_ALBEDO = 0.15
"""float: The ground reflectance used until GEN.462 gives it by band (stand-in: the Earth-like land mean)."""
EARTH_OZONE_DU = 300.0
"""float: The ozone column of the reference Earth, DU (measured: 279 to 313 DU in the models and 300 DU in the record)."""
EARTH_BOND_ALBEDO = 0.306
"""float: The Earth's Bond albedo (measured), for the temperature of its ozone layer in the reference dose."""
SUN_EFFECTIVE_K = constants.SOLAR_TEMPERATURE_K
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


def _companion_sources(body, distance_au):
    """GEN.1442: `[(star, temperature_k, flux_w_m2)]` for the other stars of the body's system that light it (`radiation_field.sources_of`; none for a lone star
    or a body with no distance), the flux each one's time-averaged part of the field at the body."""
    star = getattr(body, "star", None)
    if star is None or not distance_au:
        return []
    found = []
    for member_no, other, factor in radiation_field.sources_of(star, distance_au):
        luminosity, temperature = getattr(other, "luminosity", None), getattr(other, "temperature", None)
        if member_no is None or not luminosity or not temperature:
            continue
        found.append((other, temperature, stellar_spectrum.bolometric_flux_w_m2(luminosity, 1.0) * factor))
    return found


def _spectrum_of(body, star, temperature, flux, ozone_flux, mu0, ground_albedo, column):
    """The `SurfaceSpectrum` of one source (`star`, its `temperature` and `flux` at the body); `ozone_flux` the whole field's, which warms the ozone layer."""
    share = getattr(star, "uv_c_share", None)
    stored = aerosols.of_body(body)
    centers = uv_cross_sections.centers(uv_cross_sections.standard_edges())
    haze, cloud = aerosols.layers(stored, centers)
    return atmospheric_radiation.surface_radiation(
        column, temperature_k=temperature, flux_w_m2=flux, uv_c_share=1.0 if share is None else share, uv_shares=activity.uv_shares_of(star),
        ozone_du=getattr(body, "ozone_column_du", None) or 0.0,
        ozone_temperature_k=atmospheric_radiation.ozone_layer_temperature_k(ozone_flux, getattr(body, "bond_albedo", None)), mu0=mu0,
        ground_albedo=ground_albedo,
        haze=haze, low_aerosol=cloud, cloud_cover=stored.cloud_cover)


def surface_spectrum(body, mu0=1.0, ground_albedo=DEFAULT_GROUND_ALBEDO, distance_au=None):
    """
    The `SurfaceSpectrum` at the ground of `body` for its host star at cosine `mu0` of the zenith angle (0.2 to 1), or `None` for a body with no air,
    no star or no gravity. The other stars of the system each have a spectrum of their own (`companion_spectra`); `uv_index` and `dna_dose_w_m2` add them.
    `distance_au` is the body's distance from its star when the body does not carry it.
    """
    column = air_column.column_of(body)
    where = distance_au or _star_distance_au(body)
    lit = _star_flux_w_m2(body, where)
    if column is None or lit is None:
        return None
    temperature, flux = lit
    ozone_flux = flux + sum(other[2] for other in _companion_sources(body, where))
    return _spectrum_of(body, body.star, temperature, flux, ozone_flux, mu0, ground_albedo, column)


def companion_spectra(body, mu0=1.0, ground_albedo=DEFAULT_GROUND_ALBEDO, distance_au=None):
    """GEN.1442: the `SurfaceSpectrum` of each other star of the body's system at its ground, a list (empty for a lone star or a body with no air)."""
    where = distance_au or _star_distance_au(body)
    others = _companion_sources(body, where)
    if not others:
        return []
    column = air_column.column_of(body)
    lit = _star_flux_w_m2(body, where)
    if column is None or lit is None:
        return []
    ozone_flux = lit[1] + sum(other[2] for other in others)
    return [_spectrum_of(body, star, temperature, flux, ozone_flux, mu0, ground_albedo, column) for star, temperature, flux in others]


def uv_index(body, mu0=1.0, ground_albedo=DEFAULT_GROUND_ALBEDO, distance_au=None):
    """The CIE UV index at the ground of `body` (every star of its system adds its own), 0 where there is no air to compute it in."""
    spectrum = surface_spectrum(body, mu0, ground_albedo, distance_au)
    if spectrum is None:
        return 0.0
    return spectrum.uv_index() + sum(other.uv_index() for other in companion_spectra(body, mu0, ground_albedo, distance_au))


def _dna_dose_under_air_w_m2(body, mu0, ground_albedo, distance_au):
    """The Setlow dose at the ground under the body's air from every star of its system, or `None` where there is no air or no star to compute it from."""
    spectrum = surface_spectrum(body, mu0, ground_albedo, distance_au)
    if spectrum is None:
        return None
    return spectrum.dna_dose_w_m2() + sum(other.dna_dose_w_m2() for other in companion_spectra(body, mu0, ground_albedo, distance_au))


def dna_dose_w_m2(body, mu0=1.0, ground_albedo=DEFAULT_GROUND_ALBEDO, distance_au=None):
    """The Setlow DNA-damage weighted irradiance at the ground of `body` (every star of its system adds its own), W/m2, 0 where there is no air to compute
    it in."""
    dose = _dna_dose_under_air_w_m2(body, mu0, ground_albedo, distance_au)
    return 0.0 if dose is None else dose


def bare_dna_dose_w_m2(body, mu0=1.0, distance_au=None):
    """The Setlow DNA-damage weighted irradiance of the star's unfiltered light on a horizontal surface, W/m2, or `None` for a body with no star or no
    distance. The ground of a body with no air."""
    where = distance_au or _star_distance_au(body)
    lit = _star_flux_w_m2(body, where)
    if lit is None:
        return None
    edges = uv_cross_sections.standard_edges()
    weight = atmospheric_radiation.setlow_weight(uv_cross_sections.centers(edges))
    total = 0.0
    for star, temperature, flux in [(body.star, *lit), *_companion_sources(body, where)]:
        share = getattr(star, "uv_c_share", None)
        top = stellar_spectrum.irradiance_w_m2(edges, temperature, flux, 1.0 if share is None else share, activity.uv_shares_of(star))
        total += float(mu0 * np.sum(top * weight))      # GEN.1442: every star of the system adds its own
    return total


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
        dose = _dna_dose_under_air_w_m2(body, mu0, ground_albedo, distance_au)
    return None if dose is None else dose / earth_dna_dose_w_m2()
