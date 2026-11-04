# planetgen/physics/atmospheric_radiation.py

"""
Atmospheric radiation: what the air does to the star's light, band by band (GEN.443; docs/design/radiation-model.md R6, R7,
R17 and sections 2, 3 and 7).

`surface_radiation` takes the air above a place (`physics/air_column.AirColumn`, from the stored mixture and gravity),
the star's light at the planet (`physics/stellar_spectrum`), the ozone column (GEN.413), the aerosol layers (GEN.414) and the ground
albedo, and returns the spectrum that reaches the ground on a grid of 399 cells from 100 nm to 4 micron, with the sums over
the bands of R6 (UV-C, UV-B, UV-A, visible, near infrared), the CIE erythemal UV index and the Setlow DNA-damage dose of section 7.

**The layers**, top to bottom, are solved by `physics/two_stream` (delta-Eddington, layers added with the light bouncing between
them, a Lambertian ground below):

1. the air above the ozone layer: the mixture's Rayleigh scattering (`physics/rayleigh`) and the absorption of every gas the
   gas table has ultraviolet cross sections for (`physics/uv_cross_sections`, band means of [R1873, R2040]:
   O2, CO2, H2O, SO2, H2S, NH3, CH4 ...), a share `p_ozone / P` of the whole column;
2. the ozone (a pure absorber): the column in Dobson units times the high-resolution cross section at the layer's temperature;
3. the rest of the air, with the haze that is mixed with it;
4. the cloud and dust, near the ground, taken as one layer (`cloud_cover` of the sky overcast: the result is the area-weighted mix
   of the clear and the overcast result);
5. the ground, of albedo `ground_albedo`.

**The ozone layer's pressure** is the peak of the Chapman production layer, `6.3 kPa` per `PAL` where `PAL = p_O2 / 21.2 kPa`
(part 04 3.3: 63 hPa at 1 PAL, 630 at 0.1, 6300 at 0.01; an upper bound on the depth, as the part says). A caller that knows
the pressure of the ozone layer passes `ozone_pressure_pa`.

**Limits**, stated rather than hidden: the near infrared above 1100 nm is Rayleigh, aerosol and ground only (the H2O, CO2 and CH4
bands wait for the correlated-k grids of R2, so the near-infrared share is an upper bound on what the ground gets); a gas mixed
through the air that absorbs strongly at one height (O2 in the Schumann-Runge bands) is a band mean, which understates the
transmission of a thick column (the gas table says so for O2); the diffuse field is a two-stream field, within 1.3 percent of
a Monte Carlo run for `mu0` of 0.5 and 1 and 3.5 percent low at 0.2; there is no spherical geometry (use `mu0` of at least 0.2) and
the cross sections do not change with the temperature of the air except ozone's.

Units: nm, W/m2, Dobson units (1 DU = 2.687e16 molecules/cm2), Pa.
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import constants, rayleigh, stellar_spectrum, two_stream, uv_cross_sections

OZONE_PEAK_PA_PER_PAL = 6300.0
"""float: The pressure of the Chapman peak at 1 PAL of oxygen, Pa (part 04 3.3: 63 hPa)."""
EARTH_PO2_PA = 21200.0
"""float: 1 PAL of oxygen, Pa (radiation-model.md R8: 21.2 kPa)."""

# the CIE erythemal action spectrum (ISO 17166 / CIE 2014, part 02 2.6)
UVI_PER_W_M2 = 40.0
"""float: The UV index is 40 m2/W times the erythemally weighted irradiance."""

SETLOW_NM = (256.0, 260.0, 266.0, 268.0, 270.0, 276.0, 280.0, 286.0, 290.0, 294.0, 298.0, 300.0, 304.0, 310.0, 314.0, 320.0, 328.0,
             340.0, 350.0, 360.0)
SETLOW_WEIGHT = (1.07, 1.12, 1.23, 1.23, 1.13, 0.933, 0.791, 0.530, 0.309, 0.163, 0.0552, 0.0333, 0.00901, 0.00115, 2.5e-4, 4.3e-5,
                 1.1e-5, 4.6e-6, 2.6e-6, 1.8e-6)
"""The [R4033, R4035] DNA-damage action spectrum, per quantum (TUV-x data; the file says it is normalised at 300 nm but stores 0.0333 there,
reported not resolved: part 02 2.6). Multiplied by `lambda / 300` for an energy irradiance."""

CLOUD_OMEGA = 1.0 - 1e-6
CLOUD_ASYMMETRY = 0.85
DUST_ASYMMETRY = 0.673


def erythemal_weight(wavelength_nm):
    """The CIE erythemal action spectrum `S_er`: 1 for 250 to 298 nm, `10^(0.094 (298 - lambda))` to 328 nm, `10^(0.015 (140 - lambda))`
    to 400 nm, 0 elsewhere (part 02 2.6; the NOAA text has (139 - lambda), a typo)."""
    nm = np.asarray(wavelength_nm, dtype=float)
    return np.where((nm >= 250.0) & (nm <= 298.0), 1.0,
                    np.where((nm > 298.0) & (nm <= 328.0), 10.0 ** (0.094 * (298.0 - nm)),
                             np.where((nm > 328.0) & (nm <= 400.0), 10.0 ** (0.015 * (140.0 - nm)), 0.0)))


def setlow_weight(wavelength_nm):
    """The Setlow DNA-damage weight per unit energy: log-linear between the tabulated points, held at 1.07 below 256 nm and
    continued at the last slope above 360 nm (negligible there), times `lambda / 300`."""
    nm = np.asarray(wavelength_nm, dtype=float)
    log_weight = np.interp(nm, SETLOW_NM, np.log(SETLOW_WEIGHT))
    slope = (math.log(SETLOW_WEIGHT[-1]) - math.log(SETLOW_WEIGHT[-2])) / (SETLOW_NM[-1] - SETLOW_NM[-2])
    log_weight = np.where(nm > SETLOW_NM[-1], math.log(SETLOW_WEIGHT[-1]) + slope * (nm - SETLOW_NM[-1]), log_weight)
    return np.where(nm <= 400.0, np.exp(log_weight) * nm / 300.0, 0.0)


@dataclass(frozen=True)
class AerosolLayer:
    """
    One scattering layer of aerosol or cloud on the model's grid (the product of GEN.414's `physics/aerosols`).

    Attributes:
        tau (ndarray): Optical depth at each cell (scattering plus absorption).
        omega (ndarray): Single-scattering albedo at each cell.
        g (ndarray): Asymmetry factor at each cell.
    """

    tau: np.ndarray
    omega: np.ndarray
    g: np.ndarray


def grey_layer(tau_550, omega, g, centres_nm, exponent=0.0):
    """An `AerosolLayer` whose optical depth is `tau_550 (550 nm / lambda)^exponent` (0 for a grey cloud), with constant albedo and asymmetry."""
    tau = tau_550 * (550.0 / np.asarray(centres_nm, dtype=float)) ** exponent
    return AerosolLayer(tau, np.full_like(tau, omega), np.full_like(tau, g))


@dataclass(frozen=True)
class SurfaceSpectrum:
    """
    The light that reaches the ground.

    Attributes:
        edges_nm (ndarray): The cell edges.
        top_w_m2 (ndarray): The star's flux through a surface facing it, per cell, at the planet (normal incidence).
        ground_w_m2 (ndarray): The flux on a horizontal surface at the ground, per cell (direct plus diffuse), W/m2.
        direct_w_m2 (ndarray): The part of it that comes straight from the star.
        reflected_w_m2 (ndarray): What leaves the top of the atmosphere, per cell, per unit horizontal incident flux times the same.
        mu0 (float): The cosine of the zenith angle.
        ozone_du (float): The ozone column used.
    """

    edges_nm: np.ndarray
    top_w_m2: np.ndarray
    ground_w_m2: np.ndarray
    direct_w_m2: np.ndarray
    outgoing_w_m2: np.ndarray
    mu0: float
    ozone_du: float

    @property
    def centres_nm(self):
        return 0.5 * (self.edges_nm[:-1] + self.edges_nm[1:])

    def _band(self, values, band):
        low, high = stellar_spectrum.BAND_EDGES_NM[band]
        return float(np.sum(values[(self.centres_nm >= low) & (self.centres_nm < high)]))

    def band_top_w_m2(self, band):
        """The star's flux on a horizontal surface at the top of the atmosphere in `band`, W/m2."""
        return self.mu0 * self._band(self.top_w_m2, band)

    def band_ground_w_m2(self, band):
        """The flux on the ground in `band`, W/m2."""
        return self._band(self.ground_w_m2, band)

    def band_transmission(self, band):
        """The ground flux over the top flux in `band` (0 when the star has none there)."""
        top = self.band_top_w_m2(band)
        return self.band_ground_w_m2(band) / top if top > 0.0 else 0.0

    def band_fractions(self):
        """The transmission of each band of R6 below the thermal infrared, as a dict."""
        return {band: self.band_transmission(band) for band in ("uv_c", "uv_b", "uv_a", "visible", "near_infrared")}

    def erythemal_w_m2(self):
        """The erythemally weighted irradiance on a horizontal surface at the ground, W/m2 (CIE)."""
        return float(np.sum(self.ground_w_m2 * erythemal_weight(self.centres_nm)))

    def uv_index(self):
        """The CIE UV index at the ground: 40 m2/W times the erythemally weighted irradiance."""
        return UVI_PER_W_M2 * self.erythemal_w_m2()

    def top_uv_index(self):
        """The UV index of the star's light on a horizontal surface with no atmosphere."""
        return UVI_PER_W_M2 * self.mu0 * float(np.sum(self.top_w_m2 * erythemal_weight(self.centres_nm)))

    def dna_dose_w_m2(self):
        """The Setlow DNA-damage weighted irradiance at the ground, W/m2 (energy basis: weight times lambda / 300)."""
        return float(np.sum(self.ground_w_m2 * setlow_weight(self.centres_nm)))

    def top_dna_dose_w_m2(self):
        """The Setlow weighted irradiance with no atmosphere, W/m2."""
        return self.mu0 * float(np.sum(self.top_w_m2 * setlow_weight(self.centres_nm)))


def ozone_peak_pressure_pa(column, ozone_du=0.0):
    """The pressure of the ozone layer, Pa: the Chapman peak at `6.3 kPa / PAL` for the column's oxygen, never below 10 Pa or deeper than
    nine tenths of the surface pressure. A column with no oxygen gets nine tenths of its surface pressure (the ozone of a column without
    O2 is all made at the bottom)."""
    p_o2 = column.partials_pa.get("O2", 0.0)
    pal = p_o2 / EARTH_PO2_PA
    if pal <= 0.0:
        return 0.9 * column.pressure_pa
    return float(min(max(OZONE_PEAK_PA_PER_PAL / pal, 10.0), 0.9 * column.pressure_pa))


def mix_layers(a, b):
    """Combines two scattering species in one slab: returns (tau, omega, g) arrays from `(tau, omega, g)` triples a and b."""
    tau = a[0] + b[0]
    scattering = a[0] * a[1] + b[0] * b[1]
    safe = np.where(tau > 0.0, tau, 1.0)
    omega = np.where(tau > 0.0, scattering / safe, 1.0)
    asymmetry = np.where(scattering > 0.0, (a[0] * a[1] * a[2] + b[0] * b[1] * b[2]) / np.where(scattering > 0.0, scattering, 1.0), 0.0)
    return tau, omega, asymmetry


def gas_layers(column, centres, edges):
    """The Rayleigh and absorption optical depth of the whole column at each cell, and the absorption alone."""
    tau_rayleigh = np.zeros_like(centres)
    for formula, pressure in column.partials_pa.items():
        if pressure / column.pressure_pa < 1e-6:
            continue
        tau_rayleigh = tau_rayleigh + column.number_column_cm2(formula) * np.array(
            [rayleigh.cross_section_cm2(formula, nm) for nm in centres])
    tau_absorption = np.zeros_like(centres)
    for formula, pressure in column.partials_pa.items():
        if formula == "O3" or pressure / column.pressure_pa < 1e-12:
            continue
        tau_absorption = tau_absorption + column.number_column_cm2(formula) * uv_cross_sections.cross_section_cm2(formula, edges)
    return tau_rayleigh, tau_absorption


DEFAULT_OZONE_TEMPERATURE_K = uv_cross_sections.DEFAULT_OZONE_TEMPERATURE_K
DEFAULT_BOND_ALBEDO = 0.3
"""float: The Bond albedo of a body that has none stored, for its ozone layer's temperature (Earth's is 0.306)."""
SKIN_TEMPERATURE_FACTOR = 2.0 ** -0.25
"""float: The temperature of the thin top of a grey atmosphere over its equilibrium temperature, `2^(-1/4)` (about 0.84)."""


def ozone_layer_temperature_k(flux_w_m2, bond_albedo=None):
    """
    The temperature of a body's ozone layer, K (GEN.748): its skin temperature, the equilibrium temperature `(F (1 - A) / 4 sigma)^(1/4)`
    times `SKIN_TEMPERATURE_FACTOR`. Earth's is 214 K (its ozone layer sits at 210 to 230 K), a body at the same flux and albedo with
    a thicker or thinner air has the same one. The cross sections are held at the ends of their tables (218 to 295 K in the ultraviolet,
    220 to 298 K in the Chappuis band).
    """
    albedo = DEFAULT_BOND_ALBEDO if bond_albedo is None else min(max(bond_albedo, 0.0), 1.0)
    equilibrium = (flux_w_m2 * (1.0 - albedo) / (4.0 * constants.STEFAN_BOLTZMANN_CONSTANT)) ** 0.25
    return equilibrium * SKIN_TEMPERATURE_FACTOR


def surface_radiation(column, *, temperature_k, flux_w_m2, uv_c_share=1.0, ozone_du=0.0, ozone_temperature_k=None,
                      ozone_pressure_pa=None, mu0=1.0, ground_albedo=0.15, haze=None, low_aerosol=None, cloud_cover=1.0,
                      edges_nm=None, uv_shares=None):
    """
    The spectrum that reaches the ground below `column`.

    Args:
        column (AirColumn): The air above the place.
        temperature_k (float): The star's effective temperature, K.
        flux_w_m2 (float): The star's bolometric flux at the planet, W/m2.
        uv_c_share (float, optional): The star's UV-C share of its blackbody (R17; 1 is the blackbody).
        uv_shares (dict, optional): The star's `uv_c`, `uv_b` and `uv_a` shares of its blackbody (`stellar_library.uv_shares`); when given, `uv_c_share` is not used.
        ozone_du (float, optional): The ozone column, Dobson units (GEN.413).
        ozone_temperature_k (float, optional): The temperature of the ozone layer, K (`ozone_layer_temperature_k`; 220 when `None`).
        ozone_pressure_pa (float, optional): The pressure of the ozone layer, Pa; the Chapman peak of `ozone_peak_pressure_pa` when `None`.
        mu0 (float, optional): The cosine of the star's zenith angle, 0.2 to 1 (1 is overhead).
        ground_albedo (float or ndarray, optional): The ground's reflectance, a number or one value per cell.
        haze (AerosolLayer, optional): Aerosol mixed with the lower air (GEN.414).
        low_aerosol (AerosolLayer, optional): Cloud and dust near the ground (GEN.414).
        cloud_cover (float, optional): The share of the sky `low_aerosol` covers; the rest is clear.
        edges_nm (ndarray, optional): The grid; `uv_cross_sections.standard_edges()` when `None`.

    Returns:
        SurfaceSpectrum
    """
    if not 0.2 <= mu0 <= 1.0:
        raise ValueError("the two-stream model is made for the star between 0.2 and 1 in cosine of the zenith angle")
    edges = uv_cross_sections.standard_edges() if edges_nm is None else np.asarray(edges_nm, dtype=float)
    centres = uv_cross_sections.centres(edges)
    top = stellar_spectrum.irradiance_w_m2(edges, temperature_k, flux_w_m2, uv_c_share, uv_shares)
    tau_rayleigh, tau_absorption = gas_layers(column, centres, edges)

    pressure_ozone = ozone_peak_pressure_pa(column, ozone_du) if ozone_pressure_pa is None else ozone_pressure_pa
    above = min(max(pressure_ozone / column.pressure_pa, 0.0), 1.0)
    tau_ozone = np.zeros_like(centres)
    if ozone_du > 0.0:
        sigma = uv_cross_sections.ozone_cross_section_cm2(edges, ozone_temperature_k)
        tau_ozone = ozone_du * uv_cross_sections.DOBSON_CM2 * sigma

    def gas_slab(share):
        tau = share * (tau_rayleigh + tau_absorption)
        omega = np.where(tau > 0.0, share * tau_rayleigh / np.where(tau > 0.0, tau, 1.0), 1.0)
        return tau, omega, np.zeros_like(tau)

    upper = two_stream.scattering_layer(*gas_slab(above), mu0)
    ozone = two_stream.absorbing_layer(tau_ozone, mu0)
    lower_slab = gas_slab(1.0 - above)
    if haze is not None:
        lower_slab = mix_layers(lower_slab, (haze.tau, haze.omega, haze.g))
    lower = two_stream.scattering_layer(*lower_slab, mu0)
    albedo = np.broadcast_to(np.asarray(ground_albedo, dtype=float), centres.shape)

    def solve(with_low_layer):
        layers = [upper, ozone, lower]
        if with_low_layer:
            layers.append(two_stream.scattering_layer(low_aerosol.tau, low_aerosol.omega, low_aerosol.g, mu0))
        stack = two_stream.stack_layers(layers)
        arriving, outgoing = two_stream.ground_flux(stack, albedo)
        direct = stack.beam_transmission
        return arriving, outgoing, direct

    arriving, outgoing, direct = solve(False)
    if low_aerosol is not None and cloud_cover > 0.0:
        arriving_cloud, outgoing_cloud, direct_cloud = solve(True)
        cover = min(max(cloud_cover, 0.0), 1.0)
        arriving = (1.0 - cover) * arriving + cover * arriving_cloud
        outgoing = (1.0 - cover) * outgoing + cover * outgoing_cloud
        direct = (1.0 - cover) * direct + cover * direct_cloud
    horizontal = mu0 * top
    return SurfaceSpectrum(edges, top, horizontal * arriving, horizontal * direct, horizontal * outgoing, mu0, ozone_du)
