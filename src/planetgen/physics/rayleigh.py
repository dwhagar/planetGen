# planetgen/physics/rayleigh.py

"""
Rayleigh scattering of an atmosphere (GEN.417; docs/design/composition/atmosphere-composition-parts/02-gas-optics-uv-visible.md
section 3.1, docs/design/composition/radiation-model.md section 3).

The optical depth of the scattering by the air above the ground at wavelength `lambda` is

    tau_R(lambda) = sum_i N_i sigma_i(lambda),     N_i = p_i N_A / (m g)

with `N_i` the column of gas `i` (`AirColumn.number_column_cm2`) and `sigma_i` its cross section per molecule. The gas
table (`physics/gases.py`) gives each gas's cross section at 450, 550 and 650 nm; the wavelength law between and beyond
them is the one the research uses for any gas,

    sigma(lambda) = A lambda^-4 (1 + B / lambda^2)^2        (lambda in um, B in um2)

with `A` and `B` fitted to the gas's three values (the dispersion term `B` runs from 0.0024 um2 for helium to 0.0115 for
ammonia, and 0.0075 is the research's generic value). The fit is exact at 450, 550 and 650 nm and an extrapolation
outside 400 to 700 nm that holds from 200 nm up to 1100 nm for the gases whose first strong absorption lies below 200 nm
(the research: no use within about 80 nm of a strong absorption band, so SO2, H2S and O3 are scattering stand-ins
below 300 nm, where their absorption is a hundred times larger anyway).

Checks against the sources, in `tests/test_rayleigh.py`: dry Earth air at 550 nm gives 0.0970 (Bucholtz 1995 Table 4:
0.09721, 1 percent); Venus-like air (96.5 percent CO2) at 92 bar gives tau 16.4 at 550 nm (radiation-model.md section 11,
world C); 1 bar of 86/14 H2/He gives 0.089 at 550 nm (part 03 7.1, the giant reference level).

The cross sections are good to 2.5 to 4 percent for N2 and Ar, 7 for CO2 and 10 to 18 for CO, CH4 and N2O
(Sneep and Ubachs 2005, a measurement at 532 nm); the generic polarizability route is good to -9 to +6 percent.
`Rayleigh.state` and `Rayleigh.limit` of each row say which a gas is. Nothing here draws a value: it is a function of the
stored mixture and gravity.

Units: nm for wavelengths, cm2 per molecule for cross sections, cm-2 for columns.
"""


from planetgen.physics import gases

REFERENCE_WAVELENGTHS_NM = (450.0, 550.0, 650.0)
"""tuple: The wavelengths the gas table's cross sections are given at."""

MAXIMUM_DISPERSION_UM2 = 0.05
"""float: The largest dispersion term the fit may take; a table whose three values would need more is a power law instead."""

GENERIC_DISPERSION_UM2 = 0.0075
"""float: The research's generic dispersion term (part 02 section 3.1), used for the shape beyond the table's range when
the three values do not fix one."""

_FITS = {}


class NoRayleigh(ValueError):
    """A gas of the mixture has no Rayleigh cross section in the gas table."""


def _shape(wavelength_um, dispersion_um2):
    return wavelength_um ** -4 * (1.0 + dispersion_um2 / wavelength_um ** 2) ** 2


def _fit(formula):
    """`(A, B)` for `formula`: the amplitude (cm2 um4) and dispersion term (um2) through its three tabulated values."""
    cached = _FITS.get(formula)
    if cached is not None:
        return cached
    entry = gases.gas(formula).rayleigh
    if entry is None:
        raise NoRayleigh(f"{formula} has no Rayleigh cross section in the gas table (physics/data/gases.json)")
    short, middle, long_ = (entry.sigma_cm2[int(wavelength)] for wavelength in REFERENCE_WAVELENGTHS_NM)
    ratio = short / long_
    um = [wavelength / 1000.0 for wavelength in REFERENCE_WAVELENGTHS_NM]

    def mismatch(dispersion):
        return _shape(um[0], dispersion) / _shape(um[2], dispersion) - ratio

    low, high = 0.0, MAXIMUM_DISPERSION_UM2
    if mismatch(low) >= 0.0:           # a ratio the pure lambda^-4 law already reaches or exceeds
        dispersion = 0.0
    elif mismatch(high) <= 0.0:
        dispersion = high
    else:
        for _ in range(60):
            dispersion = 0.5 * (low + high)
            if mismatch(dispersion) < 0.0:
                low = dispersion
            else:
                high = dispersion
        dispersion = 0.5 * (low + high)
    amplitude = middle / _shape(um[1], dispersion)
    _FITS[formula] = (amplitude, dispersion)
    return _FITS[formula]


def dispersion_um2(formula):
    """The fitted dispersion term `B` of `formula`, um2."""
    return _fit(formula)[1]


def cross_section_cm2(formula, wavelength_nm):
    """
    The Rayleigh cross section of one molecule of `formula` at `wavelength_nm`, cm2.

    Raises:
        NoRayleigh: The gas table has no Rayleigh entry for `formula`.
    """
    amplitude, dispersion = _fit(formula)
    return amplitude * _shape(wavelength_nm / 1000.0, dispersion)


def optical_depth(column, wavelength_nm):
    """
    The vertical Rayleigh optical depth of `column` at `wavelength_nm`.

    Args:
        column (AirColumn): The air above the ground.
        wavelength_nm (float): The wavelength, nm.

    Raises:
        NoRayleigh: A gas of the mixture that holds more than one part in a million has no cross section.
    """
    total = 0.0
    for formula, pressure in column.partials_pa.items():
        if pressure / column.pressure_pa < 1e-6:
            continue
        total += column.number_column_cm2(formula) * cross_section_cm2(formula, wavelength_nm)
    return total


