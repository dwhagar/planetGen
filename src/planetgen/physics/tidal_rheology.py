# planetgen/physics/tidal_rheology.py

"""
The tidal Love number and quality factor from temperature, melt and tide frequency (GEN.954)
=============================================================================================

Before GEN.955 a body had one `k2` and one `Q` whatever its temperature (rock Q 100, ice `k2/Q` 3e-4) and a `k2/Q` fitted to the one body that
shows a molten state (Io). A real mantle's response depends on how fast it is flexed against how fast it
relaxes: its `Q` is hundreds when it is cold and stiff, tens when it is warm, and falls to a minimum near the temperature where the
relaxation time equals the period of the tide; partly molten it is far weaker and its `k2` rises toward the fluid 3/2. This module is that
response, as a pure function of the body's size, its mantle's temperature and melt, and the tide's frequency. `interior_heat` reads it for
rock (GEN.955: `tidal_response` for the spin and the circularisation, the tide in the heat box); GEN.960 gives ice its own shell.

**The response.** A homogeneous sphere of mean density and radius `R` under its own gravity has the tidal Love number

    k2 = (3/2) / (1 + (19/2) mu~ / (rho g R)),        rho g R = 3 G M^2 / (4 pi R^4),

with `mu~` the *complex* rigidity, the inverse of the complex compliance `J` at the tide's angular frequency `chi`
(the elastic-viscoelastic correspondence principle; with the real rigidity this is `tidal_heating.rigid_k2`). The quality factor is
`Q = |k2| / (-Im k2)`. Two compliances:

* Maxwell (Henning, O'Connell and Sasselov 2009, ApJ 707:1000, arXiv:0912.1907, whose rheology it is): `J = 1/mu - i / (eta chi)`. A spring
  and a dashpot: dissipation peaks where `eta chi / mu` is 1 and falls off at both ends, but a cold mantle's `Q` runs to 1e5 or more,
  which no measurement of rock shows.
* Andrade (the default; Efroimsky 2012, ApJ 746:150 and Castillo-Rogez et al. 2011, Icarus 211:728, after Andrade 1910):
  `J = (1/mu) [1 + Gamma(1 + a) (i chi tau_A)^(-a)] - i / (eta chi)`, `tau_A = zeta eta / mu`. The extra term is the transient creep of
  grain boundaries; with `a = 0.3` it holds `Q` of a cold mantle near 100 to 1,000 and puts the Moon at the Q of 38 to 60 that the lunar
  laser ranging measures, which the Maxwell form cannot.

**The mantle.** Viscosity is Arrhenius `eta(T) = eta0 exp(Ea / R T)` with the body's own drawn `eta0` and `Ea` (`interior_heat.Interior`,
Guimond et al. 2022), so the heat flow of the one-box mantle and the tide read the same rheology. Melt lowers both viscosity and rigidity:

* `melt_fraction`: linear in temperature between the solidus and the liquidus (the one-box model's, `tidal_heating.SOLIDUS_K` and
  `LIQUIDUS_K`; they are passed in, so the wiring can give them a pressure dependence).
* Viscosity: `exp(-B phi)` with `B = 26` below the disaggregation fraction `CRITICAL_MELT_FRACTION` (Hirth and Kohlstedt 2003, for
  diffusion creep; Mei et al. 2002 find 26 to 31), then a straight line in `log10 eta` down to a silicate liquid's 1 Pa s at
  `phi = 1`: above the critical fraction the crystals no longer touch and the mush flows as a suspension.
* Cold ceiling: the Andrade `Q` of a mantle far below its solidus is capped at `COLD_Q_CEILING` (1,000), a mantle that is all melt included (a magma ocean keeps the loss of its solid lid and core): the lag of a stiff mantle goes to zero
  in the model but not in rock (the Maxwell form, kept as the reference, is not capped).
* Rigidity: `mu0 (1 - phi / phi_c)^2` down to a floor of `MELT_RIGIDITY_FLOOR` of `mu0` (the share the mush keeps when the matrix is gone;
  seismic mushes lose half or more of their shear modulus by a few tens of percent melt).

So a body between the solidus and `phi_c` is flexed at a frequency close to its relaxation time and dissipates most; at the heat-pipe melt
fraction of Io (0.3) this gives `k2/Q` of a few hundredths, the 0.015 that Io's orbital acceleration needs (Lainey et al. 2009), with no
number fitted to Io: the fit is gone, `interior_heat` balances this heating against the loss of the heat pipe.

**Ice** is the same equations with ice's rigidity (3.3 GPa), a lower activation energy (49 kJ/mol) and a viscosity pinned at the melting
point (`ICE`; Ojakangas and Stevenson 1989, Icarus 81:220 and McCarthy and Castillo-Rogez 2013 for the figures, read second-hand). A
homogeneous ice sphere is not Europa (its shell floats on an ocean and dissipates in a warm convecting layer): GEN.960 gives the shell its
own treatment on top of this form.

**What this does not do.** One temperature and one melt fraction stand for the whole mantle; the core, the lithosphere's cold lid and the
asthenosphere's thin warm layer, which carries most of Io's tide, are not separate layers; the rigidity does not fall with temperature
below the solidus; the frequency of a spin-orbit tide is the caller's (`2 |omega - n|` for despinning, `n` for an eccentricity tide).
These are stated stand-ins (`docs/design/research/gap-filling-plan.md` section 6.15).
"""

import cmath
import math
from dataclasses import dataclass

from planetgen.physics import constants

GAS_CONSTANT = 8.314462618
"""float: J/(mol K)."""
MAX_VISCOSITY_EXPONENT = 200.0
"""float: The largest `Ea / R T` the viscosity law is evaluated at (a viscosity of 1e87 times the prefactor): a mantle colder than that does not relax in any
tide, its response is the elastic one with the cold ceiling on `Q`, and the exponential would overflow at about 700 (a body at 50 K with a mantle drawn at
300 kJ/mol). Source tier T7 (a numerical cutoff), tweak level L0."""
LIQUID_VISCOSITY_PA_S = 1.0
"""float: A silicate liquid's viscosity, Pa s (komatiite to basalt melt, 0.1 to 10). Source tier T2, tweak level L0."""
MELT_VISCOSITY_FACTOR = 26.0
"""float: `B` of `exp(-B phi)`, the viscosity of a partly molten mantle against its solid one (Hirth and Kohlstedt 2003, diffusion creep,
26; Mei et al. 2002, 26 to 31; recalled). Source tier T2, tweak level L0."""
CRITICAL_MELT_FRACTION = 0.4
"""float: The melt share at which the crystals stop touching and the mush becomes a suspension (the rheological transition, 0.2 to 0.4 in
Arzi 1978, 0.4 in Henning et al. 2009; recalled). Source tier T3, tweak level L1 (the top of the range)."""
MELT_RIGIDITY_FLOOR = 1.0e-3
"""float: The share of its rigidity a mush keeps once the crystal matrix is gone. Source tier T6 (physics fixes the form: a suspension has no
static rigidity; the size is argued), tweak level L3 (chosen so the Io state's `k2` stays under 1.5)."""
COLD_Q_CEILING = 1000.0
"""float: The largest `Q` the Andrade form is allowed to give. A mantle far below its solidus is infinitely stiff in this model and its `Q` runs to 1e14, and a mantle that is all melt dissipates nothing; real cold rock
keeps a loss of its own (frictional and extrinsic, Q of 500 to a few thousand in the Earth's lithosphere and the Moon's cold interior; Dziewonski and Anderson 1981,
Nakamura 1982, recalled). Source tier T6, tweak level L3 (chosen so that a cold small planet is not frictionless)."""
ANDRADE_EXPONENT = 0.3
"""float: The Andrade exponent `a` (laboratory olivine 0.2 to 0.4; 0.3 in Jackson and Faul 2010 and Castillo-Rogez et al. 2011; recalled).
Source tier T2, tweak level L0."""
ANDRADE_TIMESCALE_RATIO = 1.0
"""float: `zeta = tau_A / tau_M`, the Andrade timescale over the Maxwell time (1 in Efroimsky 2012; Renaud and Henning 2018 explore 0.01 to
100 for the unknown). Source tier T6, tweak level L0."""


@dataclass(frozen=True)
class Material:
    """
    What a tidally flexed body is made of.

    Attributes:
        rigidity_pa (float): The solid's shear modulus, Pa.
        andrade_exponent (float): `a`.
        andrade_ratio (float): `zeta`.
    """

    rigidity_pa: float
    andrade_exponent: float = ANDRADE_EXPONENT
    andrade_ratio: float = ANDRADE_TIMESCALE_RATIO


ROCK = Material(7.0e10)
"""Material: A rocky mantle (the PREM average of 65 to 80 GPa, `tidal_heating.ROCK_RIGIDITY_PA`). Source tier T2, tweak level L0."""
ICE = Material(3.3e9)
"""Material: Ice I (shear modulus 3.3 to 3.5 GPa near its melting point; Gammon et al. 1983, recalled). Source tier T2, tweak level L0."""
ICE_ACTIVATION_J_MOL = 4.9e4
"""float: The activation energy of ice I's creep (49 kJ/mol; Goldsby and Kohlstedt 2001 for diffusion and grain-boundary sliding, McCarthy and
Castillo-Rogez 2013). Source tier T2, tweak level L0."""
ICE_MELTING_VISCOSITY_PA_S = 1.0e14
"""float: The viscosity of ice at its melting point, Pa s (1e13 to 1e15 across grain sizes; Ojakangas and Stevenson 1989 use 1e14). Source tier
T3, tweak level L0."""
ICE_MELTING_K = 273.0
"""float: The temperature that pins it, K."""


@dataclass(frozen=True)
class Response:
    """
    A body's tidal response at one frequency.

    Attributes:
        k2 (float): The magnitude of the complex Love number.
        q (float): The quality factor, `|k2| / (-Im k2)`; `inf` for a perfectly elastic or perfectly fluid body.
    """

    k2: float
    q: float

    @property
    def k2_over_q(self):
        """`k2/Q`, which sets the tide's power (`tidal_heating.peale_power_w`)."""
        return self.k2 / self.q if math.isfinite(self.q) else 0.0


def arrhenius_prefactor(reference_viscosity_pa_s, reference_k, activation_j_mol):
    """`eta0` of `eta(T) = eta0 exp(Ea / R T)` that gives `reference_viscosity_pa_s` at `reference_k` (for ice, pinned at its melting point)."""
    return reference_viscosity_pa_s * math.exp(-activation_j_mol / (GAS_CONSTANT * reference_k))


def melt_fraction(temperature_k, solidus_k, liquidus_k):
    """The melt share of the mantle at `temperature_k`: 0 at the solidus, 1 at the liquidus, linear between."""
    if liquidus_k <= solidus_k:
        raise ValueError(f"a liquidus of {liquidus_k} K must lie above the solidus of {solidus_k} K")
    return min(max((temperature_k - solidus_k) / (liquidus_k - solidus_k), 0.0), 1.0)


def viscosity_pa_s(temperature_k, eta0_pa_s, activation_j_mol, melt=0.0):
    """
    The viscosity at `temperature_k` and melt share `melt`, Pa s: Arrhenius for the solid, times `exp(-B melt)` up to `CRITICAL_MELT_FRACTION`,
    then a straight line in `log10 eta` from there to `LIQUID_VISCOSITY_PA_S` at a melt share of 1. The solid's value is taken at the
    temperature given, so a body above its solidus must be asked about with the temperature it has.
    """
    solid = eta0_pa_s * math.exp(min(activation_j_mol / (GAS_CONSTANT * temperature_k), MAX_VISCOSITY_EXPONENT))
    if melt <= 0.0:
        return solid
    if melt <= CRITICAL_MELT_FRACTION:
        return solid * math.exp(-MELT_VISCOSITY_FACTOR * melt)
    at_critical = solid * math.exp(-MELT_VISCOSITY_FACTOR * CRITICAL_MELT_FRACTION)
    share = (min(melt, 1.0) - CRITICAL_MELT_FRACTION) / (1.0 - CRITICAL_MELT_FRACTION)
    return 10.0 ** (math.log10(at_critical) + share * (math.log10(LIQUID_VISCOSITY_PA_S) - math.log10(at_critical)))


def rigidity_pa(melt, material=ROCK):
    """The shear modulus at melt share `melt`, Pa: `mu0 (1 - melt / CRITICAL_MELT_FRACTION)^2`, never under `MELT_RIGIDITY_FLOOR` of `mu0`."""
    weakened = max(1.0 - melt / CRITICAL_MELT_FRACTION, 0.0) ** 2
    return material.rigidity_pa * max(weakened, MELT_RIGIDITY_FLOOR)


def compliance(frequency_rad_s, rigidity, viscosity, material=ROCK, maxwell=False):
    """The complex compliance `J` of a body of `rigidity` Pa and `viscosity` Pa s at the tide's angular `frequency_rad_s` (the module's two forms;
    Andrade unless `maxwell`)."""
    unrelaxed = 1.0 / rigidity
    creep = 0.0
    if not maxwell:
        timescale = material.andrade_ratio * viscosity / rigidity
        creep = math.gamma(1.0 + material.andrade_exponent) * (1j * frequency_rad_s * timescale) ** (-material.andrade_exponent)
    return unrelaxed * (1.0 + creep) - 1j / (viscosity * frequency_rad_s)


def love_number(mass_kg, radius_m, frequency_rad_s, rigidity, viscosity, material=ROCK, maxwell=False):
    """The complex `k2` of a homogeneous sphere (the module's formula)."""
    pressure = 3.0 * constants.G * mass_kg ** 2 / (4.0 * math.pi * radius_m ** 4)
    complex_rigidity = 1.0 / compliance(frequency_rad_s, rigidity, viscosity, material, maxwell)
    return 1.5 / (1.0 + 9.5 * complex_rigidity / pressure)


def response(mass_kg, radius_m, frequency_rad_s, temperature_k, eta0_pa_s, activation_j_mol, melt=0.0, material=ROCK, maxwell=False):
    """
    The `Response` of a body of `mass_kg` and `radius_m` flexed at `frequency_rad_s`, with a mantle at `temperature_k` and melt share `melt`
    whose solid viscosity is `eta0_pa_s exp(activation_j_mol / R T)`. The same inputs the one-box mantle (`interior_heat`) holds.
    """
    if frequency_rad_s <= 0.0:
        raise ValueError(f"a tide frequency of {frequency_rad_s} rad/s does not flex a body")
    k2 = love_number(mass_kg, radius_m, frequency_rad_s, rigidity_pa(melt, material),
                     viscosity_pa_s(temperature_k, eta0_pa_s, activation_j_mol, melt), material, maxwell)
    lag = -k2.imag
    q = abs(k2) / lag if lag > 0.0 else math.inf
    return Response(abs(k2), q if maxwell else min(q, COLD_Q_CEILING))


def tidal_frequency_rad_s(orbital_rate, spin_rate=None):
    """The angular frequency of the tide an orbit raises: the orbital rate for an eccentricity tide on a body locked to it (`spin_rate` omitted),
    `2 |spin - orbital|` for the semi-diurnal tide of a spinning one (Efroimsky and Lainey 2007)."""
    if spin_rate is None:
        return orbital_rate
    return 2.0 * abs(spin_rate - orbital_rate)


def rock_response(mass_kg, radius_m, frequency_rad_s, temperature_k, eta0_pa_s, activation_j_mol, solidus_k, liquidus_k, maxwell=False):
    """`response` for a rocky mantle at `temperature_k`, its melt share worked out from the solidus and liquidus."""
    return response(mass_kg, radius_m, frequency_rad_s, temperature_k, eta0_pa_s, activation_j_mol,
                    melt_fraction(temperature_k, solidus_k, liquidus_k), ROCK, maxwell)


def ice_response(mass_kg, radius_m, frequency_rad_s, temperature_k):
    """`response` for an ice body at `temperature_k` (no melt share: the ice shell's own treatment is GEN.960), the viscosity pinned at
    `ICE_MELTING_VISCOSITY_PA_S` at `ICE_MELTING_K`."""
    eta0 = arrhenius_prefactor(ICE_MELTING_VISCOSITY_PA_S, ICE_MELTING_K, ICE_ACTIVATION_J_MOL)
    return response(mass_kg, radius_m, frequency_rad_s, temperature_k, eta0, ICE_ACTIVATION_J_MOL, 0.0, ICE)
