"""
The core's own thermal history (GEN.838), run beside the mantle box of `interior_heat`.

`interior_heat.heat` follows one mantle temperature (the core sits at it and only adds its heat capacity). The dynamo needs more: the core
cools by the heat it loses across the **core-mantle boundary**, and what that boundary layer carries is set by the *difference* between the core's
temperature and the base of the mantle, not by the mantle's temperature alone. This module is that second state, a core temperature `T_c` at the
top of the core, driven by the mantle box's temperature history (the mantle box is not changed by it: its stored heat flow and Urey ratio are
what they were, a body with no core keeps no core state).

1. **Geometry.** The core is one uniform shell of the body's metal mass inside the core radius the mantle box already uses; the mantle is one
   uniform shell outside it (the shells of `interior_heat.binding_energy_j`). Pressures are those of the shells times `PRESSURE_FACTOR` (the
   compression a uniform body leaves out; Earth's core-mantle boundary at 136 GPa and its centre at 364 GPa).
2. **Boundary layer.** `Q = k_m A (T_c - T_b) / delta`, `T_b` the mantle's adiabat carried down to the boundary (`MANTLE_ADIABAT_K_GPA` per GPa
   above the box's temperature), `delta` the layer's thickness from its critical Rayleigh number at the viscosity of the layer's mean
   temperature with an activation volume (`ACTIVATION_VOLUME_M3_MOL`; the box's viscosity law is the upper mantle's, the lower mantle is stiffer).
3. **Cooling.** `C dT_c/dt = -Q`, `C` the core's heat capacity plus, once there is an inner core, the latent heat and gravitational energy of the iron
   that freezes onto it per kelvin of cooling (`LATENT_HEAT_J_KG`, `GRAVITATIONAL_ENERGY_J_KG`).
4. **Inner core.** The core's adiabat (`CORE_ADIABAT_K_GPA` per GPa below `T_c`) meets the iron liquidus (`T_0 sqrt(1 + P / A)`, anchored on Earth's inner
   core boundary, 5,400 K at 329 GPa) where the inner core stops: the centre freezes first, and the inner core's radius follows `T_c`.
5. **Dynamo.** The core conducts down its own adiabat, `Q_ad = k_c A (dT/dP)_ad rho g`; the convecting power left over, `(1 - T_b / T_c)(Q - Q_ad)`,
   plus the gravitational energy of the freezing iron, drives the dynamo (`Core.dynamo_w`). The dynamo runs while that power is above zero and the
   core still has a liquid shell. `Core.dynamo_ended_gy` is the age it last stopped (`None` while it runs).

Every number that is not a published value is a stand-in in the gap registry (docs/design/research/gap-filling-plan.md section 6.136).
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants

GAS_CONSTANT = 8.314462618
"""float: J/(mol K)."""
GY_S = 1.0e9 * constants.SECONDS_PER_YEAR
"""float: Seconds in a gigayear."""
MAX_EXPONENT = 600.0
"""float: The largest exponent of the viscosity law (`exp` overflows at 709): a box this cold has a viscosity no heat crosses (a numerical bound)."""

PRESSURE_FACTOR = 1.16
"""float: STAND-IN CC1: the compression a uniform shell leaves out, as one factor on the shells' pressures."""
MANTLE_ADIABAT_K_GPA = 8.0
"""float: STAND-IN CC2: K per GPa the mantle's adiabat gains from the box's temperature to the core-mantle boundary."""
CORE_ADIABAT_K_GPA = 6.74
"""float: STAND-IN CC3: K per GPa the core's adiabat gains below the top of the core (4,100 K at 136 GPa to 5,400 K at 329 GPa)."""
LIQUIDUS_ZERO_K = 1811.0
"""float: The melting temperature of iron at 1 bar, K."""
LIQUIDUS_ANCHOR = (329.0, 5400.0)
"""tuple: STAND-IN CC4: `(P in GPa, T in K)` of Earth's inner core boundary, which sets the liquidus's pressure scale."""
CORE_CONDUCTIVITY = 40.0
"""float: STAND-IN CC5: W/(m K), the core's thermal conductivity (the literature spans 30 to 100)."""
LATENT_HEAT_J_KG = 2.5e5
"""float: STAND-IN CC6: J/kg, the latent heat of the iron that freezes onto the inner core."""
GRAVITATIONAL_ENERGY_J_KG = 2.5e5
"""float: STAND-IN CC6: J/kg, the gravitational energy the light elements left behind release as it grows."""
ACTIVATION_VOLUME_M3_MOL = 4.7e-6
"""float: STAND-IN CC7: m3/mol, the activation volume that stiffens the mantle's lowest layer."""
PRESSURE_SCALE_GPA = 100.0
"""float: STAND-IN CC7: GPa over which the activation volume falls to half, so the pressure term of the activation enthalpy grows as a logarithm."""
INITIAL_EXCESS_K = 1000.0
"""float: STAND-IN CC8: K the core starts hotter than the base of the mantle (the sinking iron's heat)."""

MIN_LIQUID_SHELL = 0.25
"""float: STAND-IN CC9: the liquid shell's thickness over the core's radius under which no dynamo runs (a thin shell cannot sustain the flow: Mars)."""
EARTH_CORE_RADIUS_M = 3.632e6
EARTH_CORE_DENSITY_KG_M3 = 9671.0
EARTH_DYNAMO_W = 2.4e12
"""float: STAND-IN CC10: Earth's core as this model makes it (the radius and density of its uniform core, and the dynamo power it ends with at
4.568 Gy, 2.4 TW): the unit of `Core.dipole_scale`."""

BOUNDARY_RAYLEIGH = 450.0
"""float: The critical Rayleigh number of the lower boundary layer: the mantle box's ([R4339], second-hand)."""
MANTLE_CONDUCTIVITY = 4.0
MANTLE_EXPANSIVITY = 2.5e-5
MANTLE_DIFFUSIVITY = 1.0e-6
CORE_SPECIFIC_HEAT = 840.0
"""float: The mantle box's table ([R4339] via [R1716: Table 2]): W/(m K), 1/K, m2/s, J/(kg K)."""


@dataclass(frozen=True)
class Core:
    """
    A core at an age.

    Attributes:
        temperature_k (float): At the top of the core.
        boundary_k (float): The mantle's, at the other side of the boundary layer.
        cmb_w (float): The heat the core loses across the core-mantle boundary.
        adiabatic_w (float): The share of it the core conducts down its own adiabat.
        inner_fraction (float): The inner core's radius over the core's (0 while the core is all liquid, 1 when it is all solid).
        dynamo_w (float): The power left for the dynamo, 0 while it does not run.
        dynamo_ended_gy (float or None): The age the dynamo last stopped; `None` while it runs.
        radius_m, density_kg_m3 (float): The uniform core's.
    """

    temperature_k: float
    boundary_k: float
    cmb_w: float
    adiabatic_w: float
    inner_fraction: float
    dynamo_w: float
    dynamo_ended_gy: float | None
    radius_m: float
    density_kg_m3: float

    @property
    def running(self):
        """Whether the dynamo runs."""
        return self.dynamo_w > 0.0

    @property
    def dipole_scale(self):
        """The dipole moment this dynamo makes over Earth's: the buoyancy-flux scaling, a field strength of `(rho^1/6)(P / V)^(1/3)` over a core's
        volume `V`, times the core's surface `R^2` (the moment is the field times `R^3`): `(R / R_e)^2 (P / P_e)^(1/3) (rho / rho_e)^(1/6)`. 0 while
        the dynamo does not run."""
        if not self.running:
            return 0.0
        return ((self.radius_m / EARTH_CORE_RADIUS_M) ** 2 * (self.dynamo_w / EARTH_DYNAMO_W) ** (1.0 / 3.0)
                * (self.density_kg_m3 / EARTH_CORE_DENSITY_KG_M3) ** (1.0 / 6.0))


class Stepper:
    """One body's core, advanced step by step beside the mantle box (`advance`); `start` first."""

    def __init__(self, interior, metal_kg, silicate_kg, core_m, dry_m):
        self.eta0_pa_s = interior.eta0_pa_s
        self.activation_j_mol = interior.activation_kj_mol * 1000.0
        self.metal_kg, self.core_m = metal_kg, core_m
        volume = 4.0 / 3.0 * math.pi * core_m ** 3
        self.density = metal_kg / volume
        mantle_density = silicate_kg / (4.0 / 3.0 * math.pi * (dry_m ** 3 - core_m ** 3))
        self.mantle_density = mantle_density
        self.gravity = constants.G * metal_kg / core_m ** 2
        self.area = 4.0 * math.pi * core_m ** 2
        enclosed = metal_kg - 4.0 / 3.0 * math.pi * mantle_density * core_m ** 3
        pressure_pa = constants.G * mantle_density * (enclosed * (1.0 / core_m - 1.0 / dry_m)
                                                      + 4.0 / 3.0 * math.pi * mantle_density * (dry_m ** 2 - core_m ** 2) / 2.0)
        self.cmb_gpa = PRESSURE_FACTOR * pressure_pa / 1.0e9
        self.core_slope_gpa_m2 = PRESSURE_FACTOR * 2.0 * math.pi / 3.0 * constants.G * self.density ** 2 / 1.0e9
        self.centre_gpa = self.cmb_gpa + self.core_slope_gpa_m2 * core_m ** 2
        self.liquidus_scale_gpa = LIQUIDUS_ANCHOR[0] / ((LIQUIDUS_ANCHOR[1] / LIQUIDUS_ZERO_K) ** 2 - 1.0)
        self.capacity = metal_kg * CORE_SPECIFIC_HEAT
        self.pressure_enthalpy_j_mol = ACTIVATION_VOLUME_M3_MOL * PRESSURE_SCALE_GPA * 1.0e9 * math.log(1.0 + self.cmb_gpa / PRESSURE_SCALE_GPA)
        self.adiabatic_w = CORE_CONDUCTIVITY * self.area * CORE_ADIABAT_K_GPA / 1.0e9 * self.density * self.gravity
        self.temperature_k = 0.0
        self.ended_gy = None
        self._running = True

    def boundary_k(self, mantle_k):
        """The base of the mantle for a box at `mantle_k`, K."""
        return mantle_k + MANTLE_ADIABAT_K_GPA * self.cmb_gpa

    def liquidus_k(self, pressure_gpa):
        """The iron liquidus at `pressure_gpa`, K."""
        return LIQUIDUS_ZERO_K * math.sqrt(1.0 + pressure_gpa / self.liquidus_scale_gpa)

    def _inner(self, temperature_k):
        """`(radius, crossing)` of the inner core when the top of the core is at `temperature_k`: the radius is 0 while the whole core is above its
        liquidus and the core's radius when none of it is; the crossing is the pressure (GPa) where the adiabat meets the liquidus, `None` in those two cases."""
        k = CORE_ADIABAT_K_GPA

        def gap(pressure_gpa):
            return self.liquidus_k(pressure_gpa) - (temperature_k + k * (pressure_gpa - self.cmb_gpa))
        if gap(self.centre_gpa) < 0.0:
            return 0.0, None
        if gap(self.cmb_gpa) >= 0.0:
            return self.core_m, None
        # The liquidus `T0 sqrt(1 + P / s)` meets the adiabat `u + k P` (u = T - k P_cmb) where `(u + k P)^2 = T0^2 (1 + P / s)`: a quadratic in P. The gap
        # is concave, closed at one end of the core and open at the other, so exactly one root lies between them.
        u = temperature_k - k * self.cmb_gpa
        b = 2.0 * u * k - LIQUIDUS_ZERO_K ** 2 / self.liquidus_scale_gpa
        c = u * u - LIQUIDUS_ZERO_K ** 2
        root = math.sqrt(max(b * b - 4.0 * k * k * c, 0.0))
        crossing = min(((-b + sign * root) / (2.0 * k * k) for sign in (1.0, -1.0)),
                       key=lambda pressure: abs(min(max(pressure, self.cmb_gpa), self.centre_gpa) - pressure))
        crossing = min(max(crossing, self.cmb_gpa), self.centre_gpa)
        return math.sqrt(max(self.core_m ** 2 - (crossing - self.cmb_gpa) / self.core_slope_gpa_m2, 0.0)), crossing

    def inner_radius_m(self, temperature_k):
        """The inner core's radius when the top of the core is at `temperature_k`: 0 while the whole core is above its liquidus, the core's
        radius when none of it is."""
        return self._inner(temperature_k)[0]

    def _growth_per_k(self, temperature_k):
        """Inner core mass gained per kelvin of cooling, kg/K: the derivative of the inner core's mass `4/3 pi rho r^3` along the crossing, where the
        gap `liquidus - adiabat` is nil (`dP/dT = 1 / (liquidus' - adiabat')`) and the radius follows the pressure through the core's own gravity."""
        radius, crossing = self._inner(temperature_k)
        if crossing is None:
            return 0.0
        slope = (self.liquidus_k(crossing) / (2.0 * (self.liquidus_scale_gpa + crossing))) - CORE_ADIABAT_K_GPA
        if slope <= 0.0:
            return 0.0
        return 2.0 * math.pi * self.density * radius / (self.core_slope_gpa_m2 * slope)

    def _flux(self, temperature_k, mantle_k):
        """The heat across the boundary, W (negative while the mantle is the hotter)."""
        boundary_k = self.boundary_k(mantle_k)
        drop = temperature_k - boundary_k
        mean_k = max(0.5 * (temperature_k + boundary_k), 1.0)
        viscosity = self.eta0_pa_s * math.exp(min((self.activation_j_mol + self.pressure_enthalpy_j_mol)
                                                  / (GAS_CONSTANT * mean_k), MAX_EXPONENT))      # a very cold box is a solid one: no heat crosses, not an overflow
        gap = abs(drop)
        if gap < 1.0e-6:
            return 0.0
        layer = (BOUNDARY_RAYLEIGH * MANTLE_DIFFUSIVITY * viscosity / (MANTLE_EXPANSIVITY * self.mantle_density * self.gravity * gap)) ** (1.0 / 3.0)
        return math.copysign(MANTLE_CONDUCTIVITY * self.area * gap / layer, drop)

    def _cooling(self, temperature_k, mantle_k):
        """`(dT/dt in K/s, Q, growth in kg/s)` of a core at `temperature_k` under a box at `mantle_k`."""
        heat_w = self._flux(temperature_k, mantle_k)
        growth_per_k = self._growth_per_k(temperature_k)
        effective = self.capacity + (LATENT_HEAT_J_KG + GRAVITATIONAL_ENERGY_J_KG) * growth_per_k
        rate = -heat_w / effective
        return rate, heat_w, -rate * growth_per_k

    def start(self, mantle_k):
        """Sets the core's temperature at formation: `INITIAL_EXCESS_K` above the base of the mantle."""
        self.temperature_k = self.boundary_k(mantle_k) + INITIAL_EXCESS_K

    def advance(self, age_gy, step_gy, mantle_start_k, mantle_end_k):
        """Cools the core for `step_gy` from `age_gy` while the box goes from `mantle_start_k` to `mantle_end_k` (midpoint rule)."""
        h = step_gy * GY_S
        k1, _, _ = self._cooling(self.temperature_k, mantle_start_k)
        middle_k = 0.5 * (mantle_start_k + mantle_end_k)
        k2, _, _ = self._cooling(self.temperature_k + 0.5 * h * k1, middle_k)
        self.temperature_k += h * k2
        power = self._power(self.temperature_k, mantle_end_k)
        if self._running and power <= 0.0:
            self._running = False
            self.ended_gy = age_gy + step_gy
        elif not self._running and power > 0.0:
            self._running = True
            self.ended_gy = None

    def _power(self, temperature_k, mantle_k):
        """The dynamo's power, W: negative or zero while it does not run."""
        if 1.0 - self.inner_radius_m(temperature_k) / self.core_m < MIN_LIQUID_SHELL:
            return 0.0
        _, heat_w, growth_kg_s = self._cooling(temperature_k, mantle_k)
        boundary_k = self.boundary_k(mantle_k)
        efficiency = 1.0 - boundary_k / temperature_k if temperature_k > boundary_k else 0.0
        return efficiency * (heat_w - self.adiabatic_w) + GRAVITATIONAL_ENERGY_J_KG * growth_kg_s

    def state(self, mantle_k):
        """The `Core` now, under a box at `mantle_k`."""
        power = self._power(self.temperature_k, mantle_k)
        _, heat_w, _ = self._cooling(self.temperature_k, mantle_k)
        return Core(self.temperature_k, self.boundary_k(mantle_k), heat_w, self.adiabatic_w,
                    self.inner_radius_m(self.temperature_k) / self.core_m, max(power, 0.0), None if power > 0.0 else self.ended_gy,
                    self.core_m, self.density)
