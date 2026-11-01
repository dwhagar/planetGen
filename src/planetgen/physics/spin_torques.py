# planetgen/physics/spin_torques.py

"""
The torques on a planet's spin: its star's tide, each moon's tide and the thermal tide of its air (GEN.956)
=========================================================================================================

A constant-Q despinning by the tide of what a body orbits (the form this module's gravitational torque reduces to for one raiser and a constant `Q`) is the
simplest case. A planet's spin has three more things acting on it: the tide each of its **moons** raises in it (the Moon slows the Earth and drives itself outward; a moon inside the planet's
corotation radius is dragged in), the **thermal tide** of a thick atmosphere, whose heated bulge is pulled by the star the other way and can
hold a planet away from synchronous rotation (Venus, the non-synchronous planets of [R2606]), and the planet's own **rheology**,
which makes the gravitational torque depend on frequency (`tidal_rheology`). This module is the sum of those torques and the evolution of a spin
under it. `body_spin` fills it from a body's stored fields and `planets.generate_spin` runs it (GEN.957).

**Gravitational torque.** A body of radius `R` spinning at `omega` under the tide of a mass `M` on an orbit of radius `a` and rate `n`, at the
semi-diurnal frequency `sigma = 2 (omega - n)`:

    N_g = - F K_g b_g(sigma),      K_g = G M^2 R^5 / a^6,      b_g(sigma) = sign(sigma) k2/Q (|sigma|)

`F` is `GRAVITATIONAL_TORQUE_FACTOR`: 3 as the despinning has had it since GEN.731 (af[R1584]; SP7) and 3/2 in Murray and Dermott and in
[R2606] (the same constant-Q result written for a different definition of Q). `k2/Q` comes from a callable, so it can be constant
(`constant_dissipation`, today's behaviour: a torque that does not depend on the spin, flipping sign at synchronous rotation) or the
frequency-dependent rheology of the body (`rheology_dissipation`, which has no discontinuity at synchronism and, at low frequency, a torque
proportional to it, so that asynchronous states can exist).

**Thermal torque** (Leconte, Wu, [R2606]:632, [R2606], equations 1 to 3, read in full):

    N_a = + F K_a q0 x / (1 + x^2),      K_a = 3 M R^3 / (5 rho a^3),      x = (omega - n) / omega0

with `q0` the amplitude of the degree-2 surface pressure field of the tide at zero frequency (Pa) and `omega0` the inverse of the time the
atmosphere and surface take to come to thermal equilibrium: a Debye (Maxwell-like) response, positive for a prograde spin, opposite to the
gravitational torque. The paper has no scaling law for either, only its Table 1 of 3D climate runs (Venus 92 bar, 2610 W/m2: 201 Pa and 3.77e-7 /s;
Earth-like air at 1366 W/m2: 1,180 Pa and 2.30e-6 /s at 1 bar, 4,050 Pa and 1.46e-6 /s at 10 bar; at 450 W/m2: 890 and 2,960 Pa, 1.18e-6 and
7.17e-7 /s): `thermal_pressure_amplitude_pa` and `thermal_relaxation_rate` interpolate that table (stand-ins, registry rows TQ1 to TQ3). [R0199] ([R0199], read in full) give the same Debye shape, `sigma/(sigma^2 + sigma0^2)` for a convective atmosphere, and
find that a stably stratified one flattens the torque to nothing: the sum here is the convective (maximum) case, and a stratified air is the
caller's to switch off.

**Equilibria.** The spin obeys `I d omega/dt = sum of the torques`. It is a one-dimensional flow, so it runs from its birth value in the
direction of the net torque until it meets an attractor: a synchronous state (torque changing sign at the orbital rate), a spin-orbit
resonance whose window of eccentricity holds the orbit's (`tidal_spin.RESONANCES`, [R2823]) that takes the spin with the chance `tidal_spin.RESONANCE_CAPTURE_PROBABILITY` (a uniform the body carries, GEN.1172), or a zero of the smooth torque where the
thermal and gravitational tides balance (an asynchronous state, prograde or retrograde: `equilibria`). With constant `Q` those zeros are at
`x = (1 +- sqrt(1 - 4 c^2)) / (2 c)`, `c = K_g (k2/Q) / (K_a q0)`, which exist when `c <= 1/2`, the `gamma_Q >= 1` of Leconte et al.'s equation 12 (their
equation 13), the outer one stable.

**What this does not do.** The planet's obliquity is taken as zero for the torques; one orbit of each raiser, circular for the torque (the eccentricity only
selects resonances); no core-mantle friction (the second half of Venus's story, [R0895]); the thermal torque is that of a
convective atmosphere with a fitted amplitude; moons act on the spin only, and `moon_orbits` (GEN.958) moves the moons. Stand-ins are in
`docs/design/research/gap-filling-plan.md` section 6.24.
"""

import math
import warnings
from dataclasses import dataclass

from scipy import integrate, optimize

from planetgen.physics import constants, tidal_rheology, tidal_spin

GRAVITATIONAL_TORQUE_FACTOR = 3.0
"""float: `F` of the gravitational torque as the despinning has had it since GEN.731 [R1584, R2606] and Murray and
Dermott write 3/2 for the same constant-Q result. The thermal torque is scaled with it, so the ratio of the two stays the paper's. Source tier T3,
tweak level L0."""
THERMAL_TABLE = (
    # (surface pressure bar, q0 Pa at 1366 W/m2, 1 / omega0 s at 1366 W/m2)
    (1.0, 1180.0, 1.0 / 2.30e-6),
    (10.0, 4050.0, 1.0 / 1.46e-6),
    (92.0, 201.0 / 1.19, 1.62 / 3.77e-7),
)
"""tuple: [R2606: Table 1] brought to the flux of 1366 W/m2 by `FLUX_EXPONENT_Q0` and `FLUX_EXPONENT_OMEGA0`: the pressure amplitude and
the thermal time of Earth-like air at 1 and 10 bar and of Venus's 92 bar (2,610 W/m2). Source tier T2 (read in the paper), tweak level L1 (the
Venus row is moved to the common flux with the fitted exponent)."""
FLUX_EXPONENT_Q0 = 0.27
"""float: `q0 ~ F^0.27`: the table's 1,366 against 450 W/m2 at 1 and 10 bar give 0.25 and 0.28 (the paper says only that q0 changes by under 25 percent
across the habitable zone). Source tier T6, tweak level L2 (fitted to the table)."""
FLUX_EXPONENT_OMEGA0 = 0.75
"""float: `omega0 ~ F^(3/4)` (the paper: `omega0 ~ T0^3 ~ F^(3/4)`). Source tier T1, tweak level L0."""
THERMAL_REFERENCE_FLUX_W_M2 = 1366.0
"""float: The flux of the common column of `THERMAL_TABLE`."""
THERMAL_SMALL_PRESSURE_SLOPE = 1.0
"""float: Below the table's 1 bar, `q0` is proportional to the surface pressure (the paper's argument that the amplitude goes as the pressure over the scale height). Source
tier T5, tweak level L0."""
EQUILIBRIUM_TOLERANCE = 1.0e-3
"""float: How close, in units of the distance from an equilibrium to the orbital rate, the spin must come to count as having reached it."""
SCAN_POINTS = 240
"""int: Grid points per side of the orbital rate where `equilibria` looks for a sign change of the smooth torque."""
STOP_CLOSENESS = 1.0e-9
"""float: The flow to a zero of the torque ends within this share of the distance travelled (the approach is exponential, so the log of it is all the time it adds). Source tier T6, tweak level L0."""
TORQUE_FLOOR_N_M = 1.0e-30
"""float: The least torque the flow integral divides by, N m; below it the spin does not move in any age. Source tier T6, tweak level L0."""


@dataclass(frozen=True)
class Planet:
    """
    The spinning body.

    Attributes:
        mass_kg, radius_m (float): Its mass and radius.
        moment_factor (float): `alpha = I / (m R^2)` (`tidal_spin.ROCKY_MOMENT_FACTOR`).
    """

    mass_kg: float
    radius_m: float
    moment_factor: float = tidal_spin.ROCKY_MOMENT_FACTOR

    @property
    def inertia(self):
        """`I`, kg m2."""
        return self.moment_factor * self.mass_kg * self.radius_m ** 2

    @property
    def mean_density(self):
        """kg/m3."""
        return self.mass_kg / (4.0 / 3.0 * math.pi * self.radius_m ** 3)


@dataclass(frozen=True)
class Raiser:
    """
    A mass whose gravity raises a tide in the planet.

    Attributes:
        mass_kg (float): Its mass (the star's, or a moon's).
        distance_m (float): Its distance from the planet's centre of attraction (the star's: the planet's semi-major axis; a moon's: its own).
        orbital_rate (float): The rate at which it goes round, rad/s, signed (negative for a retrograde moon).
    """

    mass_kg: float
    distance_m: float
    orbital_rate: float


def orbital_rate(host_mass_kg, distance_m):
    """The mean motion of a body at `distance_m` round `host_mass_kg`, rad/s."""
    return math.sqrt(constants.G * host_mass_kg / distance_m ** 3)


def constant_dissipation(k2, q):
    """The `k2/Q` of a body that does not depend on frequency (today's rheology), as the callable `Torques` takes."""
    value = k2 / q

    def dissipation(frequency):
        return value
    dissipation.constant = value
    return dissipation


def rheology_dissipation(planet, temperature_k, eta0_pa_s, activation_j_mol, solidus_k, liquidus_k):
    """The frequency-dependent `k2/Q` of a rocky planet at a mantle temperature (`tidal_rheology.rock_response`), as the callable `Torques` takes."""
    melt = tidal_rheology.melt_fraction(temperature_k, solidus_k, liquidus_k)

    def dissipation(frequency):
        return tidal_rheology.response(planet.mass_kg, planet.radius_m, max(frequency, 1.0e-20), temperature_k, eta0_pa_s, activation_j_mol,
                                       melt).k2_over_q
    return dissipation


TABLE_RANGE_RAD_S = (1.0e-14, 1.0e-1)
"""tuple: The tide frequencies `tabulate` covers (a tide a million years long to one of a minute). Source tier T6, tweak level L0."""
TABLE_POINTS = 56
"""int: Points of `tabulate`'s table, log-spaced; a dissipation that varies by a power of the frequency between them is good to a part in a thousand. Source tier T6, tweak level L0."""


def tabulate(dissipation, frequency_range=TABLE_RANGE_RAD_S, points=TABLE_POINTS):
    """
    `dissipation` (`k2/Q` against tide frequency) as a table in log-log, for one that is costly to evaluate (an ice shell integrates over its temperatures at each
    frequency, `icy_tide.k2_over_q`); the flow of a spin evaluates it thousands of times. Outside the range the end values are continued with the end slopes. A zero
    anywhere (a body with no shell) gives the constant zero.
    """
    low, high = math.log(frequency_range[0]), math.log(frequency_range[1])
    logs = [low + (high - low) * i / (points - 1) for i in range(points)]
    values = [dissipation(math.exp(x)) for x in logs]
    if not all(v > 0.0 for v in values):
        if not any(v > 0.0 for v in values):
            return lambda frequency: 0.0
        values = [max(v, 1.0e-300) for v in values]
    ys = [math.log(v) for v in values]
    slope_low = (ys[1] - ys[0]) / (logs[1] - logs[0])
    slope_high = (ys[-1] - ys[-2]) / (logs[-1] - logs[-2])
    step = (high - low) / (points - 1)

    def table(frequency):
        x = math.log(max(frequency, 1.0e-300))
        if x <= low:
            return math.exp(ys[0] + slope_low * (x - low))
        if x >= high:
            return math.exp(ys[-1] + slope_high * (x - high))
        i = min(int((x - low) / step), points - 2)
        return math.exp(ys[i] + (ys[i + 1] - ys[i]) * (x - logs[i]) / step)
    return table


def thermal_pressure_amplitude_pa(surface_pressure_pa, flux_w_m2):
    """`q0`, Pa: the zero-frequency amplitude of the thermal tide's pressure field, from `THERMAL_TABLE` (log-log in pressure above 1 bar, linear
    in pressure below it) and `FLUX_EXPONENT_Q0`."""
    return _thermal_table(surface_pressure_pa, 1)[0] * (flux_w_m2 / THERMAL_REFERENCE_FLUX_W_M2) ** FLUX_EXPONENT_Q0


def thermal_relaxation_rate(surface_pressure_pa, flux_w_m2):
    """`omega0`, rad/s: the inverse of the atmosphere's thermal time, from `THERMAL_TABLE` (the thermal time is linear in the column mass) and
    `FLUX_EXPONENT_OMEGA0`."""
    return 1.0 / _thermal_table(surface_pressure_pa, 2)[0] * (flux_w_m2 / THERMAL_REFERENCE_FLUX_W_M2) ** FLUX_EXPONENT_OMEGA0


def _thermal_table(surface_pressure_pa, column):
    bar = surface_pressure_pa / 1.0e5
    rows = THERMAL_TABLE
    if bar <= rows[0][0]:
        if column == 1:
            return (rows[0][1] * (bar / rows[0][0]) ** THERMAL_SMALL_PRESSURE_SLOPE,)
        # the thermal time: the surface's own heat capacity is the part of the 1 bar value that does not go with the column
        return (rows[0][2] * (0.94 + 0.06 * bar / rows[0][0]),)
    for low, high in zip(rows, rows[1:]):
        if bar <= high[0] or high is rows[-1]:
            share = math.log(bar / low[0]) / math.log(high[0] / low[0])
            return (math.exp(math.log(low[column]) + share * (math.log(high[column]) - math.log(low[column]))),)
    raise AssertionError("unreachable")


@dataclass(frozen=True)
class ThermalTide:
    """
    The thermal tide of a thick air under a star.

    Attributes:
        star: The `Raiser` that heats it.
        amplitude_pa (float): `q0`.
        relaxation_rate (float): `omega0`, rad/s.
    """

    star: Raiser
    amplitude_pa: float
    relaxation_rate: float


@dataclass(frozen=True)
class Torques:
    """
    Everything that acts on one planet's spin.

    Attributes:
        planet (Planet): The planet.
        raisers (tuple): Its gravitational `Raiser`s: the star first (spin-orbit resonances are measured against the first one's orbital rate), then moons.
        dissipation: `k2/Q` of the planet as a function of the tide's angular frequency (`constant_dissipation`, `rheology_dissipation`).
        thermal (ThermalTide, optional): The air's thermal tide.
        factor (float): `GRAVITATIONAL_TORQUE_FACTOR`.
        softness (float): 0 for the sharp sign flip of a constant-Q torque at each raiser's orbital rate; above 0 the sign is replaced by
            `tanh(sigma / (2 softness |n|))` (the regularisation of [R0198] equation 77, which uses 1e-2 of the orbital rate), for
            integrations that cannot step across a jump (`moon_orbits`).
    """

    planet: Planet
    raisers: tuple
    dissipation: object
    thermal: ThermalTide = None
    factor: float = GRAVITATIONAL_TORQUE_FACTOR
    softness: float = 0.0

    def gravitational(self, spin, raiser):
        """The torque on the planet of one raiser's tide, N m (negative slows a prograde spin)."""
        sigma = 2.0 * (spin - raiser.orbital_rate)
        if sigma == 0.0:
            return 0.0
        k_g = constants.G * raiser.mass_kg ** 2 * self.planet.radius_m ** 5 / raiser.distance_m ** 6
        sign = math.copysign(1.0, sigma)
        if self.softness > 0.0:
            sign = math.tanh(sigma / (2.0 * self.softness * max(abs(raiser.orbital_rate), 1.0e-30)))
        return -sign * self.factor * k_g * self.dissipation(abs(sigma))

    def atmospheric(self, spin):
        """The thermal tide's torque, N m (positive for a prograde spin)."""
        if self.thermal is None:
            return 0.0
        star = self.thermal.star
        k_a = 3.0 * star.mass_kg * self.planet.radius_m ** 3 / (5.0 * self.planet.mean_density * star.distance_m ** 3)
        x = (spin - star.orbital_rate) / self.thermal.relaxation_rate
        return self.factor * k_a * self.thermal.amplitude_pa * x / (1.0 + x * x)

    def moon_torque(self, spin, index):
        """The torque of raiser `index` (0 is the star)."""
        return self.gravitational(spin, self.raisers[index])

    def total(self, spin):
        """The sum of all the torques at `spin`, N m."""
        return math.fsum([self.gravitational(spin, raiser) for raiser in self.raisers] + [self.atmospheric(spin)])

    def smooth(self, spin):
        """The torque with each raiser's sign flip at its own orbital rate removed (the part that has zeros): the total evaluated on the side of
        the flip the spin is on."""
        return self.total(spin)


def equilibria(torques, eccentricity=0.0):
    """
    The spin rates, rad/s, at which the net torque on the planet is zero and the spin stays: `(rate, kind)` sorted by rate, `kind` one of
    `"synchronous"` (a raiser's sign flip with the torque pushing toward it from both sides), `"asynchronous"` (a smooth zero with the torque falling through it) and
    `"unstable"` (a smooth zero with the torque rising, listed so tests can see the pair). Resonances are not listed: they are captured in passing.
    """
    found = []
    scale = max([abs(r.orbital_rate) for r in torques.raisers] + [torques.thermal.relaxation_rate if torques.thermal else 0.0])
    centres = sorted({r.orbital_rate for r in torques.raisers})
    for centre in centres:
        eps = 1.0e-9 * max(abs(centre), 1.0e-12)
        below, above = torques.total(centre - eps), torques.total(centre + eps)
        if below > 0.0 > above:
            found.append((centre, "synchronous"))
    span = 40.0 * scale
    grid = [centres[0] - span + 2.0 * span * i / (2 * SCAN_POINTS) for i in range(2 * SCAN_POINTS + 1)]
    grid = sorted(set(grid) | {c + s * span * 10.0 ** (-k / 20.0) for c in centres for s in (-1.0, 1.0) for k in range(0, 100)})
    values = [torques.total(w) for w in grid]
    flips = {c for c in centres}
    for (w0, v0), (w1, v1) in zip(zip(grid, values), zip(grid[1:], values[1:])):
        if v0 == 0.0 or v0 * v1 > 0.0:
            continue
        if any(w0 < c < w1 for c in flips) or any(math.isclose(c, w0) or math.isclose(c, w1) for c in flips):
            continue
        root = optimize.brentq(torques.total, w0, w1, xtol=1.0e-18 + 1.0e-12 * abs(w0), rtol=1.0e-12)
        found.append((root, "asynchronous" if v0 > 0.0 > v1 else "unstable"))
    return sorted(found)


@dataclass(frozen=True)
class SpinResult:
    """
    Where a spin ends.

    Attributes:
        rate (float): The spin after the age, rad/s.
        ratio (float): In orbital rates of the first raiser, if it has been caught by a synchronous state or a resonance (1, 1.5, 2), else `None`.
        kind (str): `"synchronous"`, `"resonance"`, `"corotation"` (caught at a moon's orbital rate), `"asynchronous"` (at a balance of the thermal and
            gravitational tides) or `"free"` (still on its way).
    """

    rate: float
    ratio: float
    kind: str


CAPTURE_RATIOS = tuple(sorted(tidal_spin.RESONANCE_CAPTURE_PROBABILITY, reverse=True))
"""tuple: The resonances a spin can pass without being caught (2.0, 1.5), from the highest ratio down; `capture_units` of `evolve` is one uniform for each, in this order."""
ALWAYS_CAPTURED = (0.0,) * len(CAPTURE_RATIOS)
"""tuple: `capture_units` that take every resonance whose window holds the eccentricity (the lowest uniform is under any chance above zero), for a caller that wants the window alone."""


def _window_ratios(eccentricity):
    return [target for target, low, high in tidal_spin.RESONANCES if low < eccentricity < high and target in tidal_spin.RESONANCE_CAPTURE_PROBABILITY]


def evolve(torques, spin_rate, age_s, eccentricity=0.0, capture_units=ALWAYS_CAPTURED):
    """
    A planet born spinning at `spin_rate` (rad/s) after `age_s` under `torques`: the one-dimensional flow `I d omega/dt = N(omega)`, run in
    the direction of the net torque until it is caught (a resonance whose window holds `eccentricity` and whose chance the body's uniform is under, a synchronous state, an
    asynchronous balance) or the age is up. A resonance whose window holds `eccentricity` but that does not catch the spin lets it go on down. With constant `Q` and one
    raiser the spin falls linearly at the rate `F K k2 / (I Q)` to the first stop.

    Args:
        capture_units (tuple): One uniform in [0, 1) for each of `CAPTURE_RATIOS`; the spin is caught at a resonance when its uniform is below
            `tidal_spin.RESONANCE_CAPTURE_PROBABILITY` of it. A body keeps its own (`planets.draw_birth_spin`), so running its history again gives the same state.
    """
    inertia = torques.planet.inertia
    star_rate = torques.raisers[0].orbital_rate
    # stops along the way: resonances of the star's orbital rate, each raiser's sign flip, and the smooth zeros
    stops = {}
    for ratio in _window_ratios(eccentricity):
        stops[ratio * star_rate] = ("resonance", ratio)
    for raiser in torques.raisers:
        stops.setdefault(raiser.orbital_rate, ("flip", 1.0 if raiser is torques.raisers[0] else None))
    for rate, kind in equilibria(torques, eccentricity):
        if kind == "asynchronous":
            stops.setdefault(rate, ("balance", None))
    spin, remaining = spin_rate, age_s
    for _ in range(4 * len(stops) + 8):
        push = torques.total(spin)
        if push == 0.0:
            return SpinResult(spin, None, "free")
        direction = 1.0 if push > 0.0 else -1.0
        ahead = sorted((w for w in stops if (w - spin) * direction > 1.0e-12 * max(abs(spin), 1.0e-12)), key=lambda w: abs(w - spin))
        if not ahead:
            # nothing in the way: the torque runs to the end of the age
            return SpinResult(_advance(torques, inertia, spin, math.inf * direction, remaining), None, "free")
        stop = ahead[0]
        needed = _flow_time(torques, inertia, spin, stop, direction)
        if needed > remaining:
            return SpinResult(_advance(torques, inertia, spin, stop, remaining), None, "free")
        remaining -= needed
        spin = stop
        kind, ratio = stops[stop]
        if kind == "resonance":
            if capture_units[CAPTURE_RATIOS.index(ratio)] < tidal_spin.RESONANCE_CAPTURE_PROBABILITY[ratio]:
                return SpinResult(stop, ratio, "resonance")
            spin = stop + direction * 1.0e-9 * max(abs(stop), 1.0e-12)      # not caught: on down to the next stop
            continue
        if kind == "balance":
            return SpinResult(stop, None, "asynchronous")
        eps = 1.0e-9 * max(abs(stop), 1.0e-12)
        beyond = torques.total(stop + direction * eps)
        if beyond * direction < 0.0:
            return SpinResult(stop, ratio, "synchronous" if ratio is not None else "corotation")
        spin = stop + direction * eps
    return SpinResult(spin, None, "free")


def _flow_time(torques, inertia, start, stop, direction):
    """The time for the spin to go from `start` to `stop`, s: `integral I d omega / |N|` along the way."""
    lo, hi = sorted((start, stop))
    if _is_constant(torques, lo, hi):
        return inertia * (hi - lo) / abs(torques.total(0.5 * (lo + hi)))
    # A balance (or a synchronous state under a torque that fades with the tide's frequency) is a zero of the torque at `stop`: the spin comes to it
    # exponentially and the integral diverges as the log of the distance, so the flow ends within STOP_CLOSENESS of the span of it.
    margin = STOP_CLOSENESS * (hi - lo)
    if stop >= start:
        hi -= margin
    else:
        lo += margin
    with warnings.catch_warnings():
        # the integrand climbs steeply into the zero; the value is still good to the tolerance the ages need
        warnings.simplefilter("ignore", integrate.IntegrationWarning)
        value, _ = integrate.quad(lambda w: inertia / max(abs(torques.total(w)), TORQUE_FLOOR_N_M), lo, hi, epsabs=0.0, epsrel=1.0e-6, limit=200)
    return value


def _advance(torques, inertia, start, stop, available_s):
    """The spin after `available_s` of flow from `start` toward `stop` (`stop` may be infinite: the torque is then constant by assumption)."""
    direction = 1.0 if stop > start else -1.0
    if math.isinf(stop) or _is_constant(torques, *sorted((start, stop))):
        return start + direction * abs(torques.total(start)) * available_s / inertia
    total = _flow_time(torques, inertia, start, stop, direction)
    if available_s >= total:
        return stop
    return _solve_flow(torques, inertia, start, stop, direction, available_s)


def _solve_flow(torques, inertia, start, stop, direction, available_s):
    """Bisection on the end of the flow integral (the integral is monotone in the end point)."""
    low, high = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (low + high)
        w = start + (stop - start) * mid
        if _flow_time(torques, inertia, start, w, direction) < available_s:
            low = mid
        else:
            high = mid
    return start + (stop - start) * 0.5 * (low + high)


def _is_constant(torques, lo, hi):
    """Whether the net torque is the same at three points of `[lo, hi]` (constant-Q gravitational torques are, between sign flips)."""
    if torques.thermal is not None or getattr(torques.dissipation, "constant", None) is None:
        return False
    values = [torques.total(lo + (hi - lo) * s) for s in (0.01, 0.5, 0.99)]
    return all(abs(v - values[1]) <= 1.0e-12 * abs(values[1]) for v in values)
