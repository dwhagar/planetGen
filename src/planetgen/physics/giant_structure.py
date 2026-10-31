# planetgen/physics/giant_structure.py

"""
The hydrostatic structure of a giant planet, an ice giant or a brown dwarf from the hydrogen-helium equation of state
(GEN.733, docs/design/composition/giant-interiors.md section 11, registry 6.134)
==============================================================================================================

The pressure at every mass inside a giant, and the radius there, from the equations of a self-gravitating fluid in
hydrostatic balance,

    dP/dm = -G m / (4 pi r^4),      dr/dm = 1 / (4 pi r^2 rho),

with the density of the hydrogen-helium envelope read from the downloaded tables of Chabrier and Debras
(`datasets/hhe_eos.py`) along an adiabat, a line of the tables' entropy (the temperature at each pressure is where that entropy is
met, `Grid.adiabat_log_t`; the tables' adiabatic gradient carries numerical noise). It
replaces the n = 1 polytrope the layer pressures were read from (`giant_interior.Polytrope`): a polytrope has one index for the
whole body, where the real fluid goes from an ideal molecular gas to a degenerate metal.

The body is built from its centre in up to three layers (a `Layer`): a core of the heavy elements the envelope does not hold, an
ice mantle (an ice giant), and the envelope of hydrogen and helium with the envelope's share of heavy elements mixed in by the
ideal volume law. The density of the heavy elements is the repo's own: iron's Vinet fit, the silicate mantle table and the
water isotherm of `bulk_structure`, tabulated on a grid of pressure once. The integration goes inward from the one-bar level in
the logarithm of the pressure (fourth-order Runge-Kutta) and stops where the mass inside falls to `CENTRE_MASS_FRACTION` of the
body's, closing the rest as a sphere of the density it has there.

What fixes the entropy. The radius and the one-bar temperature are both inputs of the generator (the radius from the tables of
GEN.460, the temperature from the starlight), and a model that took both would not close: an irradiated or inflated giant is not on
the adiabat its one-bar temperature starts. So the radius is kept (the caller fixed it) and the entropy of the adiabat is the
unknown, found by `fit_adiabat` so that the model closes. Its one-bar temperature is the model's own and is reported
(`Profile.adiabat_1bar_k`); it is the generator's one-bar temperature for a body on its own adiabat and differs for one that is not
(GS5). A radius no adiabat inside `ADIABAT_RANGE_K` closes (a body too light for its radius on the coldest, or too dense on the hottest)
takes that adiabat and the radius it does close at, and says so (`Profile.closed` False, `Profile.radius_m` the model's own).

Known limits (registry 6.134): the heavy elements are in a core and a uniformly mixed envelope, with no gradient or dilute core
(GS2); the helium share is one value for the whole body, the natal one (GS3); the heavy-element density is a volume sum of the
repo's iron, silicate and water curves at their fixed temperatures, and ammonia and methane ice are given water's (GS4); the
adiabat is a line of the mixture table's entropy, and the table starts at 100 K, so the first scale heights of a body colder than
that are read at 100 K (GS5); and the table's solid-phase and low-density values are unphysical, as its README says (GS6).

Use:

    subject = StructureSubject(mass_kg=1.898e27, radius_m=6.9911e7, helium_share=0.27, envelope_z=0.04, heavy_mass_kg=1.9e26)
    profile = fit_adiabat(subject)
    profile.at_mass_inside(0.9)           # (radius fraction, pressure in Pa)
    profile.at_pressure(2.0e11)           # (radius fraction, mass inside fraction)
"""

import math
from bisect import bisect_left
from dataclasses import dataclass, field

from planetgen.datasets import hhe_eos
from planetgen.physics import bulk_structure, constants

G = constants.G
ONE_BAR_PA = 1.0e5
"""float: The pressure the integration starts from (Pa), the 1 bar level the radii are given at."""

STEP_LN_PA = 0.25
"""float: The step of the integration in ln P (about 70 steps from one bar to a few terapascals). A smaller step moves the radius
a model closes at by under 0.3 percent (tested)."""

ADIABAT_RANGE_K = (30.0, 8000.0)
"""tuple: The one-bar temperatures of the adiabats the fit may use (K). The low end is under the table's 100 K edge, where the
adiabat is carried up the molecular gas's own (`hhe_eos.Grid.log_entropy_of_surface_temperature`, GS5)."""

RADIUS_TOLERANCE = 2.0e-4
"""float: How close to closing (the defect, a share of the radius) the fitted adiabat is taken."""

MAX_EVALUATIONS = 40
"""int: The most integrations one fit makes."""

STEP_FRACTION = 0.1
"""float: The most the radius or the mass inside may change in one step as a share of itself (a dense core changes them fast)."""

CENTRE_MASS_FRACTION = 2.0e-3
"""float: The share of the body's mass left inside the radius where the integration stops; the rest is closed as a uniform sphere."""

CENTRE_RADIUS_FRACTION = 0.02
"""float: The integration also stops at this share of the radius (a model too light for its radius reaches it with mass to spare)."""

CENTRAL_PRESSURE_LIMIT_PA = 1.0e22
"""float: The pressure past which the integration stops (the table's end, 1e13 GPa)."""

HEAVY_NODES_PA = (1.0e5, 1.0e17, 220)
"""tuple: The first and last pressure (Pa) and the number of nodes of the tabulated heavy-element densities."""

LAYER_KINDS = ("core", "mantle", "envelope")
"""tuple: The layers from the centre out, as a body is built."""


class StructureError(ValueError):
    """A subject that cannot describe a body."""


# ---------------------------------------------------------------------------------------------------------------------
# The heavy elements
# ---------------------------------------------------------------------------------------------------------------------

_P_FIRST, _P_LAST, _N_NODES = HEAVY_NODES_PA
_U0 = math.log(_P_FIRST)
_DU = (math.log(_P_LAST) - _U0) / (_N_NODES - 1)
_CURVES = {}


def _curve(name):
    """The specific volume (m3/kg) of one pure heavy material at the nodes of ln P: `iron`, `rock` or `water`."""
    if name not in _CURVES:
        function = {"iron": bulk_structure.core_density, "rock": bulk_structure.mantle_density,
                    "water": bulk_structure.water_density}[name]
        _CURVES[name] = tuple(1.0 / function(math.exp(_U0 + i * _DU)) for i in range(_N_NODES))
    return _CURVES[name]


def heavy_volume_table(ice_fraction, iron_share):
    """
    The specific volume (m3/kg) of a mix of the heavy elements at the nodes of ln P, by the ideal volume law: `ice_fraction` of
    the mass is water-like ice (ammonia and methane ice take water's curve, GS4), the rest refractory, of which `iron_share` is iron
    and the rest silicate.
    """
    if not (0.0 <= ice_fraction <= 1.0 and 0.0 <= iron_share <= 1.0):
        raise StructureError("the ice fraction and the iron share of the heavy elements are fractions")
    iron, rock, water = _curve("iron"), _curve("rock"), _curve("water")
    refractory = 1.0 - ice_fraction
    return tuple(ice_fraction * water[i] + refractory * (iron_share * iron[i] + (1.0 - iron_share) * rock[i]) for i in range(_N_NODES))


# ---------------------------------------------------------------------------------------------------------------------
# The subject and the profile
# ---------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Layer:
    """
    One layer of the model, from the centre out.

    Attributes:
        kind (str): `core` (heavy elements only, the ice share `StructureSubject.core_ice_fraction`), `mantle` (ice only) or
            `envelope` (hydrogen and helium with `StructureSubject.envelope_z` of heavy elements mixed in).
        mass_kg (float): The layer's mass; the outermost layer takes the rest of the body's mass, so its own is not read.
    """

    kind: str
    mass_kg: float = 0.0

    def __post_init__(self):
        if self.kind not in LAYER_KINDS or self.mass_kg < 0.0:
            raise StructureError(f"a layer is one of {LAYER_KINDS} with a mass of zero or more, not {self.kind!r} {self.mass_kg}")


@dataclass(frozen=True)
class StructureSubject:
    """
    What the structure is built from.

    Attributes:
        mass_kg, radius_m (float): The body's, the radius at the one-bar level.
        helium_share (float): Y / (X + Y) of the hydrogen and helium.
        envelope_z (float): The heavy elements' share of the envelope's mass.
        layers (tuple): The `Layer`s from the centre out, the last the envelope. Empty: a core of the heavy elements the envelope
            does not hold (`heavy_mass_kg` less what `envelope_z` of the body takes) under the envelope; with none left over the
            envelope is the whole body and holds all of `heavy_mass_kg`.
        heavy_mass_kg (float): The mass of all the heavy elements, used when `layers` is empty.
        envelope_ice_fraction, core_ice_fraction (float): The ice share of the heavy elements in the envelope and in the core.
        iron_share (float): The iron share of the refractory heavy elements.
        temperature_1bar_k (float): The one-bar temperature the fit starts its search from (the generator's, from the starlight).
    """

    mass_kg: float
    radius_m: float
    helium_share: float
    envelope_z: float = 0.0
    layers: tuple = ()
    heavy_mass_kg: float = 0.0
    envelope_ice_fraction: float = 0.0
    core_ice_fraction: float = 0.0
    iron_share: float = 0.3
    temperature_1bar_k: float = 150.0

    def __post_init__(self):
        if not (self.mass_kg > 0.0 and self.radius_m > 0.0):
            raise StructureError("a giant needs a positive mass and radius")
        if not (0.0 < self.helium_share < 1.0 and 0.0 <= self.envelope_z < 1.0):
            raise StructureError("the helium share and the envelope's heavy-element share are fractions")
        if not (0.0 <= self.heavy_mass_kg < self.mass_kg):
            raise StructureError("the heavy-element mass must be under the body's mass")
        if self.layers and self.layers[-1].kind != "envelope":
            raise StructureError("the outermost layer is the envelope")

    def resolved(self):
        """`(layers, envelope_z)`: the layers from the centre out and the envelope's heavy-element share."""
        if self.layers:
            return self.layers, self.envelope_z
        core = (self.heavy_mass_kg - self.envelope_z * self.mass_kg) / (1.0 - self.envelope_z)
        if core > 0.0:
            return (Layer("core", core), Layer("envelope")), self.envelope_z
        return (Layer("envelope"),), self.heavy_mass_kg / self.mass_kg


@dataclass(frozen=True)
class Profile:
    """
    The structure: pressure, radius and mass inside along the model, outermost first. Read through `at_pressure` and
    `at_mass_inside`, the interface of `giant_interior.Polytrope`.

    Attributes:
        mass_kg, radius_m (float): The body's mass and the radius the model has: the subject's, or (not `closed`) its own.
        adiabat_1bar_k (float): The one-bar temperature of the adiabat the model closes on (K); under 100 K it is the table's 100 K edge
            carried to one bar along the molecular gas's adiabat (GS5), good to a few percent.
        closed (bool): Whether the model closed on the subject's radius. False: no adiabat inside `ADIABAT_RANGE_K` closes it, the
            nearest was taken and `radius_m` is the radius it closes at.
        central_pressure_pa (float): The pressure at the centre.
        central_temperature_k (float): The temperature the adiabat reaches at the pressure where the integration stopped (K); a core is not on it.
        evaluations (int): How many integrations the fit made.
    """

    mass_kg: float
    radius_m: float
    adiabat_1bar_k: float
    closed: bool
    central_pressure_pa: float
    central_temperature_k: float
    evaluations: int = 0
    _ln_p: tuple = field(repr=False, default=())
    _mass: tuple = field(repr=False, default=())
    _radius: tuple = field(repr=False, default=())

    def at_pressure(self, pressure_pa):
        """`(radius fraction, mass inside fraction)` where the pressure is `pressure_pa`; (1, 1) at or below the one-bar level and
        (0, 0) at or above the centre."""
        if pressure_pa >= self.central_pressure_pa:
            return 0.0, 0.0
        if pressure_pa <= ONE_BAR_PA:
            return 1.0, 1.0
        u = math.log(pressure_pa)
        i = min(max(bisect_left(self._ln_p, u) - 1, 0), len(self._ln_p) - 2)
        t = (u - self._ln_p[i]) / (self._ln_p[i + 1] - self._ln_p[i])
        return (self._radius[i] + t * (self._radius[i + 1] - self._radius[i]),
                self._mass[i] + t * (self._mass[i + 1] - self._mass[i]))

    def at_mass_inside(self, fraction):
        """`(radius fraction, pressure in Pa)` where the mass inside is `fraction` of the whole; (1, 0) for all of it and the centre for none."""
        if fraction >= 1.0:
            return 1.0, 0.0
        if fraction <= 0.0:
            return 0.0, self.central_pressure_pa
        masses = self._mass       # falling with depth
        lo, hi = 0, len(masses) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if masses[mid] > fraction:
                lo = mid
            else:
                hi = mid
        span = masses[lo] - masses[hi]
        t = 0.0 if span == 0.0 else (masses[lo] - fraction) / span
        u = self._ln_p[lo] + t * (self._ln_p[hi] - self._ln_p[lo])
        return self._radius[lo] + t * (self._radius[hi] - self._radius[lo]), math.exp(u)

    def pressure_at_mass_inside_pa(self, fraction):
        return self.at_mass_inside(fraction)[1]

    def at_radius(self, fraction):
        """`(mass inside fraction, pressure in Pa)` at radius fraction `fraction` of the radius the model has."""
        if fraction >= 1.0:
            return 1.0, ONE_BAR_PA
        if fraction <= 0.0:
            return 0.0, self.central_pressure_pa
        radii = self._radius      # falling with depth
        lo, hi = 0, len(radii) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if radii[mid] > fraction:
                lo = mid
            else:
                hi = mid
        span = radii[lo] - radii[hi]
        t = 0.0 if span == 0.0 else (radii[lo] - fraction) / span
        return (self._mass[lo] + t * (self._mass[hi] - self._mass[lo]),
                math.exp(self._ln_p[lo] + t * (self._ln_p[hi] - self._ln_p[lo])))


# ---------------------------------------------------------------------------------------------------------------------
# The integration
# ---------------------------------------------------------------------------------------------------------------------

class _Model:
    """The densities of one subject, ready to be integrated for any one-bar temperature."""

    def __init__(self, subject):
        self.subject = subject
        self.grid = hhe_eos.grid()
        layers, self.z = subject.resolved()
        self.kinds = [layer.kind for layer in layers]
        edges, total = [], 0.0
        for layer in layers[:-1]:
            total += layer.mass_kg
            edges.append(total)
        if total >= subject.mass_kg:
            raise StructureError("the layers hold more than the body's mass")
        self.edges = edges                                  # mass inside the top of each layer but the outermost, kg
        tables = {"core": heavy_volume_table(subject.core_ice_fraction, subject.iron_share),
                  "mantle": heavy_volume_table(1.0, subject.iron_share),
                  "envelope": heavy_volume_table(subject.envelope_ice_fraction, subject.iron_share)}
        self.volumes = [tables[kind] for kind in self.kinds]
        self.y_shift = subject.helium_share - hhe_eos.TABLE_HELIUM_SHARE
        self.mass_stop = CENTRE_MASS_FRACTION * subject.mass_kg

    def _volume_m3_per_kg(self, u, lt_nodes, layer):
        """The specific volume (m3/kg) of layer `layer` at ln P `u`, the envelope's on the adiabat whose temperatures at the table's
        pressure nodes are `lt_nodes` (log10 K)."""
        x = min(max((u - _U0) / _DU, 0.0), _N_NODES - 1.000001)
        i = int(x)
        table = self.volumes[layer]
        volume_z = table[i] + (x - i) * (table[i + 1] - table[i])
        if self.kinds[layer] != "envelope":
            return volume_z
        grid = self.grid
        log_p = u / _LN10 - 9.0
        y = (log_p - grid.log_p0) / grid.d_log_p
        y = min(max(y, 0.0), grid.n_p - 1.000001)
        j = int(y)
        log_t = lt_nodes[j] + (y - j) * (lt_nodes[j + 1] - lt_nodes[j])
        log_rho, slope = grid.sample(log_t, log_p)
        volume_hhe = 10.0 ** -log_rho * (1.0 + self.y_shift * slope) * 1.0e-3
        return (1.0 - self.z) * volume_hhe + self.z * volume_z

    def integrate(self, log_s0, record=None, radius_m=None):
        """
        Integrates from the one-bar level inward along the adiabat of entropy `log_s0` (log10 MJ/kg/K), until the mass inside falls to
        `mass_stop` or the radius to the radius stop.

        Returns:
            tuple: `(defect, centre)` (`radius_m` replaces the subject's radius, for the search of a model's own). `defect` is how
                far the model is from closing as a share of the radius: positive for a model too dense for its radius (it runs out of
                mass with radius left), negative for one too light; it is the radius left less the radius a uniform sphere of the
                mass left and the density there would have. `centre` is `(u, r, m, volume, central pressure)` at the stop.
        """
        mass, radius = self.subject.mass_kg, radius_m or self.subject.radius_m
        volume_of = self._volume_m3_per_kg
        edges = self.edges
        mass_stop, radius_stop = self.mass_stop, CENTRE_RADIUS_FRACTION * radius
        four_pi = 4.0 * math.pi
        lt_nodes = self.grid.adiabat_log_t(log_s0)
        layer = len(self.kinds) - 1
        u, r, m = math.log(ONE_BAR_PA), radius, mass
        if record is not None:
            record.append((u, r, m))

        def derivative(u_, r_, m_):
            volume = volume_of(u_, lt_nodes, layer)
            drdu = -math.exp(u_) * volume * r_ * r_ / (G * m_)
            return drdu, four_pi * r_ * r_ * drdu / volume

        def advance(h, k1):
            k2 = derivative(u + h / 2, r + h / 2 * k1[0], m + h / 2 * k1[1])
            k3 = derivative(u + h / 2, r + h / 2 * k2[0], m + h / 2 * k2[1])
            k4 = derivative(u + h, r + h * k3[0], m + h * k3[1])
            return r + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), m + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])

        limit = math.log(CENTRAL_PRESSURE_LIMIT_PA)
        while True:
            bottom = edges[layer - 1] if layer > 0 else 0.0
            target = max(bottom, mass_stop)
            k1 = derivative(u, r, m)
            h = min(STEP_LN_PA, STEP_FRACTION * r / abs(k1[0]), STEP_FRACTION * m / abs(k1[1]))
            r_new, m_new = advance(h, k1)
            if m_new > target and r_new > radius_stop:
                u, r, m = u + h, r_new, m_new
                if record is not None:
                    record.append((u, r, m))
                if u > limit:
                    break
                continue
            share_m = (m - target) / (m - m_new) if m_new <= target else 2.0
            share_r = (r - radius_stop) / (r - r_new) if r_new <= radius_stop else 2.0
            share = max(min(share_m, share_r, 1.0), 1.0e-6)
            r_new, m_new = advance(h * share, k1)
            u, r, m = u + h * share, r_new, m_new
            if record is not None:
                record.append((u, r, m))
            if share_m <= share_r and target == bottom and bottom > mass_stop:
                layer -= 1                      # on the top of the next layer down
                continue
            break
        volume = volume_of(u, lt_nodes, layer)
        r_sphere = (3.0 * m * volume / four_pi) ** (1.0 / 3.0)       # the sphere the mass left makes at the density here
        central = math.exp(u) + 2.0 * math.pi / 3.0 * G * r_sphere * r_sphere / (volume * volume)      # P_c of a uniform sphere
        return (r - r_sphere) / radius, (u, r, m, volume, central)


_LN10 = math.log(10.0)


def _temperature_at(model, log_s0, u):
    """The temperature (K) of the adiabat of entropy `log_s0` at ln P `u`."""
    grid = model.grid
    nodes = grid.adiabat_log_t(log_s0)
    y = min(max((u / _LN10 - 9.0 - grid.log_p0) / grid.d_log_p, 0.0), grid.n_p - 1.000001)
    j = int(y)
    return 10.0 ** (nodes[j] + (y - j) * (nodes[j + 1] - nodes[j]))


def _profile(model, log_s0, closed, evaluations, radius_m):
    record = []
    _defect, (u, _r, _m, _volume, central) = model.integrate(log_s0, record, radius_m)
    mass = model.subject.mass_kg
    ln_p = tuple(item[0] for item in record) + (math.log(central),)
    masses = tuple(item[2] / mass for item in record) + (0.0,)
    radii = tuple(item[1] / radius_m for item in record) + (0.0,)
    return Profile(mass, radius_m, model.grid.surface_temperature_k(log_s0), closed, central, _temperature_at(model, log_s0, u),
                   evaluations, ln_p, masses, radii)


def _own_radius_m(model, log_s0):
    """The radius (m) at which the model closes on the adiabat of entropy `log_s0`, searched from the subject's radius: the nearest the
    adiabats of `ADIABAT_RANGE_K` come to a radius they cannot close."""
    radius = model.subject.radius_m
    low, high = (radius, 4.0 * radius) if model.integrate(log_s0)[0] < 0.0 else (0.25 * radius, radius)
    for _ in range(40):
        middle = 0.5 * (low + high)
        if model.integrate(log_s0, None, middle)[0] > 0.0:       # left hollow: the radius is more than the mass fills
            high = middle
        else:
            low = middle
        if high - low < RADIUS_TOLERANCE * radius:
            break
    return 0.5 * (low + high)


def fit_adiabat(subject):
    """
    The profile of `subject` on the adiabat whose entropy makes the model close inside the subject's radius (see the module note).

    Raises:
        fetch.DataMissing: The hydrogen-helium tables were not downloaded.
        StructureError: The subject cannot describe a body.
    """
    model = _Model(subject)
    grid = model.grid
    low, high = (grid.log_entropy_of_surface_temperature(t) for t in ADIABAT_RANGE_K)
    seed_k = min(max(subject.temperature_1bar_k, ADIABAT_RANGE_K[0] * 1.5), ADIABAT_RANGE_K[1] / 1.5)
    seed = grid.log_entropy_of_surface_temperature(seed_k)
    step = 0.04
    s_lo, s_hi = max(seed - step, low), min(seed + step, high)
    f_lo, f_hi = model.integrate(s_lo)[0], model.integrate(s_hi)[0]
    evaluations = 2
    # the defect falls as the adiabat warms (a warmer body is less dense): widen until it brackets zero
    while f_lo < 0.0 and s_lo > low:
        s_hi, f_hi = s_lo, f_lo
        step *= 1.8
        s_lo = max(s_lo - step, low)
        f_lo = model.integrate(s_lo)[0]
        evaluations += 1
    while f_hi > 0.0 and s_hi < high:
        s_lo, f_lo = s_hi, f_hi
        step *= 1.8
        s_hi = min(s_hi + step, high)
        f_hi = model.integrate(s_hi)[0]
        evaluations += 1
    if f_lo < 0.0:                                            # too light for the radius even on the coldest adiabat
        return _profile(model, s_lo, False, evaluations, _own_radius_m(model, s_lo))
    if f_hi > 0.0:                                            # too dense for it even on the hottest
        return _profile(model, s_hi, False, evaluations, _own_radius_m(model, s_hi))
    side = 0                                                  # Illinois regula falsi on log S
    while evaluations < MAX_EVALUATIONS:
        x = (s_lo * f_hi - s_hi * f_lo) / (f_hi - f_lo)
        f = model.integrate(x)[0]
        evaluations += 1
        if abs(f) < RADIUS_TOLERANCE:
            return _profile(model, x, True, evaluations, subject.radius_m)
        if f > 0.0:
            s_lo, f_lo = x, f
            if side == 1:
                f_hi *= 0.5
            side = 1
        else:
            s_hi, f_hi = x, f
            if side == -1:
                f_lo *= 0.5
            side = -1
    return _profile(model, 0.5 * (s_lo + s_hi), True, evaluations, subject.radius_m)
