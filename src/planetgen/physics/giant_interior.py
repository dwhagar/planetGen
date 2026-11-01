# planetgen/physics/giant_interior.py

"""
The interior layers of a giant planet, an ice giant or a brown dwarf (GEN.642, docs/design/composition/giant-interiors.md)
==============================================================================================================

What a giant is made of, layer by layer from the top: a molecular hydrogen and helium envelope, a metallic hydrogen
envelope under it, a helium-rich zone where helium has rained out, a dilute ("fuzzy") core where the heavy elements are
mixed with metallic hydrogen, and a compact core of rock and iron; an ice giant has a molecular envelope over a mantle of
water, ammonia and methane over a rock and iron core. Each layer carries its pressures, its mass, the mass of each
material in it (rows of the shared materials table, GEN.427) and its state.

This module places and fills the layers from the quantities the radius and heavy-element laws (GEN.460) and the bulk
draw (GEN.428) already give: the mass, the radius, the heavy-element mass M_Z, the envelope metallicity, the natal helium
share, the one-bar temperature and the split of the heavy elements into iron, rock and ices. It draws only what those do
not decide, with `util.draw`: the molecular-to-metallic pressure, whether a compact core exists and its mass, and the mass
of the dilute core. It changes no mass, radius or fraction the caller fixed: the layer masses add up to the body's mass
and the layers hold exactly M_Z of heavy elements and the natal helium (tested to 1e-9).

How the layers are placed
-------------------------
The pressure at a layer's boundary is read from a hydrostatic profile at the mass inside it (`_structure`). With the downloaded
hydrogen-helium tables (GEN.733, `datasets/hhe_eos.py`) the profile is the structure of `giant_structure.fit_adiabat`: the
equations of hydrostatic balance integrated inward from one bar with the tables' density along the adiabat that closes the model
on the radius the caller fixed (a gas giant and a brown dwarf: a core and an envelope of the heavy elements, from M_Z and the
envelope metallicity; an ice giant: the rock core, the ice mantle and the envelope the draw placed). A run told to disregard
missing data (GEN.599) takes the profiles this note was written with: for a gas giant the n = 1 polytrope the repo already uses for
its central pressure (`P_c = pi G M^2 / (8 R^4)`, `P / P_c = (sin xi / xi)^2`, mass inside `(sin xi - xi cos xi) / pi`,
`xi = pi r / R`), for a brown dwarf the n = 3/2 polytrope (Lane-Emden, integrated here), and for an ice giant three layers of
constant density scaled to its radius (a stand-in, design note decision D7). A layer's mass fixes its pressures and not the other
way round, except for the two levels that are pressures: the molecular-to-metallic pressure `P_MM` (drawn 0.4 to 2 Mbar, D1) and
the onset of helium rain (2 Mbar, D2).

Where the heavy elements go (D5)
--------------------------------
Rock and iron fill the compact core first, then the dilute core, then the deep metallic envelope. Ice (water, ammonia,
methane), which the molecular envelope can hold as vapour, fills the molecular envelope first, then the metallic
envelope, then the dilute core. Whatever one material cannot fill, the other takes, so the layers always hold M_Z.

Left out (design note section 8): temperature, neon, helium rain in a brown dwarf. Each number
that has no source is named in `RANGES` and the table `giant-compact-core` with its state.

Use:

    subject = GiantSubject(mass_earth=317.8, radius_km=71492.0, heavy_mass_earth=42.0, envelope_z=0.04,
                           helium_share=0.2745, temperature_1bar_k=166.0, iron_share_of_refractory=0.29,
                           ice_shares={"water_liquid": 0.3})
    interior = draw_interior(subject)
    interior.layer("dilute_core").mass_earth
    interior.whole_body_shares()           # metal, rock, water, ice, envelope, envelope_helium, envelope_z
"""

import math
from dataclasses import dataclass, field

from planetgen.datasets import fetch, hhe_eos
from planetgen.generation import possibility_table as pt
from planetgen.physics import constants, giant_structure, unit_registry as units
from planetgen.util import draw

_mearth_to_kg = units.converter("Mearth", "kg")

# ---------------------------------------------------------------------------------------------------------------------
# What is known and what is a placeholder. Every entry names its source state (docs/design/composition/giant-interiors.md section 7).
# ---------------------------------------------------------------------------------------------------------------------

METALLIZATION_PRESSURE_RANGE_GPA = (40.0, 200.0)
"""tuple: The molecular-to-metallic hydrogen pressure, drawn log-uniform (D1). 0.4 to 0.6 Mbar from quantum Monte Carlo at
Jupiter's conditions ([R2966], READ), 1 to 2 Mbar from the older liquid-liquid transition work
([R1038], READ); the widest range wins (Boss's rule 3)."""

HELIUM_RAIN_ONSET_GPA = 200.0
"""float: Where helium becomes immiscible in metallic hydrogen for a proto-solar adiabat, 2 Mbar ([R2852: Figure 1 and 3], READ). The helium gradient starts here and the envelope above it keeps one depleted helium share."""

HELIUM_RAIN_POINTS_K = ((140.0, 0.26), (166.0, 0.85))
"""tuple: STAND-IN (D2, replaced by the [R3973] curve and a stored adiabat, GEN.641 follow-up): the
one-bar temperature and the molecular envelope's helium share against the natal share where it was measured. Saturn's
Y = 0.07 against 0.27 ([R2852], READ) and Jupiter's 0.234 against 0.2745 (Galileo, via Helled and Howard
2024, READ). A straight line through them, held at the Saturn value below it and at 1 from where it reaches 1."""

DEEP_ENVELOPE_MAX_Z = 0.9
"""float: STAND-IN: the most heavy elements the deep metallic envelope (with no dilute core, or over it) carries as a mass
fraction; what does not fit stays in the dilute core."""

MIN_ENVELOPE_FRACTION = 0.1
"""float: STAND-IN: a dilute core never takes more than all but this share of the mass."""

MIN_GRADIENT_FRACTION = 0.05
"""float: STAND-IN: when helium rain is active the layers between the rain level and the dilute core keep at least this share
of the mass, so the helium that rained out has somewhere to go."""

HELIUM_SHELL_THRESHOLD = 0.5
"""float: STAND-IN (D2): below this depletion half of the helium that rained out forms a helium shell above the core (the
Saturn picture, [R2852] section 2, READ for the existence, not for the share)."""

HELIUM_SHELL_FRACTION = 0.5
SHELL_HELIUM_SHARE = 0.9
"""float: STAND-IN (D2): the share of the removed helium that settles in the shell, and the shell's helium share of its
hydrogen and helium."""

COMPACT_CORE_RANGE_EARTH = (0.5, 20.0)
"""tuple: The compact core mass, log-uniform (D3). A few Earth masses at most for Jupiter ([R1895], READ),
15 to 18 for Saturn (via [R1893], 2ND-HAND); the widest range wins."""

DILUTE_CORE_RANGE_EARTH = (20.0, 100.0)
"""tuple: The dilute core mass (hydrogen and helium with heavy elements), log-uniform (D4). About 20 Earth masses where
evaporated solids mix in (2ND-HAND through [R1895]) to about 100 (Helled and Stevenson
2024, READ)."""

DILUTE_CORE_MAX_MASS_FRACTION = 0.65
"""float: Saturn's stable region is 55 of 95 Earth masses, 0.58 ([R2851], READ); the cap is a little above."""

DILUTE_CORE_MAX_HEAVY_FRACTION = 0.8
"""float: The heavy-element share of the dilute core's inner part reaches 0.8 ([R2851], READ)."""

DILUTE_CORE_MIN_BODY_EARTH = 30.0
"""float: Below this mass the heavy elements are a gradient through the whole body and there is no separate dilute core (D4)."""

ICE_GIANT_ENVELOPE_LIMIT = 0.5
"""float: A body whose hydrogen and helium are under this share of its mass is an ice giant ([R1892]: 10 to
20 percent for Uranus and Neptune)."""

ICE_GIANT_ROCK_FRACTION_RANGE = (0.05, 0.67)
"""tuple: The rock and iron share of an ice giant's heavy elements (D4.4, OPEN): 0.05 from the non-adiabatic empirical models, 0.67 from [R4880] (rock to water at most 2)."""

ICE_GIANT_RELATIVE_DENSITY = {"rock_core": 2.0, "ice_mantle": 1.0, "molecular_envelope": 0.25}
"""dict: STAND-IN (D7): the relative densities of the three constant-density layers of an ice giant. Only the ratios matter:
the scale is set so the layers fill the radius."""

SUPERIONIC_PRESSURE_GPA = 100.0
"""float: Water above this pressure is superionic ice (ice XVIII at 100 to 400 [R3081] abstract, READ)."""

# ---------------------------------------------------------------------------------------------------------------------

LAYER_KINDS = ("molecular_envelope", "metallic_envelope", "helium_gradient", "helium_shell", "dilute_core",
               "compact_core", "ice_mantle", "rock_core")
"""tuple: Every layer kind, the keys of the child table that stores them (DB.89)."""

INTERIOR_KINDS = ("gas_giant", "ice_giant", "brown_dwarf")
"""tuple: The kinds of interior (`GiantInterior.kind`), the values `planets.interior_kind` and `moons.interior_kind` hold (DB.89)."""

STATES = ("molecular_fluid", "metallic_fluid", "ionic_fluid", "superionic", "molten")
"""tuple: The states a layer can have."""

ENVELOPE_KINDS = ("molecular_envelope", "metallic_envelope", "helium_gradient", "helium_shell")
"""tuple: The layers whose heavy elements count as the envelope's metallicity; the heavy elements of every other layer are
counted in the whole-body metal, rock, water and ice fractions (the contract with GEN.428, section 5 of the design note)."""

HEAVY_MATERIALS = ("iron", "rock", "water_liquid", "ice_x", "nh3_ice", "ch4_ice")
GAS_MATERIALS = ("hydrogen", "helium")
WATER_MATERIALS = ("water_liquid", "ice_x")
OTHER_ICE_MATERIALS = ("nh3_ice", "ch4_ice")

MASS_TOLERANCE = 1e-9

COMPACT_CORE_TABLE = pt.REGISTRY.register(pt.Table(
    "giant-compact-core", 1,
    [pt.Row("compact core", (pt.ANY_CLASS,), pt.OPEN, "open", state=pt.OPEN_STATE,
            notes="A pure rock and iron core of 0.5 to 20 Earth masses. Formation models predict one (Helled and Stevenson "
                  "2024); gravity cannot see one under about 10 Earth masses; the share is a placeholder."),
     pt.Row("no compact core", (pt.ANY_CLASS,), pt.OPEN, "open", state=pt.OPEN_STATE,
            notes="Convective mixing or a giant impact removed it, or it never formed [R1895].")],
    description="Whether a giant has a compact core of heavy elements."))


class InteriorError(ValueError):
    """A subject or a set of layers that cannot describe a body."""


# ---------------------------------------------------------------------------------------------------------------------
# Polytropes
# ---------------------------------------------------------------------------------------------------------------------

def _lane_emden_table(index, steps=4000):
    """`(xi, theta, mass)` samples of the polytrope of `index`, mass as `xi^2 |theta'|`, out to the first zero of theta."""
    h = 1e-3
    xi, theta, phi = h, 1.0 - h * h / 6.0, -h / 3.0  # phi = theta'

    def derivative(x, t, p):
        return p, -(max(t, 0.0) ** index) - 2.0 * p / x

    rows = [(xi, theta, xi * xi * -phi)]
    while theta > 0.0 and len(rows) < steps * 10:
        k1 = derivative(xi, theta, phi)
        k2 = derivative(xi + h / 2, theta + h / 2 * k1[0], phi + h / 2 * k1[1])
        k3 = derivative(xi + h / 2, theta + h / 2 * k2[0], phi + h / 2 * k2[1])
        k4 = derivative(xi + h, theta + h * k3[0], phi + h * k3[1])
        new_theta = theta + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        new_phi = phi + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        if new_theta <= 0.0:
            # land on the zero by linear interpolation
            fraction = theta / (theta - new_theta)
            xi_zero = xi + fraction * h
            phi_zero = phi + fraction * (new_phi - phi)
            rows.append((xi_zero, 0.0, xi_zero * xi_zero * -phi_zero))
            return rows
        xi, theta, phi = xi + h, new_theta, new_phi
        rows.append((xi, theta, xi * xi * -phi))
    raise InteriorError(f"the polytrope of index {index} did not reach zero")


class Polytrope:
    """
    A polytrope scaled to a mass and a radius: pressure, radius fraction and mass inside as functions of each other.

    Args:
        index (float): The polytropic index, 1 for a giant planet and 1.5 for a degenerate brown dwarf.
        mass_kg, radius_m (float): The body's.
    """

    def __init__(self, index, mass_kg, radius_m):
        self.index = float(index)
        self.mass_kg = mass_kg
        self.radius_m = radius_m
        self._rows = _polytrope_rows(self.index)
        xi1, _theta, mass1 = self._rows[-1]
        self._xi1 = xi1
        self._mass1 = mass1
        # P_c = G M^2 / R^4 * xi1^4 / (4 pi (n + 1) (xi1^2 |theta'(xi1)|)^2)
        self.central_pressure_pa = (constants.G * mass_kg ** 2 / radius_m ** 4 * xi1 ** 4
                                    / (4.0 * math.pi * (self.index + 1.0) * mass1 ** 2))

    def _at_index(self, i):
        xi, theta, mass = self._rows[i]
        return xi / self._xi1, theta, mass / self._mass1

    def at_pressure(self, pressure_pa):
        """`(radius fraction, mass inside fraction)` where the pressure is `pressure_pa`; (1, 1) at or below the surface
        and (0, 0) at or above the centre."""
        if pressure_pa >= self.central_pressure_pa:
            return 0.0, 0.0
        if pressure_pa <= 0.0:
            return 1.0, 1.0
        target = (pressure_pa / self.central_pressure_pa) ** (1.0 / (self.index + 1.0))  # theta
        return self._interpolate(lambda row: row[1], target, descending=True)

    def at_mass_inside(self, fraction):
        """`(radius fraction, pressure in Pa)` where the mass inside is `fraction` of the whole."""
        if fraction >= 1.0:
            return 1.0, 0.0
        if fraction <= 0.0:
            return 0.0, self.central_pressure_pa
        radius, theta, _mass = self._interpolate_to(lambda row: row[2], fraction * self._mass1)
        return radius, self.central_pressure_pa * theta ** (self.index + 1.0)

    def pressure_at_mass_inside_pa(self, fraction):
        return self.at_mass_inside(fraction)[1]

    def _interpolate(self, key, target, descending):
        rows = self._rows
        lo, hi = 0, len(rows) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if (key(rows[mid]) > target) == descending:
                lo = mid
            else:
                hi = mid
        a, b = rows[lo], rows[hi]
        span = key(b) - key(a)
        t = 0.0 if span == 0.0 else (target - key(a)) / span
        xi = a[0] + t * (b[0] - a[0])
        mass = a[2] + t * (b[2] - a[2])
        return xi / self._xi1, mass / self._mass1

    def _interpolate_to(self, key, target):
        rows = self._rows
        lo, hi = 0, len(rows) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if key(rows[mid]) < target:
                lo = mid
            else:
                hi = mid
        a, b = rows[lo], rows[hi]
        span = key(b) - key(a)
        t = 0.0 if span == 0.0 else (target - key(a)) / span
        return (a[0] + t * (b[0] - a[0])) / self._xi1, a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])


_POLYTROPE_ROWS = {}


def _polytrope_rows(index):
    if index not in _POLYTROPE_ROWS:
        _POLYTROPE_ROWS[index] = _lane_emden_table(index)
    return _POLYTROPE_ROWS[index]




# ---------------------------------------------------------------------------------------------------------------------
# The result
# ---------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Layer:
    """
    One layer.

    Attributes:
        kind (str): One of `LAYER_KINDS`.
        state (str): One of `STATES`.
        top_pressure_gpa, bottom_pressure_gpa (float): Pressure at the layer's top and bottom.
        mass_earth (float): The layer's mass in Earth masses.
        components (dict): Each material's mass in the layer in Earth masses (rows of the shared materials table); they add
            up to `mass_earth`.
        top_radius_fraction, bottom_radius_fraction (float): Where the layer's top and bottom lie as a fraction of the
            body's radius (1 at the 1 bar surface, 0 at the centre).
    """

    kind: str
    state: str
    top_pressure_gpa: float
    bottom_pressure_gpa: float
    mass_earth: float
    components: dict
    top_radius_fraction: float = 1.0
    bottom_radius_fraction: float = 0.0

    def component_fractions(self):
        """Each material's share of the layer's mass."""
        return {name: mass / self.mass_earth for name, mass in self.components.items()} if self.mass_earth > 0 else {}

    def heavy_mass_earth(self):
        return math.fsum(m for name, m in self.components.items() if name in HEAVY_MATERIALS)

    def helium_mass_earth(self):
        return self.components.get("helium", 0.0)


@dataclass(frozen=True)
class GiantInterior:
    """
    The layers of one body, outside in, and the numbers that decided them.

    Attributes:
        layers (tuple): `Layer`s from the top down.
        total_mass_earth (float): The body's mass.
        kind (str): `gas_giant`, `ice_giant` or `brown_dwarf`.
        metallization_pressure_gpa (float): The drawn molecular-to-metallic pressure.
        helium_depletion (float): The molecular envelope's helium share against the natal share; 1 with no rain.
    """

    layers: tuple
    total_mass_earth: float
    kind: str
    metallization_pressure_gpa: float
    helium_depletion: float

    def layer(self, kind):
        """The layer of `kind`, or `None`."""
        for layer in self.layers:
            if layer.kind == kind:
                return layer
        return None

    @property
    def core_mass_earth(self):
        """The compact core's mass (or an ice giant's rock core), 0 with none."""
        compact = self.layer("compact_core") or self.layer("rock_core")
        return compact.mass_earth if compact else 0.0

    @property
    def dilute_core_mass_earth(self):
        layer = self.layer("dilute_core")
        return layer.mass_earth if layer else 0.0

    @property
    def has_helium_shell(self):
        return self.layer("helium_shell") is not None

    def heavy_mass_earth(self):
        return math.fsum(layer.heavy_mass_earth() for layer in self.layers)

    def helium_mass_earth(self):
        return math.fsum(layer.helium_mass_earth() for layer in self.layers)

    def whole_body_shares(self):
        """
        The totals GEN.428 stores for the body (`bulk_composition.BulkFractions`): the heavy elements of the core layers
        as `metal`, `rock`, `water` and `ice` shares of the whole body, everything else (all the hydrogen and helium, and the
        heavy elements of the envelope layers) as `envelope`, with the envelope's `envelope_helium` and `envelope_z`.
        """
        total = self.total_mass_earth
        iron = rock = water = ice = 0.0
        envelope_heavy = 0.0
        for layer in self.layers:
            for name, mass in layer.components.items():
                if name not in HEAVY_MATERIALS:
                    continue
                if layer.kind in ENVELOPE_KINDS:
                    envelope_heavy += mass
                elif name == "iron":
                    iron += mass
                elif name == "rock":
                    rock += mass
                elif name in WATER_MATERIALS:
                    water += mass
                else:
                    ice += mass
        envelope = total - iron - rock - water - ice
        return {"metal": iron / total, "rock": rock / total, "water": water / total, "ice": ice / total,
                "envelope": envelope / total, "envelope_helium": self.helium_mass_earth() / envelope,
                "envelope_z": envelope_heavy / envelope}

    def validate(self, heavy_mass_earth=None, helium_share=None):
        """Checks the invariants of design note section 6: layer masses add up to the body, pressures never fall going
        down, each layer's materials add up to its mass, and (given) the layers hold the stated heavy elements and the
        natal helium. Raises `InteriorError`."""
        if abs(math.fsum(layer.mass_earth for layer in self.layers) - self.total_mass_earth) > MASS_TOLERANCE * self.total_mass_earth:
            raise InteriorError("the layer masses do not add up to the body's mass")
        previous = 0.0
        for layer in self.layers:
            if layer.kind not in LAYER_KINDS or layer.state not in STATES:
                raise InteriorError(f"layer {layer.kind!r} state {layer.state!r} is unknown")
            if layer.top_pressure_gpa < previous - 1e-6 * max(previous, 1.0) or layer.bottom_pressure_gpa < layer.top_pressure_gpa - 1e-9:
                raise InteriorError(f"the pressures of layer {layer.kind!r} do not rise downward")
            previous = layer.bottom_pressure_gpa
            if abs(math.fsum(layer.components.values()) - layer.mass_earth) > MASS_TOLERANCE * max(layer.mass_earth, 1.0):
                raise InteriorError(f"the materials of layer {layer.kind!r} do not add up to its mass")
            if any(m < -MASS_TOLERANCE for m in layer.components.values()):
                raise InteriorError(f"layer {layer.kind!r} has a negative component")
        if heavy_mass_earth is not None and abs(self.heavy_mass_earth() - heavy_mass_earth) > MASS_TOLERANCE * self.total_mass_earth:
            raise InteriorError("the layers do not hold the stated heavy-element mass")
        if helium_share is not None:
            expected = helium_share * (self.total_mass_earth - self.heavy_mass_earth())
            if abs(self.helium_mass_earth() - expected) > MASS_TOLERANCE * self.total_mass_earth:
                raise InteriorError("the layers do not hold the natal helium")


BODY_FIELDS = ("interior_kind", "core_mass_earth", "dilute_core_mass_earth", "metallization_pressure_gpa", "helium_depletion",
               "has_helium_shell")
"""tuple: The summary attributes (and `planets`/`moons` columns, DB.89) of a giant's interior; `None` on a body that has none."""

LAYER_FIELDS = ("kind", "state", "top_pressure_gpa", "bottom_pressure_gpa", "mass_earth", "top_radius_fraction",
                "bottom_radius_fraction")
"""tuple: A stored layer's scalar fields; its materials (`HEAVY_MATERIALS` and `GAS_MATERIALS`, Earth masses) follow."""

LAYERS_ATTRIBUTE = "interior_layers"
"""str: The body attribute that holds the layers as a list of flat dicts (one per layer, outside in) and the child table that stores them."""


def body_values(interior):
    """
    The attributes a body holds for `interior` (DB.89): `BODY_FIELDS`, and `interior_layers`, a list of flat dicts (one per
    layer, outside in) with `LAYER_FIELDS` and the mass of every material (0.0 where the layer has none), so a save and a
    load return them unchanged.
    """
    layers = [{**{name: getattr(layer, name) for name in LAYER_FIELDS},
               **{material: layer.components.get(material, 0.0) for material in (*GAS_MATERIALS, *HEAVY_MATERIALS)}}
              for layer in interior.layers]
    return {"interior_kind": interior.kind, "core_mass_earth": interior.core_mass_earth,
            "dilute_core_mass_earth": interior.dilute_core_mass_earth,
            "metallization_pressure_gpa": interior.metallization_pressure_gpa, "helium_depletion": interior.helium_depletion,
            "has_helium_shell": interior.has_helium_shell, LAYERS_ATTRIBUTE: layers}


def clear(body):
    """Removes a stored interior from `body` (a body regenerated as a class with none)."""
    for name in (*BODY_FIELDS, LAYERS_ATTRIBUTE):
        setattr(body, name, None)


def interior_of(body):
    """
    The `GiantInterior` a body's stored attributes describe (`body_values` read back), or `None` for a body with none.
    """
    rows = getattr(body, LAYERS_ATTRIBUTE, None)
    if not rows or getattr(body, "interior_kind", None) is None:
        return None
    layers = tuple(Layer(
        kind=row["kind"], state=row["state"], top_pressure_gpa=row["top_pressure_gpa"],
        bottom_pressure_gpa=row["bottom_pressure_gpa"], mass_earth=row["mass_earth"],
        components={material: row[material] for material in (*GAS_MATERIALS, *HEAVY_MATERIALS) if row[material] > 0.0},
        top_radius_fraction=row["top_radius_fraction"], bottom_radius_fraction=row["bottom_radius_fraction"]) for row in rows)
    return GiantInterior(layers, math.fsum(layer.mass_earth for layer in layers), body.interior_kind,
                         body.metallization_pressure_gpa, body.helium_depletion)


@dataclass(frozen=True)
class GiantSubject:
    """
    What the interior is built from (all fixed by the radius and bulk draws before this module runs).

    Attributes:
        mass_earth, radius_km (float): The body's mass and radius (GEN.429, GEN.460).
        heavy_mass_earth (float): M_Z, all elements heavier than helium (GEN.460), less than the mass.
        envelope_z (float): The heavy-element mass fraction of the molecular envelope (the atmosphere's metallicity).
        helium_share (float): Y / (X + Y) of the natal gas.
        temperature_1bar_k (float): The one-bar temperature, which decides helium rain.
        iron_share_of_refractory (float): The iron share of the rock and iron (the star's, as for a terrestrial planet).
        ice_shares (dict): The ices' shares of the heavy elements when the body formed beyond the snow line, by material
            (`water_liquid`, `nh3_ice`, `ch4_ice`); empty when it formed inside it. They add up to less than 1.
        is_brown_dwarf (bool): A brown dwarf: one envelope over no core, on the n = 3/2 polytrope.
    """

    mass_earth: float
    radius_km: float
    heavy_mass_earth: float
    envelope_z: float
    helium_share: float
    temperature_1bar_k: float
    iron_share_of_refractory: float = 0.3
    ice_shares: dict = field(default_factory=dict)
    is_brown_dwarf: bool = False

    def __post_init__(self):
        if not (self.mass_earth > 0 and self.radius_km > 0):
            raise InteriorError("a giant needs a positive mass and radius")
        if not (0.0 <= self.heavy_mass_earth < self.mass_earth):
            raise InteriorError("the heavy-element mass must be under the body's mass")
        if not (0.0 <= self.envelope_z < 1.0 and 0.0 < self.helium_share < 1.0 and 0.0 <= self.iron_share_of_refractory <= 1.0):
            raise InteriorError("an envelope metallicity, helium share or iron share is outside its range")
        if any(v < 0 for v in self.ice_shares.values()) or sum(self.ice_shares.values()) >= 1.0 + 1e-12:
            raise InteriorError("the ice shares of the heavy elements must be positive and add up to at most 1")

    @property
    def hydrogen_helium_fraction(self):
        return (self.mass_earth - self.heavy_mass_earth) / self.mass_earth


# ---------------------------------------------------------------------------------------------------------------------
# Helium
# ---------------------------------------------------------------------------------------------------------------------

def helium_depletion(temperature_1bar_k):
    """The molecular envelope's helium share against the natal share (D2, a stand-in): 1 above the temperature where the
    line through the two measured points reaches 1, falling to the cold point's value below it."""
    (t_cold, y_cold), (t_warm, y_warm) = HELIUM_RAIN_POINTS_K
    value = y_cold + (y_warm - y_cold) * (temperature_1bar_k - t_cold) / (t_warm - t_cold)
    return min(1.0, max(y_cold, value))


# ---------------------------------------------------------------------------------------------------------------------
# The draw
# ---------------------------------------------------------------------------------------------------------------------

def _log_uniform(low, high):
    return math.exp(draw.uniform(math.log(low), math.log(high)))


def _heavy_composition(subject, refractory_share):
    """The materials of the heavy elements as shares of M_Z: `{iron, rock, water_liquid, nh3_ice, ch4_ice}` where the ices
    take `1 - refractory_share` in the proportions of `subject.ice_shares`."""
    ice_total = 1.0 - refractory_share
    shares = {"iron": refractory_share * subject.iron_share_of_refractory,
              "rock": refractory_share * (1.0 - subject.iron_share_of_refractory)}
    weights = {name: value for name, value in subject.ice_shares.items() if value > 0}
    weight_total = math.fsum(weights.values())
    for name in ("water_liquid", "nh3_ice", "ch4_ice"):
        shares[name] = ice_total * weights.get(name, 0.0) / weight_total if weight_total > 0 else 0.0
    return shares


class _Bin:
    """The heavy elements still to place: refractory (iron and rock) and ice, in Earth masses, taken by layers."""

    def __init__(self, subject, refractory_share):
        shares = _heavy_composition(subject, refractory_share)
        total = subject.heavy_mass_earth
        self.iron = shares["iron"] * total
        self.rock = shares["rock"] * total
        self.ices = {name: shares[name] * total for name in ("water_liquid", "nh3_ice", "ch4_ice")}

    @property
    def refractory(self):
        return self.iron + self.rock

    @property
    def ice(self):
        return math.fsum(self.ices.values())

    @property
    def total(self):
        return self.refractory + self.ice

    def take(self, amount, prefer):
        """Takes `amount` of heavy elements, `prefer` ("refractory" or "ice") first; returns `{material: mass}`."""
        amount = min(amount, self.total)
        order = ("refractory", "ice") if prefer == "refractory" else ("ice", "refractory")
        taken = {}
        for group in order:
            available = self.refractory if group == "refractory" else self.ice
            part = min(amount, available)
            if part <= 0:
                continue
            if group == "refractory":
                share = part / self.refractory
                for name in ("iron", "rock"):
                    value = getattr(self, name)
                    taken[name] = taken.get(name, 0.0) + value * share
                    setattr(self, name, value - value * share)
            else:
                share = part / self.ice
                for name in tuple(self.ices):
                    value = self.ices[name]
                    taken[name] = taken.get(name, 0.0) + value * share
                    self.ices[name] = value - value * share
            amount -= part
        return taken


def _layer(kind, state, top_pa, bottom_pa, heavy, helium, hydrogen, top_radius=1.0, bottom_radius=0.0):
    components = dict(heavy)
    if helium:
        components["helium"] = helium
    if hydrogen:
        components["hydrogen"] = hydrogen
    components = {name: value for name, value in components.items() if value > 0.0}
    mass = math.fsum(components.values())
    return Layer(kind, state, top_pa / 1e9, bottom_pa / 1e9, mass, components, top_radius, bottom_radius)


def _water_material(bottom_pressure_gpa):
    return "ice_x" if bottom_pressure_gpa >= SUPERIONIC_PRESSURE_GPA else "water_liquid"


STAND_IN_POLYTROPE = "giant layer pressures from the polytrope"
"""str: What the run records when the hydrogen-helium tables were not downloaded and the profile is the polytrope's (GEN.599)."""


def _structure(subject, polytrope_index):
    """
    The hydrostatic profile layer pressures are read from (GEN.733): the structure of `giant_structure` from the downloaded
    hydrogen-helium tables, closed on the subject's radius, or, where the tables are missing and the run disregards that, the
    polytrope of `polytrope_index` (the profile before GEN.733).
    """
    mass_kg = _mearth_to_kg(subject.mass_earth)
    radius_m = subject.radius_km * 1e3
    if fetch.locate(hhe_eos.DATASET).use_stand_in(STAND_IN_POLYTROPE):
        return Polytrope(polytrope_index, mass_kg, radius_m)
    heavy_kg = subject.heavy_mass_earth * constants.EARTH_MASS_TO_KG
    if subject.is_brown_dwarf:
        z, core, env_ice, core_ice = subject.heavy_mass_earth / subject.mass_earth, 0.0, 0.0, 0.0
    else:
        z = subject.envelope_z
        core = max((heavy_kg - z * mass_kg) / (1.0 - z), 0.0)
        env_heavy = heavy_kg - core if core > 0.0 else heavy_kg
        ice_kg = sum(subject.ice_shares.values()) * heavy_kg
        env_ice = min(1.0, ice_kg / env_heavy) if env_heavy > 0.0 else 0.0
        core_ice = min(max((ice_kg - env_ice * env_heavy) / core, 0.0), 1.0) if core > 0.0 else 0.0
    return giant_structure.fit_adiabat(giant_structure.StructureSubject(
        mass_kg, radius_m, subject.helium_share, envelope_z=z, heavy_mass_kg=heavy_kg, envelope_ice_fraction=env_ice,
        core_ice_fraction=core_ice, iron_share=subject.iron_share_of_refractory, temperature_1bar_k=subject.temperature_1bar_k))


def draw_interior(subject):
    """
    The layers of `subject`, drawn with the ambient stream (`util.draw`). Gas giants, ice giants and brown dwarfs are told
    apart by the subject (design note section 3). The draw order is fixed (a change is an epoch change): P_MM, then the
    compact core's presence and mass, then the dilute core's mass.

    Returns:
        GiantInterior: validated against the subject's heavy-element mass and natal helium.
    """
    if subject.is_brown_dwarf:
        interior = _brown_dwarf(subject)
    elif subject.hydrogen_helium_fraction < ICE_GIANT_ENVELOPE_LIMIT:
        interior = _ice_giant(subject)
    else:
        interior = _gas_giant(subject)
    interior.validate(subject.heavy_mass_earth, subject.helium_share)
    return interior


def _gas_giant(subject):
    mass = subject.mass_earth
    profile = _structure(subject, 1.0)
    p_mm_gpa = _log_uniform(*METALLIZATION_PRESSURE_RANGE_GPA)
    ratio_mm = profile.at_pressure(p_mm_gpa * 1e9)[1]  # mass inside the P_MM level
    ratio_rain = profile.at_pressure(HELIUM_RAIN_ONSET_GPA * 1e9)[1]

    # The compact core
    row = COMPACT_CORE_TABLE.draw("giant").row.name
    core_mass = 0.0
    if row == "compact core":
        core_mass = _log_uniform(*COMPACT_CORE_RANGE_EARTH)
    # The dilute core
    dilute_drawn = _log_uniform(*DILUTE_CORE_RANGE_EARTH)

    heavy = subject.heavy_mass_earth
    ice_available = sum(subject.ice_shares.values())
    bin_ = _Bin(subject, refractory_share=1.0 - ice_available)

    core_mass = min(core_mass, 0.8 * heavy, bin_.refractory)

    mass_mol = mass * (1.0 - ratio_mm)
    mass_to_rain = mass * (1.0 - ratio_rain)           # everything above the helium rain level
    z_out = min(subject.envelope_z, 0.9 * (heavy - core_mass) / (mass - core_mass))

    # The dilute core: at least what the heavy elements need at DILUTE_CORE_MAX_HEAVY_FRACTION, at most the cap
    depletion_now = helium_depletion(subject.temperature_1bar_k)
    dilute_cap = max(0.0, DILUTE_CORE_MAX_MASS_FRACTION * mass - core_mass)
    excess = heavy - z_out * (mass - core_mass) - core_mass     # heavy elements left for the dilute core above the envelope's
    dilute = 0.0
    if excess > 0.0:
        # A small body draws a dilute core only when its heavy elements have nowhere else to go
        needed = excess / max(DILUTE_CORE_MAX_HEAVY_FRACTION - z_out, 1e-6)
        dilute = min(dilute_cap, max(dilute_drawn if mass >= DILUTE_CORE_MIN_BODY_EARTH else 0.0, needed))
        dilute = min(dilute, mass - core_mass - MIN_ENVELOPE_FRACTION * mass)
        if depletion_now < 1.0:  # leave the helium gradient region its room above the dilute core
            dilute = max(min(dilute, mass - core_mass - max(mass_mol, mass_to_rain) - MIN_GRADIENT_FRACTION * mass),
                         min(needed, dilute_cap))

    if mass_mol > mass - core_mass - dilute:  # a body with little pressure to spare metallizes higher up
        mass_mol = mass - core_mass - dilute
        p_mm_gpa = profile.pressure_at_mass_inside_pa(1.0 - mass_mol / mass) / 1e9

    # Helium
    depletion = helium_depletion(subject.temperature_1bar_k)
    y_nat = subject.helium_share
    rain_active = depletion < 1.0

    # Layer masses (Earth masses, hydrogen, helium and heavy elements together)
    metallic_total = mass - core_mass - dilute - mass_mol
    if metallic_total < 0.0:
        raise InteriorError("the layers do not fit in the body")
    shell_mass = 0.0
    gradient_mass = 0.0
    upper_metallic = metallic_total
    if rain_active and mass_to_rain > mass_mol:
        upper_metallic = min(metallic_total, mass_to_rain - mass_mol)
        gradient_mass = metallic_total - upper_metallic
    elif rain_active:
        gradient_mass = metallic_total
        upper_metallic = 0.0

    # Heavy elements: the envelope layers hold z_out of their mass, the dilute core takes what is left.
    # What the dilute core cannot hold (or all of it, with no dilute core) enriches the deep envelope, never the atmosphere.
    deep_mass = upper_metallic + gradient_mass
    z_deep = z_out
    if deep_mass > 0.0:
        to_place = heavy - core_mass - DILUTE_CORE_MAX_HEAVY_FRACTION * dilute - z_out * mass_mol
        z_deep = min(max(z_out, to_place / deep_mass), DEEP_ENVELOPE_MAX_Z)
    env_heavy = {"molecular_envelope": z_out * mass_mol, "metallic_envelope": z_deep * upper_metallic,
                 "helium_gradient": z_deep * gradient_mass}
    env_total = math.fsum(env_heavy.values())
    dilute_heavy = heavy - core_mass - env_total
    if dilute_heavy < -1e-12:
        raise InteriorError("the envelope holds more heavy elements than M_Z")
    if dilute == 0.0 and dilute_heavy > 0.0:
        # no separate dilute core: the deep metallic envelope carries the rest
        if gradient_mass > 0.0:
            env_heavy["helium_gradient"] += dilute_heavy
        elif upper_metallic > 0.0:
            env_heavy["metallic_envelope"] += dilute_heavy
        else:
            env_heavy["molecular_envelope"] += dilute_heavy
        dilute_heavy = 0.0
    dilute_heavy = max(dilute_heavy, 0.0)

    # Place the heavy elements: refractory to the compact core first, ice to the envelope first
    core_take = bin_.take(core_mass, "refractory")
    layer_heavy = {}
    for name in ("molecular_envelope", "metallic_envelope", "helium_gradient"):
        layer_heavy[name] = bin_.take(env_heavy[name], "ice")
    layer_heavy["dilute_core"] = bin_.take(dilute_heavy, "refractory")
    leftover = bin_.take(bin_.total, "refractory")
    for name, value in leftover.items():  # rounding dust only
        target = layer_heavy["dilute_core"] if dilute > 0 else layer_heavy["molecular_envelope"]
        target[name] = target.get(name, 0.0) + value

    # Helium: natal share everywhere, the molecular envelope and the envelope above the rain level depleted
    def gas_of(layer_mass, heavy_mass, y_share):
        hhe = layer_mass - heavy_mass
        return hhe * y_share, hhe * (1.0 - y_share)

    def heavy_of(kind):
        return math.fsum(layer_heavy.get(kind, {}).values())

    gradient_hhe = gradient_mass - heavy_of("helium_gradient")
    hhe_above = (mass_mol - heavy_of("molecular_envelope")) + (upper_metallic - heavy_of("metallic_envelope"))
    y_out = depletion * y_nat
    if rain_active and hhe_above > 0.0:
        # The rained-out helium has to fit in the gradient and shell below, no more than SHELL_HELIUM_SHARE of their gas:
        # a thin region can hold only so much, so the envelope above keeps the rest.
        dilute_hhe = dilute - heavy_of("dilute_core")
        capacity_gradient = max(0.0, (SHELL_HELIUM_SHARE - y_nat) * gradient_hhe)
        capacity_dilute = max(0.0, (SHELL_HELIUM_SHARE - y_nat) * dilute_hhe)
        y_out = max(y_out, y_nat - (capacity_gradient + capacity_dilute) / hhe_above)
        depletion = y_out / y_nat
    gas = {}
    for kind, layer_mass in (("molecular_envelope", mass_mol), ("metallic_envelope", upper_metallic)):
        gas[kind] = gas_of(layer_mass, heavy_of(kind), y_out)
    removed = (y_nat - y_out) * hhe_above
    to_dilute = 0.0  # the rained-out helium sinks into the gradient region first and fills the dilute core's gas with the rest
    if removed > 0.0:
        to_dilute = removed * capacity_dilute / (capacity_gradient + capacity_dilute)
        removed -= to_dilute
    dilute_hhe = dilute - heavy_of("dilute_core")
    gas["dilute_core"] = (y_nat * dilute_hhe + to_dilute, dilute_hhe - y_nat * dilute_hhe - to_dilute)
    shell_helium = 0.0
    shell_hydrogen = 0.0
    if rain_active and depletion < HELIUM_SHELL_THRESHOLD and gradient_hhe > 0.0:
        shell_helium = HELIUM_SHELL_FRACTION * removed
        shell_mass = min(shell_helium / SHELL_HELIUM_SHARE, 0.8 * gradient_hhe)
        shell_helium = shell_mass * SHELL_HELIUM_SHARE
        shell_hydrogen = shell_mass - shell_helium
    gradient_helium = y_nat * gradient_hhe + removed - shell_helium
    gradient_helium = min(gradient_helium, gradient_hhe - shell_mass)
    if gradient_mass > 0.0:
        # the shell is carved out of the gradient region's hydrogen and helium
        gas["helium_gradient"] = (gradient_helium, gradient_hhe - shell_mass - gradient_helium)
    # Anything the cap could not place stays natal-conserving: fix the total by returning it to the dilute core
    total_helium = math.fsum(g[0] for g in gas.values()) + shell_helium
    target_helium = y_nat * (mass - heavy)
    correction = target_helium - total_helium
    if abs(correction) > 0.0:
        sink = "dilute_core" if dilute > 0 else ("helium_gradient" if gradient_mass > 0 else "metallic_envelope" if upper_metallic > 0 else "molecular_envelope")
        h, hy = gas[sink]
        gas[sink] = (h + correction, hy - correction)

    # Pressures from the polytrope at the mass inside each boundary
    layers = []
    above = 0.0
    rows = [("molecular_envelope", mass_mol, "molecular_fluid"), ("metallic_envelope", upper_metallic, "metallic_fluid"),
            ("helium_gradient", gradient_mass - shell_mass, "metallic_fluid"), ("helium_shell", shell_mass, "metallic_fluid"),
            ("dilute_core", dilute, "metallic_fluid"), ("compact_core", core_mass, "molten")]
    for kind, layer_mass, state in rows:
        if layer_mass <= 0.0:
            continue
        top_radius, top = profile.at_mass_inside(1.0 - above / mass)
        bottom_radius, bottom = profile.at_mass_inside(1.0 - (above + layer_mass) / mass)
        if kind == "compact_core":
            heavy_part, helium, hydrogen = core_take, 0.0, 0.0
        elif kind == "helium_shell":
            heavy_part, helium, hydrogen = {}, shell_helium, shell_hydrogen
        else:
            heavy_part = layer_heavy.get(kind, {})
            helium, hydrogen = gas[kind][0], gas[kind][1]
        top = max(top, 0.0)
        layers.append(_layer(kind, state, top if kind != "molecular_envelope" else 1e5, bottom, heavy_part, helium, hydrogen,
                             top_radius, bottom_radius))
        above += layer_mass
    layers = _fix_water(layers)
    return GiantInterior(tuple(layers), mass, "gas_giant", p_mm_gpa, depletion)


def _fix_water(layers):
    """Water at or below 100 GPa is ordinary liquid, above it superionic ice (the materials table's `ice_x`)."""
    fixed = []
    for layer in layers:
        components = dict(layer.components)
        if "water_liquid" in components and layer.bottom_pressure_gpa >= SUPERIONIC_PRESSURE_GPA and layer.kind in (
                "ice_mantle", "dilute_core"):
            components["ice_x"] = components.pop("water_liquid")
        fixed.append(Layer(layer.kind, layer.state, layer.top_pressure_gpa, layer.bottom_pressure_gpa, layer.mass_earth,
                           components, layer.top_radius_fraction, layer.bottom_radius_fraction))
    return fixed


class _ConstantDensityProfile:
    """Pressure against radius for layers of constant density stacked from the centre, integrated inward from a surface
    pressure of 1 bar (the stand-in closure of an ice giant, design note D7). `layers` is `[(mass_kg, density)]` from the
    centre out; the radii follow from the masses and densities, so the outer radius is the sum of the layer thicknesses."""

    def __init__(self, layers, shells_per_layer=300):
        self.outer_radii = []
        self.inner_mass = []
        radius = 0.0
        mass = 0.0
        for layer_mass, density in layers:
            radius = (radius ** 3 + 3.0 * layer_mass / (4.0 * math.pi * density)) ** (1.0 / 3.0)
            mass += layer_mass
            self.outer_radii.append(radius)
            self.inner_mass.append(mass)
        self.radius = radius
        self.density = [density for _mass, density in layers]
        self.layer_masses = [layer_mass for layer_mass, _density in layers]
        self.samples = []  # (radius, pressure, mass inside), surface first
        pressure = 1e5
        for index in range(len(layers) - 1, -1, -1):
            r_out = self.outer_radii[index]
            r_in = self.outer_radii[index - 1] if index > 0 else 0.0
            m_below = self.inner_mass[index - 1] if index > 0 else 0.0
            rho = self.density[index]
            dr = (r_out - r_in) / shells_per_layer
            self.samples.append((r_out, pressure, self.inner_mass[index]))
            for step in range(shells_per_layer):
                r_mid = r_out - (step + 0.5) * dr
                m_mid = m_below + 4.0 / 3.0 * math.pi * rho * (r_mid ** 3 - r_in ** 3)
                pressure += rho * constants.G * m_mid / (r_mid ** 2) * dr
                self.samples.append((r_out - (step + 1) * dr, pressure, m_mid))
        self.central_pressure_pa = pressure

    def boundary_pressures(self, index):
        """`(top, bottom)` pressure of layer `index` (0 is the innermost)."""
        r_out = self.outer_radii[index]
        r_in = self.outer_radii[index - 1] if index > 0 else 0.0
        return self.pressure_at(r_out), self.pressure_at(r_in)

    def pressure_at(self, radius):
        if radius <= 0.0:
            return self.central_pressure_pa
        for (r_a, p_a, _ma), (r_b, p_b, _mb) in zip(self.samples, self.samples[1:]):
            if r_b <= radius <= r_a:
                t = 0.0 if r_a == r_b else (r_a - radius) / (r_a - r_b)
                return p_a + t * (p_b - p_a)
        return self.samples[0][1]

    def radius_at_pressure(self, pressure):
        """The radius where the pressure is `pressure`, or `None` outside the body's range."""
        if not self.samples[0][1] <= pressure <= self.central_pressure_pa:
            return None
        for (r_a, p_a, _ma), (r_b, p_b, _mb) in zip(self.samples, self.samples[1:]):
            if p_a <= pressure <= p_b:
                t = 0.0 if p_a == p_b else (pressure - p_a) / (p_b - p_a)
                return r_a + t * (r_b - r_a)
        return None

    def mass_inside(self, radius):
        for index, r_out in enumerate(self.outer_radii):
            if radius <= r_out:
                r_in = self.outer_radii[index - 1] if index > 0 else 0.0
                m_below = self.inner_mass[index - 1] if index > 0 else 0.0
                return m_below + 4.0 / 3.0 * math.pi * self.density[index] * (radius ** 3 - r_in ** 3)
        return self.inner_mass[-1]


class _LayeredProfile:
    """The structure of `giant_structure` seen the way `_ConstantDensityProfile` is: layers from the centre out, `masses` their masses
    (kg), with the pressures at their boundaries and the radius and mass at a pressure, so the ice giant reads either."""

    def __init__(self, structure, masses):
        self.structure = structure
        self.radius = structure.radius_m
        self.layer_masses = list(masses)
        self.inner_mass = [math.fsum(masses[:i + 1]) for i in range(len(masses))]
        self.outer_radii = [structure.at_mass_inside(m / structure.mass_kg)[0] * self.radius for m in self.inner_mass]
        self.central_pressure_pa = structure.central_pressure_pa

    def boundary_pressures(self, index):
        """`(top, bottom)` pressure of layer `index` (0 is the innermost)."""
        top = giant_structure.ONE_BAR_PA if index == len(self.inner_mass) - 1 else \
            self.structure.pressure_at_mass_inside_pa(self.inner_mass[index] / self.structure.mass_kg)
        bottom = self.structure.pressure_at_mass_inside_pa(self.inner_mass[index - 1] / self.structure.mass_kg) if index > 0 \
            else self.central_pressure_pa
        return top, bottom

    def radius_at_pressure(self, pressure):
        """The radius where the pressure is `pressure`, or `None` outside the body's range."""
        if not giant_structure.ONE_BAR_PA <= pressure <= self.central_pressure_pa:
            return None
        return self.structure.at_pressure(pressure)[0] * self.radius

    def mass_inside(self, radius):
        return self.structure.at_radius(radius / self.radius)[0] * self.structure.mass_kg


def _ice_giant(subject):
    mass = subject.mass_earth
    heavy = subject.heavy_mass_earth
    z_out = min(subject.envelope_z, 0.5)
    envelope = (mass - heavy) / (1.0 - z_out)
    p_mm_gpa = _log_uniform(*METALLIZATION_PRESSURE_RANGE_GPA)
    rock_fraction = draw.uniform(*ICE_GIANT_ROCK_FRACTION_RANGE) if subject.ice_shares else 1.0
    bin_ = _Bin(subject, refractory_share=rock_fraction)
    envelope_take = bin_.take(z_out * envelope, "ice")
    mantle_take = bin_.take(bin_.ice, "ice")
    core_take = bin_.take(bin_.total, "refractory")
    mantle = math.fsum(mantle_take.values())
    core = math.fsum(core_take.values())

    kg = constants.EARTH_MASS_TO_KG
    names = [kind for kind, m in (("rock_core", core), ("ice_mantle", mantle), ("molecular_envelope", envelope)) if m > 0]
    sizes = {"rock_core": core, "ice_mantle": mantle, "molecular_envelope": envelope}
    if fetch.locate(hhe_eos.DATASET).use_stand_in(STAND_IN_POLYTROPE):
        # Constant-density closure scaled so the layers fill the radius: the densities are the stand-in ratios times one scale.
        ratios = ICE_GIANT_RELATIVE_DENSITY
        volume = 4.0 / 3.0 * math.pi * (subject.radius_km * 1e3) ** 3
        scale = volume / math.fsum(m * kg / ratios[kind] for m, kind in
                                   ((core, "rock_core"), (mantle, "ice_mantle"), (envelope, "molecular_envelope")))
        profile = _ConstantDensityProfile([(sizes[kind] * kg, ratios[kind] / scale) for kind in names])
    else:
        envelope_ice = math.fsum(value for name, value in envelope_take.items() if name in bin_.ices)
        profile = _LayeredProfile(giant_structure.fit_adiabat(giant_structure.StructureSubject(
            mass * kg, subject.radius_km * 1e3, subject.helium_share, envelope_z=z_out,
            layers=tuple(giant_structure.Layer(kind, sizes[name] * kg) for name, kind in
                         (("rock_core", "core"), ("ice_mantle", "mantle"), ("molecular_envelope", "envelope")) if sizes[name] > 0),
            envelope_ice_fraction=min(envelope_ice / math.fsum(envelope_take.values()), 1.0) if envelope_take else 0.0,
            iron_share=subject.iron_share_of_refractory, temperature_1bar_k=subject.temperature_1bar_k)),
            [sizes[kind] * kg for kind in names])
    index = {kind: i for i, kind in enumerate(names)}

    layers = []
    env_top, env_bottom = profile.boundary_pressures(index["molecular_envelope"])
    y = subject.helium_share

    def envelope_part(layer_mass):
        heavy_part = {name: value * layer_mass / envelope for name, value in envelope_take.items()}
        hhe = layer_mass - math.fsum(heavy_part.values())
        return heavy_part, hhe * y, hhe * (1.0 - y)

    outer = profile.outer_radii
    radius_1bar = outer[-1]

    def inner_of(kind):
        return outer[index[kind] - 1] / radius_1bar if index[kind] > 0 else 0.0

    radius_mm = profile.radius_at_pressure(p_mm_gpa * 1e9)
    envelope_inner = inner_of("molecular_envelope")
    if radius_mm is not None and radius_mm > envelope_inner * radius_1bar:
        upper = (profile.inner_mass[index["molecular_envelope"]] - profile.mass_inside(radius_mm)) / kg
        mm_fraction = radius_mm / radius_1bar
        parts = [("molecular_envelope", "molecular_fluid", upper, env_top, p_mm_gpa * 1e9, 1.0, mm_fraction),
                 ("metallic_envelope", "metallic_fluid", envelope - upper, p_mm_gpa * 1e9, env_bottom, mm_fraction, envelope_inner)]
    elif env_top >= p_mm_gpa * 1e9:
        parts = [("metallic_envelope", "metallic_fluid", envelope, env_top, env_bottom, 1.0, envelope_inner)]
    else:
        parts = [("molecular_envelope", "molecular_fluid", envelope, env_top, env_bottom, 1.0, envelope_inner)]
    for kind, state, layer_mass, top, bottom, top_radius, bottom_radius in parts:
        heavy_part, helium, hydrogen = envelope_part(layer_mass)
        layers.append(_layer(kind, state, top, bottom, heavy_part, helium, hydrogen, top_radius, bottom_radius))
    if mantle > 0:
        top, bottom = profile.boundary_pressures(index["ice_mantle"])
        layers.append(_layer("ice_mantle", "superionic" if bottom / 1e9 >= SUPERIONIC_PRESSURE_GPA else "ionic_fluid",
                             top, bottom, mantle_take, 0.0, 0.0, outer[index["ice_mantle"]] / radius_1bar, inner_of("ice_mantle")))
    if core > 0:
        top, bottom = profile.boundary_pressures(index["rock_core"])
        layers.append(_layer("rock_core", "molten", top, bottom, core_take, 0.0, 0.0, outer[index["rock_core"]] / radius_1bar, 0.0))
    return GiantInterior(tuple(_fix_water(layers)), mass, "ice_giant", p_mm_gpa, 1.0)


def _brown_dwarf(subject):
    mass = subject.mass_earth
    profile = _structure(subject, 1.5)
    p_mm_gpa = _log_uniform(*METALLIZATION_PRESSURE_RANGE_GPA)
    ratio_mm = profile.at_pressure(p_mm_gpa * 1e9)[1]
    mass_mol = mass * (1.0 - ratio_mm)
    heavy = subject.heavy_mass_earth
    bin_ = _Bin(subject, refractory_share=1.0)
    z = heavy / mass
    layers = []
    y = subject.helium_share
    for kind, layer_mass, state, top, bottom in (
            ("molecular_envelope", mass_mol, "molecular_fluid", 1e5, p_mm_gpa * 1e9),
            ("metallic_envelope", mass - mass_mol, "metallic_fluid", p_mm_gpa * 1e9, profile.central_pressure_pa)):
        if layer_mass <= 0:
            continue
        heavy_part = bin_.take(z * layer_mass, "refractory")
        hhe = layer_mass - math.fsum(heavy_part.values())
        layers.append(_layer(kind, state, top, bottom, heavy_part, hhe * y, hhe * (1.0 - y),
                             1.0 if kind == "molecular_envelope" else profile.at_pressure(p_mm_gpa * 1e9)[0],
                             profile.at_pressure(p_mm_gpa * 1e9)[0] if kind == "molecular_envelope" else 0.0))
    rest = bin_.take(bin_.total, "refractory")
    if rest:  # rounding dust
        last = layers[-1]
        components = dict(last.components)
        for name, value in rest.items():
            components[name] = components.get(name, 0.0) + value
        layers[-1] = Layer(last.kind, last.state, last.top_pressure_gpa, last.bottom_pressure_gpa,
                           math.fsum(components.values()), components, last.top_radius_fraction,
                           last.bottom_radius_fraction)
    return GiantInterior(tuple(layers), mass, "brown_dwarf", p_mm_gpa, 1.0)


# ---------------------------------------------------------------------------------------------------------------------
# Description
# ---------------------------------------------------------------------------------------------------------------------

LAYER_NAMES = {
    "molecular_envelope": "molecular hydrogen and helium envelope",
    "metallic_envelope": "metallic hydrogen envelope",
    "helium_gradient": "helium-enriched metallic zone",
    "helium_shell": "helium-rich shell",
    "dilute_core": "dilute core",
    "compact_core": "compact core of rock and iron",
    "ice_mantle": "mantle of water, ammonia and methane",
    "rock_core": "core of rock and iron",
}
"""dict: What each layer kind is called in prose."""

MATERIAL_NAMES = {"hydrogen": "hydrogen", "helium": "helium", "iron": "iron", "rock": "rock", "water_liquid": "water",
                  "ice_x": "superionic water", "nh3_ice": "ammonia", "ch4_ice": "methane"}
"""dict: What each material of the shared table is called in prose."""

STATE_NAMES = {"molecular_fluid": "molecular fluid", "metallic_fluid": "metallic fluid", "ionic_fluid": "ionic fluid",
               "superionic": "superionic", "molten": "molten"}
"""dict: What each state is called in prose."""

EARTH_MASS_JUPITER = 317.83
"""float: Earth masses in one Jupiter mass, for the prose."""


def _mass_text(mass_earth):
    if mass_earth >= 0.5 * EARTH_MASS_JUPITER:
        return f"{mass_earth / EARTH_MASS_JUPITER:.2f} Jupiter masses"
    return f"{mass_earth:.1f} Earth masses" if mass_earth >= 0.1 else f"{mass_earth:.2f} Earth masses"


def _pressure_text(pressure_gpa):
    if pressure_gpa < 0.1:
        return f"{pressure_gpa * 1e4:.0f} bar"
    if pressure_gpa >= 1e3:
        return f"{pressure_gpa / 100:.0f} Mbar" if pressure_gpa >= 1e4 else f"{pressure_gpa / 100:.1f} Mbar"
    return f"{pressure_gpa / 100:.2f} Mbar" if pressure_gpa >= 10 else f"{pressure_gpa:.1f} GPa"


def describe_interior(body, min_share=0.02):
    """
    The layers of `body` (a `GiantInterior`) as prose, one sentence per layer from the top down: its name, state,
    mass, the pressures it spans, and the materials that make up at least `min_share` of it, largest first.

    Returns:
        list of str: The sentences.
    """
    sentences = []
    for layer in body.layers:
        parts = sorted(layer.component_fractions().items(), key=lambda item: -item[1])
        made = [f"{MATERIAL_NAMES.get(name, name)} {share * 100:.0f}%" for name, share in parts if share >= min_share]
        name = LAYER_NAMES[layer.kind]
        sentences.append(
            f"The {name} ({STATE_NAMES[layer.state]}) holds {_mass_text(layer.mass_earth)}, from "
            f"{_pressure_text(layer.top_pressure_gpa)} down to {_pressure_text(layer.bottom_pressure_gpa)}; it is "
            f"{', '.join(made)}.")
    return sentences


INTERIOR_KIND_WORDING = {
    "gas_giant": "Soundings read a gas giant: a deep envelope of hydrogen and helium over a compact core",
    "ice_giant": "Soundings read an ice giant: a thin envelope of hydrogen and helium over a deep mantle of hot water, ammonia and methane",
    "brown_dwarf": "Soundings read a brown dwarf: a body too massive to stay a planet, hydrogen all the way down to a small core",
}
"""dict: The opening clause for each kind of interior (`INTERIOR_KINDS`); a test pins that every kind has one."""

INTERIOR_FEATURE_WORDING = {
    "metallic_envelope": "hydrogen squeezed to a liquid metal below {pressure}, where the planet's magnetic field is made",
    "dilute_core": "a dilute core, rock and ice mixed into the hydrogen instead of sitting in a ball at the centre",
    "helium_shell": "a shell of helium the rain left behind",
    "helium_rain": "an envelope short of helium, which has sunk into the deep",
    "ice_mantle": "a conducting mantle of ionic water where the planet's magnetic field is made",
}
"""dict: The clause for each feature a giant's interior can have, in the voice of the pages. A test pins the keys."""

HELIUM_RAIN_BELOW = 0.97
"""float: A `helium_depletion` under this reads as an envelope the helium has rained out of."""


def interior_wording_keys(body):
    """
    Which wording of `INTERIOR_KIND_WORDING` and `INTERIOR_FEATURE_WORDING` a `GiantInterior` uses, the kind first and then its features in the order
    of the layers (the helium rain last).
    """
    keys = [body.kind]
    keys += [layer.kind for layer in body.layers if layer.kind in INTERIOR_FEATURE_WORDING]
    if body.helium_depletion < HELIUM_RAIN_BELOW:
        keys.append("helium_rain")
    return keys


def summarize_interior(body):
    """
    One sentence that says what a giant's interior is: its kind, then the metallic layer, the dilute core, the helium shell and the helium rain it has.

    Returns:
        str: The sentence.
    """
    keys = interior_wording_keys(body)
    kind, features = keys[0], keys[1:]
    pressure = body.metallization_pressure_gpa
    clauses = [INTERIOR_FEATURE_WORDING[key].format(pressure=_pressure_text(pressure)) for key in features]
    if not clauses:
        return INTERIOR_KIND_WORDING[kind] + "."
    joined = clauses[0] if len(clauses) == 1 else ", ".join(clauses[:-1]) + " and " + clauses[-1]
    return f"{INTERIOR_KIND_WORDING[kind]}, with {joined}."
