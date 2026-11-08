# planetgen/physics/core_composition.py

"""
What a rocky body's metal core is made of (GEN.1364): Fe-Ni with sulfur, silicon, oxygen, carbon and hydrogen, and what that does to the inner core.

`bulk_draw` draws the metal's light elements as one lump (`core_light_fraction`, a share of the metal's mass) and the star's Ni per Fe gives
`core_ni_fraction`. This module splits the lump and gives the core thermal history (`core_heat`) the three things the split decides.

1. **The split** (`split`). How much of each light element goes into the metal is set by how oxidised the silicate it formed against was. The redox is
   `delta_iw`, log10 units below the iron-wustite buffer, from the FeO the rock holds: the reaction `Fe + 1/2 O2 = FeO` puts `fO2` at `(x_FeO / x_Fe)^2`
   times a constant, so `delta_IW = 2 log10(x_FeO / x_Fe)` with `x_FeO` the mole fraction of FeO in the silicate and `x_Fe` of Fe in the metal (`METAL_X_FE`,
   0.85). A reduced body (little FeO) puts silicon and carbon in its core and keeps sulfur and oxygen in the rock; an oxidised one puts sulfur, oxygen and
   hydrogen in its core and no silicon. The mean shares are those of three bodies, interpolated on `delta_iw` (`ANCHORS`); the pressure at the base of the
   mantle tilts the oxygen and silicon (`PRESSURE_EXPONENT`); `CONCENTRATION` is the Dirichlet scatter about those means (one stream per body, keyed on the
   body, so no other draw moves).
2. **The liquidus.** The iron liquidus of `core_heat` is Earth's, which is already an alloy's. Light elements the solid does not take (`SOLID_TO_LIQUID`)
   lower the freezing point in proportion to the mole fraction they leave in the liquid, `LIQUIDUS_DEPRESSION_K`, so a core with more of them than Earth's
   freezes later (`liquidus_factor`).
3. **The compositional energy.** Freezing inner core rejects the light elements it does not take. The energy that releases is scaled from Earth's
   (`core_heat.GRAVITATIONAL_ENERGY_J_KG`) by the mass rejected per kilogram frozen (`rejected_per_kg`, `composition_energy_factor`). It is not scaled by the core's size: the
   constant of CC6 never was, and the census (reports/gen1364-core-split-census-2026-11-08.md) shows that scaling by `G M / R` would stop most small bodies' dynamos.

A body from before GEN.1364 stored no split: its core takes the Earth's reference mix (`EARTH`), which leaves its liquidus and its energy as they were.

Every number that is not a published value is a stand-in in the gap registry (docs/design/research/gap-filling-plan.md section 6.213).
"""

import math
from dataclasses import dataclass

from planetgen.util import draw

ELEMENTS = ("s", "si", "o", "c", "h")
"""tuple: The light elements, in the order of every list here."""
MOLAR_MASS = {"fe": 55.845, "ni": 58.693, "s": 32.06, "si": 28.086, "o": 15.999, "c": 12.011, "h": 1.008}
"""dict: g/mol (IUPAC standard atomic weights)."""
BODY_FIELDS = ("core_s_fraction", "core_si_fraction", "core_o_fraction", "core_c_fraction", "core_h_fraction")
"""tuple: The attributes (and columns) the split sets, mass fractions of the core in the order of `ELEMENTS`; `core_ni_fraction` and `core_light_fraction`
(their sum) come with the bulk composition."""
FIELDS_OF = dict(zip(ELEMENTS, BODY_FIELDS))
"""dict: Element to attribute."""

METAL_X_FE = 0.85
"""float: STAND-IN CX1: the mole fraction of Fe in the metal the silicate equilibrates with (the rest Ni and the light elements)."""
FEO_MOLAR_MASS = 71.844
"""float: g/mol of FeO."""
MGO_MOLAR_MASS, SIO2_MOLAR_MASS, AL2O3_MOLAR_MASS = 40.304, 60.083, 101.961
"""float: g/mol of the other major oxides of the rock."""
PRESSURE_EXPONENT = 0.3
"""float: STAND-IN CX2: the exponent of the pressure at the base of the mantle over that of the anchor that scales the oxygen and silicon shares (both
dissolve more in iron at high pressure and temperature)."""
PRESSURE_FACTOR_LIMITS = (0.5, 2.0)
"""tuple: STAND-IN CX2: the factor on a share of oxygen or silicon stays inside these."""
CONCENTRATION = 6.0
"""float: STAND-IN CX3: the Dirichlet concentration of the split about its mean shares (the 15 models of Earth's core spread silicon 1.9 to 7.4 and oxygen
0.8 to 6 wt percent about means of 5.0 and 2.7)."""
ALPHA_FLOOR = 0.02
"""float: STAND-IN CX3: the least Dirichlet parameter, so an element the mean excludes still has a chance to appear, rarely and small."""

SOLID_TO_LIQUID = {"s": 0.02, "si": 0.67, "o": 0.02, "c": 0.5, "h": 0.3}
"""dict: STAND-IN CX4: the mass concentration of each light element in the solid inner core over that in the liquid it freezes from. Silicon is 1 over 1.5 (a
liquid/solid partition coefficient of 1.5 +- 0.3, as used by [R5146], read second-hand); oxygen is highly incompatible in solid iron ([R5146], read, no
number: 0.02 is a stand-in); sulfur, carbon and hydrogen are recalled (sulfur stays in the liquid, carbon and hydrogen go in partly)."""
LIQUIDUS_DEPRESSION_K = 3700.0
"""float: STAND-IN CX5: K the freezing point falls at the inner core boundary (5,400 K) for a mole fraction of 1 of light elements left in the liquid;
the Fe-FeS eutectic falls 550 K from iron's 1,811 K for 0.44 of sulfur (1,250 K per mole fraction), scaled with the temperature."""
LIQUIDUS_FACTOR_LIMITS = (0.7, 1.15)
"""tuple: STAND-IN CX5: the liquidus' factor stays inside these (pure iron's liquidus at the inner core boundary is about 15 percent above the alloy's, recalled; the sulfur-rich
end is a stand-in)."""
ICB_TEMPERATURE_K = 5400.0
"""float: The inner core boundary the liquidus is anchored on (`core_heat.LIQUIDUS_ANCHOR`), K."""
ENERGY_FACTOR_LIMITS = (0.0, 4.0)
"""tuple: STAND-IN CX6: the compositional energy is at most this many times the Earth's per kilogram."""


@dataclass(frozen=True)
class Mix:
    """
    A core's metal as mass fractions of its mass (`fe` is what is left).

    Attributes:
        ni, s, si, o, c, h (float): Mass fractions.
    """

    ni: float
    s: float
    si: float
    o: float
    c: float
    h: float

    @property
    def light(self):
        """The light elements' share of the mass."""
        return self.s + self.si + self.o + self.c + self.h

    @property
    def fe(self):
        """The Fe share."""
        return 1.0 - self.ni - self.light

    def _moles(self, scale=None):
        """Moles per 100 g of the metal: `(metal moles, {element: moles})`; `scale` is a factor on each light element's mass (the solid's)."""
        metal = self.fe / MOLAR_MASS["fe"] + self.ni / MOLAR_MASS["ni"]
        light = {name: getattr(self, name) * (1.0 if scale is None else scale[name]) / MOLAR_MASS[name] for name in ELEMENTS}
        return metal, light

    @property
    def liquid_light_mole_fraction(self):
        """The light elements' mole fraction of the liquid metal."""
        metal, light = self._moles()
        total = metal + sum(light.values())
        return sum(light.values()) / total if total > 0.0 else 0.0

    @property
    def solid_light_mole_fraction(self):
        """The light elements' mole fraction of the solid that freezes from this liquid (`SOLID_TO_LIQUID` of each)."""
        metal, light = self._moles(SOLID_TO_LIQUID)
        total = metal + sum(light.values())
        return sum(light.values()) / total if total > 0.0 else 0.0

    @property
    def rejected_per_kg(self):
        """The mass of light elements a kilogram of the freezing solid leaves in the liquid, kg: `sum x (1 - D)`."""
        return sum(getattr(self, name) * (1.0 - SOLID_TO_LIQUID[name]) for name in ELEMENTS)

    @property
    def mole_fraction_left(self):
        """The mole fraction of light elements the liquid keeps over the solid's, the depression's driver."""
        return max(self.liquid_light_mole_fraction - self.solid_light_mole_fraction, 0.0)


_EARTH_CORE_WT_PERCENT = {"ni": 5.06, "s": 1.83, "si": 4.96, "o": 2.66, "c": 0.79, "h": 0.057}
_EARTH_CORE_FE_WT_PERCENT = 82.8
_earth_total = _EARTH_CORE_FE_WT_PERCENT + sum(_EARTH_CORE_WT_PERCENT.values())
EARTH = Mix(**{name: value / _earth_total for name, value in _EARTH_CORE_WT_PERCENT.items()})
"""Mix: Earth's core, the mean of 15 models ([R4608: Tables 1 and 2], read), the seven elements normalised to 1. It is the
reference of the liquidus and of the compositional energy and the mix of a body with no stored split."""
_MARS_WT_PERCENT = {"ni": 7.4, "s": 6.6, "o": 5.2, "h": 0.9}
_MARS_FE_WT_PERCENT = 79.5
_mars_total = _MARS_FE_WT_PERCENT + sum(_MARS_WT_PERCENT.values())
MARS = Mix(ni=7.4 / _mars_total, s=6.6 / _mars_total, si=0.0, o=5.2 / _mars_total, c=0.0, h=0.9 / _mars_total)
"""Mix: Mars' core, the model of [R4833: Table 8] (read; Co, P, Cu and Zn left out), no silicon or carbon."""
MERCURY_SHARES = {"s": 0.08, "si": 0.76, "o": 0.02, "c": 0.12, "h": 0.02}
"""dict: STAND-IN CX7: the shares of Mercury's light elements: silicon the most (5 to 25 wt percent in the favoured reduced models), carbon 1 to 4, sulfur at most
1.5 to 2 ([R5145], ranges read; there is no consensus)."""


def _shares(mix):
    light = mix.light
    return {name: getattr(mix, name) / light for name in ELEMENTS}


EARTH_ROCK = {"feo": 0.081, "al2o3": 0.045, "mg_si": 1.24}
MARS_ROCK = {"feo": 0.18, "al2o3": 0.036, "mg_si": 1.17}
"""dict: STAND-IN CX7: the FeO, Al2O3 (mass fractions) and molar Mg/Si of the bulk silicate Earth and Mars that fix their `delta_iw` here."""


def delta_iw(rock_feo, rock_al2o3=None, rock_mg_si=None):
    """
    The redox the rock sets, log10 units against the iron-wustite buffer: `2 log10(x_FeO / x_Fe)` (see the module text), `x_FeO` the mole fraction of FeO among
    the rock's FeO, MgO, SiO2 and Al2O3 (CaO is left out, a few percent of the moles), `x_Fe` `METAL_X_FE`; a rock with no FeO is at -8.

    Args:
        rock_feo (float): FeO as a mass fraction of the rock.
        rock_al2o3 (float or None): Al2O3, a mass fraction (0.045 where `None`).
        rock_mg_si (float or None): Molar Mg/Si (1.24 where `None`).
    """
    al = 0.045 if rock_al2o3 is None else rock_al2o3
    ratio = 1.24 if rock_mg_si is None else rock_mg_si
    rest = max(1.0 - rock_feo - al, 1e-6)
    silica_mol = rest / (ratio * MGO_MOLAR_MASS + SIO2_MOLAR_MASS)
    moles = {"feo": rock_feo / FEO_MOLAR_MASS, "mgo": ratio * silica_mol, "sio2": silica_mol, "al2o3": al / AL2O3_MOLAR_MASS}
    x_feo = moles["feo"] / sum(moles.values())
    if x_feo <= 0.0:
        return -8.0
    return max(2.0 * math.log10(x_feo / METAL_X_FE), -8.0)


EARTH_DELTA_IW = delta_iw(EARTH_ROCK["feo"], EARTH_ROCK["al2o3"], EARTH_ROCK["mg_si"])
MARS_DELTA_IW = delta_iw(MARS_ROCK["feo"], MARS_ROCK["al2o3"], MARS_ROCK["mg_si"])
MERCURY_DELTA_IW = -5.4
ANCHORS = (
    (MERCURY_DELTA_IW, MERCURY_SHARES, 40.0),
    (EARTH_DELTA_IW, _shares(EARTH), 136.0),
    (MARS_DELTA_IW, _shares(MARS), 20.0),
)
"""tuple: `(delta_iw, shares, pressure_gpa)` of the three bodies the mean split is interpolated between: the core formed against the redox shown (Mercury's
is IW-5.4, read; Earth's and Mars' are the model's own `delta_iw` of their rock) and the pressure at the base of their mantle (STAND-IN CX7)."""


def mean_shares(redox, pressure_gpa):
    """
    `{element: share}` of the light elements' mass, the mean of a core formed against `redox` (`delta_iw`) and a mantle of `pressure_gpa` at its base:
    the anchors' shares interpolated on the redox (held at the ends), the oxygen and silicon scaled by `(P / P_anchor)^PRESSURE_EXPONENT` (the anchor's
    pressure interpolated the same way), then added up to 1.
    """
    low, middle, high = ANCHORS
    if redox <= low[0]:
        a, b, t = low, low, 0.0
    elif redox >= high[0]:
        a, b, t = high, high, 0.0
    else:
        a, b = (low, middle) if redox < middle[0] else (middle, high)
        t = (redox - a[0]) / (b[0] - a[0])
    shares = {name: (1.0 - t) * a[1][name] + t * b[1][name] for name in ELEMENTS}
    reference = (1.0 - t) * a[2] + t * b[2]
    factor = min(max((max(pressure_gpa, 1.0e-3) / reference) ** PRESSURE_EXPONENT, PRESSURE_FACTOR_LIMITS[0]), PRESSURE_FACTOR_LIMITS[1])
    for name in ("o", "si"):
        shares[name] *= factor
    total = sum(shares.values())
    return {name: value / total for name, value in shares.items()}


def split(light, redox, pressure_gpa, key):
    """
    `{element: mass fraction of the metal}` of the light-element total `light`, the Dirichlet scatter about `mean_shares` (one stream keyed on `key`: the same
    body gets the same split and no other draw moves).
    """
    means = mean_shares(redox, pressure_gpa)
    alphas = [max(means[name] * CONCENTRATION, ALPHA_FLOOR) for name in ELEMENTS]
    shares = draw.Stream(f"core-split:{key}").dirichlet(alphas)
    return {name: light * share for name, share in zip(ELEMENTS, shares)}


def mix_of(body):
    """
    The `Mix` stored on `body`; `EARTH` for a body that stores none (a row from before GEN.1364), or a body whose stored metal has no Ni or light elements
    recorded.
    """
    ni = getattr(body, "core_ni_fraction", None)
    parts = {name: getattr(body, FIELDS_OF[name], None) for name in ELEMENTS}
    if ni is None or any(value is None for value in parts.values()):
        return EARTH
    return Mix(ni=float(ni), **{name: float(value) for name, value in parts.items()})


def liquidus_factor(mix):
    """The factor on the model's iron liquidus (Earth's, alloyed) for a core of `mix`: 1 for the Earth's, below 1 for a core with more light elements left
    in its liquid and above 1 for less (`LIQUIDUS_DEPRESSION_K`)."""
    shift_k = LIQUIDUS_DEPRESSION_K * (mix.mole_fraction_left - EARTH.mole_fraction_left)
    return min(max(1.0 - shift_k / ICB_TEMPERATURE_K, LIQUIDUS_FACTOR_LIMITS[0]), LIQUIDUS_FACTOR_LIMITS[1])


def composition_energy_factor(mix):
    """The compositional energy per kilogram of inner core frozen over the Earth's, from the mass of light elements left in the liquid per kilogram frozen
    (`Mix.rejected_per_kg`) against the Earth's."""
    low, high = ENERGY_FACTOR_LIMITS
    return min(max(mix.rejected_per_kg / EARTH.rejected_per_kg, low), high)
