# planetgen/physics/radiogenic.py

"""
One radiogenic module (GEN.459, docs/design/radiation-model.md R15 and section 6, design-drafts/parts/E section 2.8 and 3.3)
===========================================================================================================================

The potassium, thorium, uranium and aluminium-26 a rocky body carries, the heat they give and the activity they leave in its
surface rock. The heat of the interior (`planets.internal_heat_flux_w_m2`, `rogue_surface`) and the ground's dose
(`radiation.ground_dose_msv_yr`) read the same stored concentrations, so a body cannot be rich in thorium for its dose and
average for its heat: there is no independent dose draw (`a_rad` is gone) and no second set of isotope constants (the half-lives,
decay energies and masses are `physics.elements`, GEN.438).

What a body stores (`BODY_FIELDS`, in `planets` and `moons`):

* `c_k0_ppm`, `c_th0_ppm`, `c_u0_ppm`: the element's mass in the bulk silicate at formation, parts per million by mass. At formation
  because the isotopes decay: the amount of each at the body's age is the formation amount times `exp(-lambda t)`, with the isotope
  split of the element at formation (`formation_share`) in place of today's.
* `al26_ratio`: the 26Al/27Al atom ratio at formation.
* `urey_ratio`: the present radiogenic share of the heat flowing out.

The draws (`draw_content`), in a fixed order so none moves another:

1. **Thorium and uranium** follow the host star's Th/Mg: the star's own relation (`stellar_abundances.birth_radiogenic`, the Eu
   relation of GEN.439 and the [Th/Eu] scatter, always the same for one star because the draws come from the star's own stream
   `abundance:{key}`), taken against the solar-scale star and truncated at plus or minus 0.5 dex (R15; verified for the range, partly
   verified for Th away from solar [Fe/H], which is inferred through Eu). A body with no host composition (a rogue planet, a star
   with no stored [Fe/H]) draws the offset from a normal of 0.22 dex truncated at 0.5 dex (R15, the 1 to 99 percent range of the
   Hypatia catalogue in Nimmo et al. 2020; a stand-in for the unknown host). Uranium takes the thorium offset and an extra
   N(0, 0.05) dex (Dauphas 2005 via Wang 2020, partly verified; a placeholder until a source gives the Th/U scatter).
2. **Potassium** is not fitted to the star (it is volatile; Th and U are refractory lithophiles, condensation temperatures in part E
   2.8b): `d_K`, log-uniform 0.03 to 1 of the CI K/U of 67407, times the body's uranium (R15; a placeholder tied to the large K
   differences of Solar System bodies, Nimmo et al. 2020, until the K depletion survey of the research list is read).
3. **Aluminium-26** is log-uniform from 0.1 to 10 times the Solar System ratio 5.23e-5 (Lugaro et al. 2018 Table 2; the range of
   Lichtenberg et al. 2019; a stand-in until GEN.476 derives it from the formation environment).
4. **The Urey ratio** is uniform from 0.1 to 0.9 (part E 3.3; Earth's 0.47 to 0.52 from Nimmo et al. 2020 is a test target, not a
   preset).

Reference abundances are Palme and O'Neill 2003 edition, Table 4 for the bulk silicate Earth (the 2014 values of our files cannot be
confirmed from a paper; Rule applied: verified) and Rudnick and Gao 2014 Table 11 for the upper crust. Solar System bodies are targets
only: the bulk silicate Earth is the reference the star's offset scales, never a preset.

Stand-ins this module labels (until the item that replaces each): the bulk silicate share of a body's mass
(the body's own rock share, `rock_inputs`; `tuning.RADIOGENIC_SILICATE_MASS_FRACTION`, the Earth's 1 minus its core fraction, for a body with none stored), the aluminium share
of that silicate (`BULK_SILICATE_EARTH_AL_MASS_FRACTION`, scaled by the body's Al2O3 in `rock_inputs`), the isotope split at formation (the Solar System's birth mix for every star; an
age law for the r-process contribution is not read), the surface rock (`SURFACE_RESERVOIRS`, a basaltic crust with a drawn felsic
share; GEN.441 gives the ground type), and 235U's gamma series (left out of the dose; 4 percent on Earth, more on a young body).
"""

import functools
import math
import statistics
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import bulk_composition, elements
from planetgen.physics import stellar_abundances as sa
from planetgen.physics import stellar_composition as sc
from planetgen.util import draw

GY_S = sa.GY_S
"""float: Seconds in a gigayear (Julian), the year the half-lives of `elements` are counted in."""

BODY_FIELDS = ("c_k0_ppm", "c_th0_ppm", "c_u0_ppm", "al26_ratio", "urey_ratio")
"""tuple: What `generate_content` sets on a planet or moon, as stored in `planets` and `moons`."""

PARENTS = {"K": ("K-39", "K-40", "K-41"), "Th": ("Th-232",), "U": ("U-238", "U-235")}
"""dict: The natural isotopes of each heat-producing element the module follows. 234U (5.4e-5 of the atoms, a daughter of 238U in
equilibrium) is left out: its share of the mass is below a part in ten thousand."""

HEAT_PARENTS = ("K-40", "Th-232", "U-238", "U-235")
"""tuple: The radioactive parents whose heat is summed; the 26Al term is added separately."""

REFERENCE_AGE_GY = tuning.STELLAR_SOLAR_SYSTEM_AGE_GY
"""float: How long ago the Solar System formed (the age of the oldest solids): the reference tables below are today's values, and
the isotope split at formation is today's decayed back this far."""

# --- Reference abundances, ppm by mass of the element, today (Solar System age) --------------------------------------------

ELEMENTS = ("K", "Th", "U")
"""tuple: The heat-producing elements, in the order of `Reservoir.ppm`."""

BULK_SILICATE_EARTH = (260.0, 0.0834, 0.0218)
"""tuple: `(K, Th, U)` ppm of the bulk silicate Earth, Palme and O'Neill 2003 edition Table 4 (260 ppm, 83.4 ppb, 21.8 ppb; read twice,
radiation-model.md section 6). The Th and U of the 2014 edition in our files (84.9 and 22.9 ppb) cannot be confirmed from a paper."""

CI_CHONDRITE_K_OVER_U = 546.0 / 0.0081
"""float: The mass ratio of potassium to uranium in CI chondrites, 67407 (part E 2.8b, Palme and O'Neill via pyrolite, 2ND-HAND); the
bulk silicate Earth has 11927 with the values above, a depletion of 0.177."""

BULK_SILICATE_EARTH_AL_MASS_FRACTION = 0.0238
"""float: Aluminium in the bulk silicate Earth, 2.38 weight percent (part E 2.8f). The stand-in for a body with no drawn Al2O3
(`rock_inputs` scales it by the body's own)."""

SOLAR_AL26_RATIO = 5.23e-5
"""float: 26Al/27Al at formation of the Solar System, (5.23 +- 0.13)e-5 (Lugaro et al. 2018, arXiv:1808.00233 Table 2; READ)."""

# --- Surface rock: a stand-in for the ground type of GEN.441 ----------------------------------------------------------------

@dataclass(frozen=True)
class Reservoir:
    """
    A rock whose K, Th and U are known today.

    Attributes:
        name (str): What it is.
        ppm (tuple): `(K, Th, U)` ppm by mass (`ELEMENTS`).
        source (str): Where the numbers come from and how sure they are.
    """

    name: str
    ppm: tuple
    source: str

    def over_bulk_silicate(self):
        """`(K, Th, U)` of this rock over the bulk silicate Earth's (`BULK_SILICATE_EARTH`): how many times richer than the
        undifferentiated silicate it is."""
        return tuple(value / bulk for value, bulk in zip(self.ppm, BULK_SILICATE_EARTH))


SURFACE_RESERVOIRS = {
    "basalt": Reservoir("N-MORB basalt", (1162.0, 0.252, 0.083), "Gale et al. 2013 via part E 2.8b (2ND-HAND)"),
    "felsic": Reservoir("upper continental crust", (23244.0, 10.5, 2.7),
                        "Rudnick and Gao 2014 Table 11, the 2003 text (K2O 2.80 weight percent); verified"),
}
"""dict: The two rocks the surface is a mix of until GEN.441 gives the ground type: a basaltic crust and a felsic share. A body with
a stagnant lid and no tectonic reprocessing is all basalt; a mobile lid reprocesses some of it into felsic crust (part E 3.2 item 5)."""

# --- Dose coefficients: air kerma at 1 m over a flat infinite source ---------------------------------------------------

GAMMA_NGY_H_PER_BQ_KG = {"K-40": 0.0417, "U-238": 0.462, "Th-232": 0.604}
"""dict: nGy/h in air at 1 m per Bq/kg of the isotope in the ground; the U-238 and Th-232 values are for the whole series in secular
equilibrium (UNSCEAR 2000, read in part 07 2.10). With today's isotope mix they give 5.70, 2.44 and 13.1 nGy/h per ppm U, ppm Th and
weight percent K (radiation-model.md section 6)."""


class RadiogenicError(ValueError):
    """A concentration or an age this module cannot use."""


@dataclass(frozen=True)
class Content:
    """
    What a body carries (`BODY_FIELDS`).

    Attributes:
        k0_ppm (float): Potassium in the bulk silicate at formation, ppm by mass of the element.
        th0_ppm (float): Thorium, the same.
        u0_ppm (float): Uranium, the same.
        al26_ratio (float): 26Al/27Al at formation.
        urey_ratio (float): Radiogenic heat over the heat flowing out, today.
    """

    k0_ppm: float
    th0_ppm: float
    u0_ppm: float
    al26_ratio: float
    urey_ratio: float

    def ppm(self, symbol):
        """The element's formation concentration, ppm."""
        return {"K": self.k0_ppm, "Th": self.th0_ppm, "U": self.u0_ppm}[symbol]


def content_of(body):
    """The `Content` a body stores, or `None` when it stores none (a giant, or a row from before GEN.459)."""
    values = [getattr(body, name, None) for name in BODY_FIELDS]
    return None if any(value is None for value in values) else Content(*values)


# --- Isotopes at formation and at the body's age ------------------------------------------------------------------------

@functools.cache
def formation_mix(symbol):
    """
    How the element's mass is split among the isotopes in `PARENTS` at formation, and how much heavier the element was then than
    today: `({isotope: share of the element's mass at formation}, formation mass of the element / its mass today)`. The mix today
    is the natural abundances; each radioactive isotope is decayed back `REFERENCE_AGE_GY` (the Solar System's birth mix, a stand-in
    for every star).
    """
    entries = []
    for name in PARENTS[symbol]:
        found = elements.isotope(name)
        now = found.abundance * found.mass_u
        then = now * math.exp(found.decay_constant * REFERENCE_AGE_GY * GY_S) if not found.stable else now
        entries.append((name, now, then))
    total_now = math.fsum(now for _name, now, _then in entries)
    total_then = math.fsum(then for _name, _now, then in entries)
    return {name: then / total_then for name, _now, then in entries}, total_then / total_now


def isotope_mass_fractions(content, age_gy):
    """`{isotope: kg per kg of the bulk silicate}` of the parents at `age_gy` after formation: the element's formation concentration
    times the isotope's share of it, decayed. 26Al is not in it (`aluminium_26_mass_fraction`)."""
    if age_gy < 0:
        raise RadiogenicError(f"an age of {age_gy} Gy is before formation")
    found = {}
    for symbol in ELEMENTS:
        mix, _factor = formation_mix(symbol)
        for name in PARENTS[symbol]:
            isotope = elements.isotope(name)
            decay = 1.0 if isotope.stable else math.exp(-isotope.decay_constant * age_gy * GY_S)
            found[name] = content.ppm(symbol) * 1e-6 * mix[name] * decay
    return found


def aluminium_26_mass_fraction(content, age_gy, aluminium=BULK_SILICATE_EARTH_AL_MASS_FRACTION):
    """kg of 26Al per kg of the bulk silicate at `age_gy`: the aluminium share, the atom ratio at formation (times the mass ratio of
    the two isotopes), decayed with the 0.72 Myr half-life."""
    found = elements.isotope("Al-26")
    mass_ratio = found.mass_u / elements.isotope("Al-27").mass_u
    return aluminium * content.al26_ratio * mass_ratio * math.exp(-found.decay_constant * age_gy * GY_S)


def heat_w_per_kg(content, age_gy, aluminium=BULK_SILICATE_EARTH_AL_MASS_FRACTION):
    """The heat of a kilogram of the bulk silicate at `age_gy`, W/kg: the sum over the four long-lived parents (`elements.
    heat_per_kg_of_isotope`, Dye 2012 equation 6, with the whole chain's energy for U and Th) and 26Al."""
    masses = isotope_mass_fractions(content, age_gy)
    heat = math.fsum(masses[name] * elements.heat_per_kg_of_isotope(name) for name in HEAT_PARENTS)
    return heat + aluminium_26_mass_fraction(content, age_gy, aluminium) * elements.heat_per_kg_of_isotope("Al-26")


def heat_w(content, silicate_mass_kg, age_gy, aluminium=BULK_SILICATE_EARTH_AL_MASS_FRACTION):
    """The radiogenic heat of `silicate_mass_kg` of bulk silicate at `age_gy`, W."""
    return silicate_mass_kg * heat_w_per_kg(content, age_gy, aluminium)


def radiogenic_flux_w_m2(content, mass_kg, radius_km, age_gy, silicate_fraction=None, aluminium=BULK_SILICATE_EARTH_AL_MASS_FRACTION):
    """The radiogenic heat reaching the surface, W/m2: the bulk silicate's heat (`silicate_fraction` of the mass, the stand-in
    `tuning.RADIOGENIC_SILICATE_MASS_FRACTION` for a body with no stored composition) over 4 pi R^2."""
    if silicate_fraction is None:
        silicate_fraction = tuning.RADIOGENIC_SILICATE_MASS_FRACTION
    area_m2 = 4.0 * math.pi * (radius_km * 1000.0) ** 2
    return heat_w(content, mass_kg * silicate_fraction, age_gy, aluminium) / area_m2


def heat_flow_w_m2(content, mass_kg, radius_km, age_gy, silicate_fraction=None, aluminium=BULK_SILICATE_EARTH_AL_MASS_FRACTION):
    """The total heat flowing out, W/m2: the radiogenic flux over the Urey ratio, the rest being the cooling of the interior."""
    return radiogenic_flux_w_m2(content, mass_kg, radius_km, age_gy, silicate_fraction, aluminium) / content.urey_ratio


BULK_SILICATE_EARTH_AL2O3_FRACTION = 0.045
"""float: Al2O3 in the bulk silicate Earth, 4.5 weight percent (McDonough and Sun 1995): the rock whose aluminium is
`BULK_SILICATE_EARTH_AL_MASS_FRACTION`."""


def rock_inputs(body):
    """`(silicate_fraction, aluminium)` of `body`'s own rock for `heat_flow_w_m2` (GEN.428): the rock's share of its mass and its
    aluminium, the bulk silicate Earth's scaled by the body's Al2O3 over Earth's. `(None, the Earth's)` (the stand-ins) for a body with
    no stored composition or whose Al2O3 is not drawn."""
    fractions = bulk_composition.bulk_fractions(body)
    if fractions is None:
        return None, BULK_SILICATE_EARTH_AL_MASS_FRACTION
    aluminium = BULK_SILICATE_EARTH_AL_MASS_FRACTION
    if fractions.rock_al2o3 is not None:
        aluminium *= fractions.rock_al2o3 / BULK_SILICATE_EARTH_AL2O3_FRACTION
    return fractions.rock, aluminium


# --- The draws --------------------------------------------------------------------------------------------------------------

_NORMAL = statistics.NormalDist()


def truncated_offset_dex(u, sigma=None, limit=None):
    """A normal of `sigma` dex truncated at plus or minus `limit` dex from one uniform `u` in [0, 1) (inverse transform, so one
    draw whatever the value). Defaults: `tuning.RADIOGENIC_OFFSET_SIGMA_DEX` and `RADIOGENIC_OFFSET_LIMIT_DEX`."""
    sigma = tuning.RADIOGENIC_OFFSET_SIGMA_DEX if sigma is None else sigma
    limit = tuning.RADIOGENIC_OFFSET_LIMIT_DEX if limit is None else limit
    low, high = _NORMAL.cdf(-limit / sigma), _NORMAL.cdf(limit / sigma)
    return sigma * _NORMAL.inv_cdf(low + min(max(u, 1e-12), 1.0 - 1e-12) * (high - low))


@functools.lru_cache(maxsize=65536)
def host_offset_dex(fe_h, age_gy, key):
    """
    The thorium and uranium offset of a star's rock from the solar scale, dex: log10 of the star's birth Th/Mg over the Th/Mg of the
    solar-scale star (`stellar_abundances.birth_radiogenic`, the star's draws from `abundance:{key}`), truncated at plus or minus
    `tuning.RADIOGENIC_OFFSET_LIMIT_DEX`. Th and U are r-process products that scatter from star to star; Mg is the rock former they
    are reckoned against (Nimmo et al. 2020).

    Args:
        fe_h (float): The star's [Fe/H].
        age_gy (float): Its age, Gy (the branch of the thin and thick disc depends on it).
        key (str): `stellar_composition.abundance_key(...)`: the same star always has the same offset.
    """
    draws = sc.key_draws(key)
    branch = sa.branch_from_uniform(fe_h, age_gy, draws.u_branch)
    offsets = sa.abundance_offsets(fe_h, branch, draws.normals, draws.u_enhanced, draws.u_carbon)
    star = sa.composition(fe_h, offsets)
    found = sa.birth_radiogenic(star, branch, draws.eu_normal, draws.th_eu_normal)
    th_over_mg = found["mass_fractions"]["Th"] / star.elements["Mg"]
    return max(-tuning.RADIOGENIC_OFFSET_LIMIT_DEX, min(tuning.RADIOGENIC_OFFSET_LIMIT_DEX,
                                                         math.log10(th_over_mg / _solar_th_over_mg())))


@functools.cache
def _solar_th_over_mg():
    """Th/Mg by mass of the solar-scale star: [Fe/H] 0 on the thin branch with every residual and scatter at its mean."""
    star = sa.composition(0.0, sa.abundance_offsets(0.0, sa.BRANCH_THIN, (0.0, 0.0, 0.0, 0.0)))
    found = sa.birth_radiogenic(star, sa.BRANCH_THIN)
    return found["mass_fractions"]["Th"] / star.elements["Mg"]


def star_offset_dex(star):
    """`host_offset_dex` for a star object (`generation.star.Star`, a binary's proxy), or `None` when it has no [Fe/H] or no age to
    read it from. The key is the one the star's [Fe/H] stream uses (its mass, age and distance from the centre), so every planet of
    a star gets the same offset."""
    fe_h = getattr(star, "metallicity_fe_h", None)
    age = getattr(star, "age", None)
    if fe_h is None or age is None:
        return None
    mass, distance = getattr(star, "mass", None), getattr(star, "galactic_center_dist_ly", None)
    # a star made with no place in the galaxy (a single system) has no distance; its [Fe/H] stream is keyed on `None` there
    key = f"{mass!r}:{age!r}:{distance!r}" if distance is None else sc.abundance_key(mass=mass, age_gy=age, distance_ly=distance)
    return host_offset_dex(float(fe_h), float(age), key)


def draw_content(offset_dex=None, rng=draw, refractory_enrichment=1.0):
    """
    A body's `Content`. Every draw is made whatever the branch, in this order: the host-less offset, the uranium extra, the
    potassium depletion, the aluminium-26 factor and the Urey ratio, so the stream after it is the same for every body.

    Args:
        offset_dex (float or None): The host star's Th/U offset (`star_offset_dex`); `None` takes the host-less draw.
        rng: The draw source (`planetgen.util.draw` by default).
        refractory_enrichment (float): How many times the Earth's share the body's refractory elements (Al, Th, U) are in its
            silicate; the body's Al2O3 over the bulk silicate Earth's (part E 2.8b: Th and U scale with it; `rock_inputs`), 1 for a body with none.

    Returns:
        Content: At formation, in the bulk silicate.
    """
    host_less = truncated_offset_dex(rng.random())
    uranium_extra = rng.gauss(0.0, tuning.RADIOGENIC_URANIUM_EXTRA_SIGMA_DEX)
    k_depletion = math.exp(rng.uniform(math.log(tuning.RADIOGENIC_K_DEPLETION_RANGE[0]), math.log(tuning.RADIOGENIC_K_DEPLETION_RANGE[1])))
    al26_factor = math.exp(rng.uniform(math.log(tuning.RADIOGENIC_AL26_FACTOR_RANGE[0]), math.log(tuning.RADIOGENIC_AL26_FACTOR_RANGE[1])))
    urey = rng.uniform(*tuning.RADIOGENIC_UREY_RANGE)
    offset = host_less if offset_dex is None else offset_dex
    _bulk_k, bulk_th, bulk_u = BULK_SILICATE_EARTH
    th_now = bulk_th * 10.0 ** offset * refractory_enrichment
    u_now = bulk_u * 10.0 ** (offset + uranium_extra) * refractory_enrichment
    k_now = k_depletion * CI_CHONDRITE_K_OVER_U * u_now
    return Content(
        k0_ppm=k_now * formation_mix("K")[1], th0_ppm=th_now * formation_mix("Th")[1], u0_ppm=u_now * formation_mix("U")[1],
        al26_ratio=SOLAR_AL26_RATIO * al26_factor, urey_ratio=urey)


def generate_content(body, refractory_enrichment=1.0):
    """Sets `BODY_FIELDS` on a rocky planet or moon from its host star's composition (`star_offset_dex`) and its own draws; a gas
    giant has none (every field `None`)."""
    if body.body_type == "g":
        for name in BODY_FIELDS:
            setattr(body, name, None)
        return
    content = draw_content(star_offset_dex(body.star), refractory_enrichment=refractory_enrichment)
    body.c_k0_ppm, body.c_th0_ppm, body.c_u0_ppm = content.k0_ppm, content.th0_ppm, content.u0_ppm
    body.al26_ratio, body.urey_ratio = content.al26_ratio, content.urey_ratio


# --- The surface rock and its activity --------------------------------------------------------------------------------------

def surface_content(content, felsic_fraction):
    """The bulk silicate's `content` with its K, Th and U raised to the surface rock of `SURFACE_RESERVOIRS`: a basaltic crust with
    `felsic_fraction` of it felsic, each element as many times the bulk silicate's as the rock is over the bulk silicate Earth
    (`Reservoir.over_bulk_silicate`). A stand-in for the ground type of GEN.441. The formation concentrations scale together, so
    every isotope is as rich as its element."""
    if not 0.0 <= felsic_fraction <= 1.0:
        raise RadiogenicError(f"a felsic fraction of {felsic_fraction} is outside 0 to 1")
    basalt = SURFACE_RESERVOIRS["basalt"].over_bulk_silicate()
    felsic = SURFACE_RESERVOIRS["felsic"].over_bulk_silicate()
    k, th, u = ((1.0 - felsic_fraction) * low + felsic_fraction * high for low, high in zip(basalt, felsic))
    return Content(content.k0_ppm * k, content.th0_ppm * th, content.u0_ppm * u, content.al26_ratio, content.urey_ratio)


def activity_bq_per_kg(content, age_gy):
    """`{isotope: Bq/kg of the rock}` of the parents the dose reads (K-40, U-238, Th-232) at `age_gy`: the isotope's mass fraction
    times its specific activity (`elements.specific_activity`). The 226Ra that makes radon is in equilibrium with the 238U."""
    masses = isotope_mass_fractions(content, age_gy)
    return {name: masses[name] * elements.specific_activity(name) for name in GAMMA_NGY_H_PER_BQ_KG}


def gamma_ngy_per_hour(content, age_gy):
    """The dose rate in air at 1 m over `content` rock, nGy/h (infinite flat source; `GAMMA_NGY_H_PER_BQ_KG`)."""
    return math.fsum(GAMMA_NGY_H_PER_BQ_KG[name] * value for name, value in activity_bq_per_kg(content, age_gy).items())


def present_ppm(content, age_gy):
    """`{"K": ppm, "Th": ppm, "U": ppm}` of the element in `content` rock at `age_gy`: the sum of the element's isotopes then."""
    masses = isotope_mass_fractions(content, age_gy)
    return {symbol: math.fsum(masses[name] for name in PARENTS[symbol]) * 1e6 for symbol in ELEMENTS}
