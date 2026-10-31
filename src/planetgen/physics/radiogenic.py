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

The draws (`draw_content`), in a fixed order so none moves another:

1. **Thorium and uranium** follow the host star's Th/Mg: the star's own relation (`stellar_abundances.birth_radiogenic`, the Eu
   relation of GEN.439 and the [Th/Eu] scatter, always the same for one star because the draws come from the star's own stream
   `abundance:{key}`), taken against the solar-scale star and truncated at plus or minus 0.5 dex (R15; verified for the range, partly
   verified for Th away from solar [Fe/H], which is inferred through Eu). A body with no host composition (a rogue planet, a star
   with no stored [Fe/H]) draws the offset from a normal of 0.22 dex truncated at 0.5 dex (R15, the 1 to 99 percent range of the
   Hypatia catalogue in Nimmo et al. 2020; a stand-in for the unknown host). Uranium is then worked out from the
   thorium by the r-process (GEN.759, `rprocess_enrichment`): the star's Th/U at its birth is the production ratio of the r-process
   (one draw from 1.667 to 2.0, Frebel 2007) times the factor `F` the decay of the galaxy's earlier ejecta gives (the disc's star formation
   history, GEN.469, and the host's age; a body with no host takes the typical star of the galaxy), and the ratio today is that
   decayed to the Solar System's age. The old uranium scatter of 0.05 dex is gone.
2. **Potassium** is not fitted to the star (it is volatile; Th and U are refractory lithophiles, condensation temperatures in part E
   2.8b): `d_K`, log-uniform 0.03 to 1 of the CI K/U of 67407, times the body's uranium (R15; a placeholder tied to the large K
   differences of Solar System bodies, Nimmo et al. 2020, until the K depletion survey of the research list is read).
3. **Aluminium-26** is the star's birth ratio (GEN.476, `birth_radioisotopes`): the steady-state mean of its galaxy's gas (the yield of
   massive stars from the gamma-ray mass, over the disc's gas depletion time and the star's aluminium) as a sum of discrete decayed
   additions, so most stars are born with far less than the mean and a few with a great deal. The Solar System's 5.23e-5 (Lugaro et al.
   2018 Table 2) is a test point.

The Urey ratio is not drawn: it is the radiogenic share of the heat an interior gives, an output of `interior_heat` (GEN.699).

Reference abundances are Palme and O'Neill 2003 edition, Table 4 for the bulk silicate Earth (the 2014 values of our files cannot be
confirmed from a paper; Rule applied: verified) and Rudnick and Gao 2014 Table 11 for the upper crust. Solar System bodies are targets
only: the bulk silicate Earth is the reference the star's offset scales, never a preset.

Stand-ins this module labels (until the item that replaces each): the bulk silicate share of a body's mass
(the body's own rock share, `rock_of`; `tuning.RADIOGENIC_SILICATE_MASS_FRACTION`, the Earth's 1 minus its core fraction, for a body with none stored), the aluminium share
of that silicate (`BULK_SILICATE_EARTH_AL_MASS_FRACTION`, scaled by the body's Al2O3 in `rock_of`), the isotope split at formation (the Solar System's birth mix for every star; an
age law for the r-process contribution is not read), and 235U's gamma series (left out of the dose; 4 percent on Earth, more on a young body).
"""

import functools
import math
import statistics
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import birth_radioisotopes, bulk_composition, elements, rprocess_enrichment
from planetgen.physics import stellar_abundances as sa
from planetgen.physics import stellar_composition as sc
from planetgen.util import draw

GY_S = sa.GY_S
"""float: Seconds in a gigayear (Julian), the year the half-lives of `elements` are counted in."""

BODY_FIELDS = ("c_k0_ppm", "c_th0_ppm", "c_u0_ppm", "al26_ratio")
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
"""tuple: The heat-producing elements, in the order of `BULK_SILICATE_EARTH`."""

BULK_SILICATE_EARTH = (260.0, 0.0834, 0.0218)
"""tuple: `(K, Th, U)` ppm of the bulk silicate Earth, Palme and O'Neill 2003 edition Table 4 (260 ppm, 83.4 ppb, 21.8 ppb; read twice,
radiation-model.md section 6). The Th and U of the 2014 edition in our files (84.9 and 22.9 ppb) cannot be confirmed from a paper."""

CI_CHONDRITE_K_OVER_U = 546.0 / 0.0081
"""float: The mass ratio of potassium to uranium in CI chondrites, 67407 (part E 2.8b, Palme and O'Neill via pyrolite, 2ND-HAND); the
bulk silicate Earth has 11927 with the values above, a depletion of 0.177."""

BULK_SILICATE_EARTH_AL_MASS_FRACTION = 0.0238
"""float: Aluminium in the bulk silicate Earth, 2.38 weight percent (part E 2.8f). The stand-in for a body with no drawn Al2O3
(`rock_of` scales it by the body's own)."""

SOLAR_AL26_RATIO = 5.23e-5
"""float: 26Al/27Al at formation of the Solar System, (5.23 +- 0.13)e-5 (Lugaro et al. 2018, arXiv:1808.00233 Table 2; READ). A test point and a
reference for the tests: the generator does not draw from it."""

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
    """

    k0_ppm: float
    th0_ppm: float
    u0_ppm: float
    al26_ratio: float

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
    times the isotope's share of it, decayed. 26Al is not in it (its formation ratio is `Content.al26_ratio`)."""
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


def heat_terms(content, aluminium=BULK_SILICATE_EARTH_AL_MASS_FRACTION):
    """The heat of a kilogram of the bulk silicate as a sum of decays, `((W/kg at formation, decay constant per Gy), ...)`: the four
    long-lived parents (`elements.heat_per_kg_of_isotope`, Dye 2012 equation 6, with the whole chain's energy for U and Th) and 26Al.
    `heat_w_per_kg` is `sum(a * exp(-lambda t))` over them, and the thermal history (`interior_heat`) integrates the same terms."""
    terms = []
    for name in HEAT_PARENTS:
        isotope = elements.isotope(name)
        symbol = name.split("-")[0]
        at_formation = content.ppm(symbol) * 1e-6 * formation_mix(symbol)[0][name]
        terms.append((at_formation * elements.heat_per_kg_of_isotope(name), isotope.decay_constant * GY_S))
    al26 = elements.isotope("Al-26")
    mass_ratio = al26.mass_u / elements.isotope("Al-27").mass_u
    terms.append((aluminium * content.al26_ratio * mass_ratio * elements.heat_per_kg_of_isotope("Al-26"), al26.decay_constant * GY_S))
    return tuple(terms)


def heat_w_per_kg(content, age_gy, aluminium=BULK_SILICATE_EARTH_AL_MASS_FRACTION):
    """The heat of a kilogram of the bulk silicate at `age_gy`, W/kg: the sum of `heat_terms` decayed to then."""
    if age_gy < 0:
        raise RadiogenicError(f"an age of {age_gy} Gy is before formation")
    return math.fsum(amplitude * math.exp(-rate * age_gy) for amplitude, rate in heat_terms(content, aluminium))


def heat_w(content, silicate_mass_kg, age_gy, aluminium=BULK_SILICATE_EARTH_AL_MASS_FRACTION):
    """The radiogenic heat of `silicate_mass_kg` of bulk silicate at `age_gy`, W."""
    return silicate_mass_kg * heat_w_per_kg(content, age_gy, aluminium)


BULK_SILICATE_EARTH_AL2O3_FRACTION = 0.045
"""float: Al2O3 in the bulk silicate Earth, 4.5 weight percent (McDonough and Sun 1995): the rock whose aluminium is
`BULK_SILICATE_EARTH_AL_MASS_FRACTION`."""


@dataclass(frozen=True)
class Rock:
    """
    The parts of a body the heat flows through (`rock_of`).

    Attributes:
        silicate (float): The bulk silicate's share of the body's mass: the heat producer, and the convecting mantle.
        metal (float): The core's share of the body's mass.
        wet (float): The share that is water, ice or envelope (the rock's own radius is smaller than the body's by it).
        aluminium (float): Aluminium in the bulk silicate, mass fraction (26Al's parent).
    """

    silicate: float
    metal: float
    wet: float
    aluminium: float


def rock_of(body):
    """The `Rock` of `body` (GEN.428): its stored bulk composition, with the bulk silicate Earth's aluminium scaled by the body's
    Al2O3 over the Earth's. The stand-in (`tuning.RADIOGENIC_SILICATE_MASS_FRACTION` of silicate, the rest metal, no water, the Earth's
    aluminium) for a body with no stored composition; the Earth's aluminium when its Al2O3 is not drawn."""
    fractions = bulk_composition.bulk_fractions(body)
    if fractions is None:
        silicate = tuning.RADIOGENIC_SILICATE_MASS_FRACTION
        return Rock(silicate, 1.0 - silicate, 0.0, BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    aluminium = BULK_SILICATE_EARTH_AL_MASS_FRACTION
    if fractions.rock_al2o3 is not None:
        aluminium *= fractions.rock_al2o3 / BULK_SILICATE_EARTH_AL2O3_FRACTION
    return Rock(fractions.rock, fractions.metal, fractions.water + fractions.ice + fractions.envelope, aluminium)


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
    branch = sa.branch_of(fe_h, age_gy)
    offsets = sa.abundance_offsets(fe_h, branch, draws.normals, draws.u_enhanced, draws.u_carbon)
    star = sa.composition(fe_h, offsets)
    found = sa.birth_radiogenic(star, branch, draws.eu_normal, draws.th_eu_normal)
    th_over_mg = found["mass_fractions"]["Th"] / star.elements["Mg"]
    return max(-tuning.RADIOGENIC_OFFSET_LIMIT_DEX, min(tuning.RADIOGENIC_OFFSET_LIMIT_DEX,
                                                         math.log10(th_over_mg / _solar_th_over_mg())))


@functools.cache
def _solar_th_over_mg():
    """Th/Mg by mass of the solar-scale star: [Fe/H] 0 on the thin branch with every residual and scatter at its mean."""
    star = sa.composition(0.0, sa.abundance_offsets(0.0, sa.BRANCH_THIN, (0.0,) * len(sa.RESIDUAL_ELEMENTS)))
    found = sa.birth_radiogenic(star, sa.BRANCH_THIN)
    return found["mass_fractions"]["Th"] / star.elements["Mg"]


def _star_key(star):
    """The key a star's streams are named by (its mass and distance from the centre; `stellar_composition.abundance_key`), or `None` when it has no [Fe/H] or no age."""
    fe_h = getattr(star, "metallicity_fe_h", None)
    age = getattr(star, "age", None)
    if fe_h is None or age is None:
        return None
    mass, distance = getattr(star, "mass", None), getattr(star, "galactic_center_dist_ly", None)
    # a star made with no place in the galaxy (a single system) has no distance; its [Fe/H] stream is keyed on `None` there
    return f"{mass!r}:{distance!r}" if distance is None else sc.abundance_key(mass=mass, distance_ly=distance)


def star_offset_dex(star):
    """`host_offset_dex` for a star object (`generation.star.Star`, a binary's proxy), or `None` when it has no [Fe/H] or no age to
    read it from. The key is the star's (`stellar_composition.abundance_key`: its mass and distance from the centre), so every planet
    of a star gets the same offset."""
    key = _star_key(star)
    return None if key is None else host_offset_dex(float(star.metallicity_fe_h), float(star.age), key)


def star_al26_ratio(star, shape=None):
    """
    The 26Al/27Al a star was born with (GEN.476): one sum of decayed additions on a stream named by the star, about the mean of its galaxy's gas (`al26_mean` at its [Fe/H]). Every
    body of a star shares it, being born from the same gas. `None` when the star has no [Fe/H] or no age (a body with no host draws its own).
    """
    key = _star_key(star)
    if key is None:
        return None
    return birth_radioisotopes.sample_ratio("Al26", al26_mean(shape, float(star.metallicity_fe_h)), draw.Stream(f"al26-birth:{key}"))


def th_u_ratio_now(production_ratio, factor):
    """The mass ratio of thorium to uranium in rock of the Solar System's age (the element ratios `BULK_SILICATE_EARTH` is given in) for a
    host whose gas had the r-process `production_ratio` (thorium-232 to uranium-238 by number) and the enrichment `factor` `F`:
    `m(Th-232)/m(U-238) * PR * F * exp[(l238 - l232) t]`, over the share of uranium's mass that is 238 today."""
    thorium, uranium = elements.isotope("Th-232"), elements.isotope("U-238")
    decay = math.exp((uranium.decay_constant - thorium.decay_constant) * REFERENCE_AGE_GY * GY_S)
    share_238 = uranium.abundance * uranium.mass_u / math.fsum(elements.isotope(name).abundance * elements.isotope(name).mass_u
                                                               for name in PARENTS["U"])
    return thorium.mass_u / uranium.mass_u * production_ratio * factor * decay / share_238


def default_th_u_factor():
    """`F` of the typical star of the galaxy a stand-alone command or a test takes (`galaxy.shape_draw.default_shape`)."""
    from planetgen.galaxy import shape_draw
    return rprocess_enrichment.typical_th_u_factor(shape_draw.default_shape(None))


def star_th_u_factor(body):
    """`F` for the host of `body`: its age in the disc of its galaxy (`system_config.GALAXY_SHAPE`, set by the sector that made it; the default
    galaxy's when there is none), or `None` when the host has no age (the typical star of the default galaxy)."""
    from planetgen.galaxy import shape_draw
    age = getattr(getattr(body, "star", None), "age", None)
    if age is None:
        return None
    shape = getattr(getattr(body, "system_config", None), "GALAXY_SHAPE", None) or shape_draw.default_shape(None)
    return rprocess_enrichment.th_u_factor_of_shape(shape, round(float(age), 6))


@functools.lru_cache(maxsize=64)
def _depletion_time_yr(shape):
    """The gas depletion time of the thin disc of a galaxy, years: its gas over its star formation, the rings weighted by their areas as `rprocess_enrichment` weights them
    (the present value stands for every star's birth: the gas at an earlier time is not kept)."""
    from planetgen.galaxy import disk_history
    history = disk_history.history_of(shape)
    gas = float(sum(w * g for w, g in zip(history.radius_scale_lengths, history.gas_now)))
    rate = float(sum(w * r for w, r in zip(history.radius_scale_lengths, history.rate[:, -1])))
    return 1.0e9 * gas / rate


def al26_mean(shape=None, fe_h=0.0):
    """The steady-state 26Al/27Al in the gas of a star's galaxy at its [Fe/H] (`birth_radioisotopes.mean_ratio`), the default galaxy's when there is no shape."""
    from planetgen.galaxy import shape_draw
    shape = shape_draw.default_shape(None) if shape is None else shape
    aluminium = sa.PROTO_MASS_FRACTION["Al"] * 10.0 ** (0.0 if fe_h is None else float(fe_h))
    return birth_radioisotopes.mean_ratio("Al26", _depletion_time_yr(shape), aluminium)


def star_al26_mean(body):
    """`al26_mean` for the host of `body`: its galaxy (`system_config.GALAXY_SHAPE`, the default galaxy's when there is none) and its [Fe/H]."""
    star = getattr(body, "star", None)
    shape = getattr(getattr(body, "system_config", None), "GALAXY_SHAPE", None)
    return al26_mean(shape, getattr(star, "metallicity_fe_h", None))


def draw_content(offset_dex=None, rng=draw, refractory_enrichment=1.0, th_u_factor=None, al26_mean_ratio=None, al26_ratio=None):
    """
    A body's `Content`. Every draw is made whatever the branch, in this order: the host-less offset, the r-process production ratio, the
    potassium depletion and the seed of the aluminium-26 additions, so the stream after it is the same for every body.

    Args:
        offset_dex (float or None): The host star's Th/U offset (`star_offset_dex`); `None` takes the host-less draw.
        rng: The draw source (`planetgen.util.draw` by default).
        th_u_factor (float or None): The enrichment `F` of the host's birth gas (`rprocess_enrichment`, at least 1); `None` takes that of
            the typical star of the default galaxy (a stand-alone call with no galaxy).
        al26_mean_ratio (float or None): The steady-state 26Al/27Al of the host's gas (`al26_mean`); `None` takes the default galaxy's at solar [Fe/H].
        al26_ratio (float or None): The host star's birth ratio (`star_al26_ratio`), shared by its bodies; `None` draws the body's own about the mean (the seed draw is made either way).
        refractory_enrichment (float): How many times the Earth's share the body's refractory elements (Al, Th, U) are in its
            silicate; the body's Al2O3 over the bulk silicate Earth's (part E 2.8b: Th and U scale with it; `rock_of`), 1 for a body with none.

    Returns:
        Content: At formation, in the bulk silicate.
    """
    host_less = truncated_offset_dex(rng.random())
    low, high = rprocess_enrichment.PRODUCTION_RATIO_RANGE
    production_ratio = rng.uniform(low, high)
    k_depletion = math.exp(rng.uniform(math.log(tuning.RADIOGENIC_K_DEPLETION_RANGE[0]), math.log(tuning.RADIOGENIC_K_DEPLETION_RANGE[1])))
    al26_seed = rng.random()                  # one draw whatever the star: the additions are made on a stream of their own
    offset = host_less if offset_dex is None else offset_dex
    _bulk_k, bulk_th, _bulk_u = BULK_SILICATE_EARTH
    th_now = bulk_th * 10.0 ** offset * refractory_enrichment
    if th_u_factor is None:
        th_u_factor = default_th_u_factor()
    u_now = th_now / th_u_ratio_now(production_ratio, th_u_factor)
    if al26_ratio is None:
        al26_ratio = birth_radioisotopes.sample_ratio("Al26", al26_mean() if al26_mean_ratio is None else al26_mean_ratio, draw.Stream(f"al26-birth:{al26_seed!r}"))
    k_now = k_depletion * CI_CHONDRITE_K_OVER_U * u_now
    return Content(
        k0_ppm=k_now * formation_mix("K")[1], th0_ppm=th_now * formation_mix("Th")[1], u0_ppm=u_now * formation_mix("U")[1],
        al26_ratio=al26_ratio)


def generate_content(body, refractory_enrichment=1.0):
    """Sets `BODY_FIELDS` on a rocky planet or moon from its host star's composition (`star_offset_dex`) and its own draws; a gas
    giant has none (every field `None`)."""
    if body.body_type == "g":
        for name in BODY_FIELDS:
            setattr(body, name, None)
        return
    content = draw_content(star_offset_dex(body.star), refractory_enrichment=refractory_enrichment, th_u_factor=star_th_u_factor(body),
                           al26_mean_ratio=star_al26_mean(body), al26_ratio=star_al26_ratio(body.star, getattr(getattr(body, "system_config", None), "GALAXY_SHAPE", None)))
    body.c_k0_ppm, body.c_th0_ppm, body.c_u0_ppm = content.k0_ppm, content.th0_ppm, content.u0_ppm
    body.al26_ratio = content.al26_ratio


# --- The surface rock and its activity --------------------------------------------------------------------------------------

def surface_content(content, rock_ppm):
    """The bulk silicate's `content` with its K, Th and U raised to the surface rock `rock_ppm`, `(K, Th, U)` ppm today (the area-weighted
    reservoirs of the body's ground, `surface_composition.land_rock_ppm`): each element as many times the bulk silicate's as the rock is
    over the bulk silicate Earth. The formation concentrations scale together, so every isotope is as rich as its element."""
    k, th, u = (value / bulk for value, bulk in zip(rock_ppm, BULK_SILICATE_EARTH))
    return Content(content.k0_ppm * k, content.th0_ppm * th, content.u0_ppm * u, content.al26_ratio)


def activity_bq_per_kg(content, age_gy):
    """`{isotope: Bq/kg of the rock}` of the parents the dose reads (K-40, U-238, Th-232) at `age_gy`: the isotope's mass fraction
    times its specific activity (`elements.specific_activity`). The 226Ra that makes radon is in equilibrium with the 238U."""
    masses = isotope_mass_fractions(content, age_gy)
    return {name: masses[name] * elements.specific_activity(name) for name in GAMMA_NGY_H_PER_BQ_KG}


def gamma_ngy_per_hour(content, age_gy):
    """The dose rate in air at 1 m over `content` rock, nGy/h (infinite flat source; `GAMMA_NGY_H_PER_BQ_KG`)."""
    return math.fsum(GAMMA_NGY_H_PER_BQ_KG[name] * value for name, value in activity_bq_per_kg(content, age_gy).items())


