# planetgen/physics/stellar_composition.py

"""
What a star is made of, derived from the values already stored (GEN.440, docs/design/composition/stellar-composition.md sections 3, 4, 7 and 9)
======================================================================================================================================

`star_makeup` returns the composition of any star, brown dwarf or remnant from its kind, its stored [Fe/H], its age, its mass and a
key (and a white dwarf's temperature). It adds no column to the 75 million stars and changes no stored value: the branch, the four offsets and the rest of the
star's draws come from a stream of their own keyed on the star (`abundance:{key}`, separate from every other draw, so no earlier
draw moves), always consumed in the same fixed order, so the same star has the same composition at any worker count and an
edited [Fe/H], age or mass re-derives it with nothing to repair. The relations are those of `stellar_abundances` (GEN.439); this
module chooses the star's kind, draws, and applies the phase.

The key, in the order the first that applies wins (stellar-composition.md sections 3 and 9):

* `group` (a multi-star system's primary, or a cluster's seed): all members share the group's draws, so they share the branch
  and the offsets; each member still takes its own phase (a giant sibling differs from a dwarf one) and its own remnant.
* `bright_seed`: a bright star's own seed.
* `uid`: the star's 80-bit identifier.
* otherwise the quantities that already seed the [Fe/H] stream: mass, age and distance from the centre.

The kinds and what each returns:

=================== ===================================================================================================
protostar, PMS       the natal composition of the cloud, as a main-sequence star
main sequence        the birth composition (a fully convective star is homogeneous; a star's surface is its birth value);
                     the MIST luminosity, temperature, radius and lifetime at its mass, age and [Fe/H] (downloaded)
subgiant             the birth composition
red giant            the birth composition with the first dredge-up applied (`stellar_abundances.after_first_dredge_up`)
agb, supergiant      the birth composition; the shifts are open (the dredge-up grid stops at 2.5 solar masses)
wolf-rayet           the surface fractions of its subtype (hydrogen-free; carbon drawn within the subtype's range for WC)
brown dwarf          the birth composition; the deuterium, lithium and hydrogen limits at its mass and Z; the radius (stand-in)
white dwarf          the core type, an atmosphere draw (DA with probability 0.80), the pollution flag of a cool one, the radius and
                     the layers from the Bedard et al. sequences at its temperature (downloaded)
neutron star         the outer-crust nuclei of one of 20 nuclear mass models (downloaded); no composition of its own
black hole           none (no hair): mass, spin and the horizon quantities by formula
=================== ===================================================================================================

What the planet draws (GEN.428), the description layer and the radiation model read: `StarMakeup.birth` (X, Y, Z, the element
fractions and, through `stellar_abundances.mass_ratio`, Fe/Si, Mg/Si, C/O and O/Si). The birth composition is what is returned
and not the photosphere (stellar-composition.md D14); an age law for the photospheric correction (0.06 to 0.07 dex over the Sun's life)
has not been read.

The downloaded files (`planetgen.datasets`, GEN.548) are required (GEN.590, Boss's answer on the card): a kind that reads them
raises `fetch.DataMissing` when they were not downloaded, and bulk generation does not start without them
(`run_galaxy.require_datasets`), so a derived makeup never differs between a machine with the files and one without. The one
exception is a run told in so many words to disregard missing data (GEN.599, `fetch.set_disregard_missing_data`): each kind
then takes its dataset's labelled fallback, says so in its `open_items` and is listed in `fetch.stand_ins_used()`. A main-sequence
star has no fallback value (no luminosity, temperature, radius or lifetime), a white dwarf the radius of the zero-temperature
Chandrasekhar table and no layers, a neutron star the [R0722] outer crust.
"""

import math
import re
from dataclasses import dataclass, field

from planetgen import tuning
from planetgen.datasets import baraffe_isochrones, fetch, isochrones, neutron_star_matter, white_dwarf_models
from planetgen.physics import alpha_tracks, constants, elements, main_sequence, wolf_rayet
from planetgen.physics import stellar_abundances as sa
from planetgen.util import draw

PROTOSTAR = "protostar"
PRE_MAIN_SEQUENCE = "pre-main-sequence"
MAIN_SEQUENCE = "main-sequence"
SUBGIANT = "subgiant"
RED_GIANT = "red-giant"
ASYMPTOTIC_GIANT = "asymptotic-giant"
SUPERGIANT = "supergiant"
WOLF_RAYET = "wolf-rayet"
BROWN_DWARF = "brown-dwarf"
WHITE_DWARF = "white-dwarf"
NEUTRON_STAR = "neutron-star"
BLACK_HOLE = "black-hole"
KINDS = (PROTOSTAR, PRE_MAIN_SEQUENCE, MAIN_SEQUENCE, SUBGIANT, RED_GIANT, ASYMPTOTIC_GIANT, SUPERGIANT, WOLF_RAYET,
         BROWN_DWARF, WHITE_DWARF, NEUTRON_STAR, BLACK_HOLE)
"""tuple: Every kind of star, brown dwarf and remnant this module gives a composition."""

WITH_BIRTH_COMPOSITION = tuple(kind for kind in KINDS if kind not in (NEUTRON_STAR, BLACK_HOLE))
"""tuple: The kinds that carry a metallicity, so a birth composition. Neutron stars and black holes store none."""

SUPERGIANT_CLASSES = ("0", "IA+", "IA", "IAB", "IB")
"""tuple: The Yerkes classes of the supergiants and hypergiants (`generation.star`)."""

DA_SHARE = tuning.STELLAR_WD_DA_SHARE
"""float: The share of white dwarfs with a hydrogen atmosphere (DA), Xu, [R4798]."""

POLLUTED_SHARE = sum(tuning.STELLAR_WD_POLLUTED_SHARE) / 2.0
"""float: The share of cool white dwarfs polluted by accreted rocky debris (GEN.253): the middle of the 25 to 50 percent the
observations give (`tuning.STELLAR_WD_POLLUTED_SHARE`), which is also the registry's default for `white-dwarf-pollution`."""

CARBON_STAR_MINIMUM_MASS_SOL = ((0.010, 1.4), (0.020, 1.5))
"""tuple: `(Z, least initial mass in solar masses)` for third dredge-up to make C/O above 1 ([R2244], second-hand,
[R0925]); linear in Z between, 1.425 at Z 0.0125. Whether a given AGB star has run enough thermal pulses is the carbon-star
flag's (GEN.249), and the C/O it reaches is open."""


class StarMakeupError(ValueError):
    """A kind or a request this module does not cover."""


@dataclass(frozen=True)
class StarDraws:
    """
    A star's draws, in the fixed order they are consumed in; every one is always drawn.

    Attributes:
        normals (tuple): Five standard normals for the C, O, Mg, Si and Ca residuals (the calcium one is drawn last of all, so no earlier draw moved
            when calcium gained its scatter, GEN.988).
        u_enhanced (float): Decides a carbon-enhanced metal-poor star.
        u_carbon (float): The size of its carbon excess.
        eu_normal (float): The [Eu/Fe] residual in units of its scatter.
        th_eu_normal (float): The [Th/Eu] residual in units of its scatter.
        u_atmosphere (float): Decides a white dwarf's atmosphere.
        u_wolf_rayet (float): Places a Wolf-Rayet star's carbon inside its subtype's range.
        u_crust_model (float): Picks the nuclear mass model of a neutron star's outer crust.
        u_pollution (float): Decides whether a cool white dwarf is polluted (`white_dwarf_polluted`).
    """

    normals: tuple
    u_enhanced: float
    u_carbon: float
    eu_normal: float
    th_eu_normal: float
    u_atmosphere: float
    u_wolf_rayet: float
    u_crust_model: float
    u_pollution: float


def abundance_key(uid=None, mass=None, distance_ly=None, bright_seed=None, group=None):
    """
    The key of a star's draws (the stream is `abundance:` plus this): `group` (a system's primary, or `cluster:{seed}`) if
    given, else `bright:{seed}`, else the star's `uid`, else the star's birth mass and distance from the galactic centre (as their `repr`). The age is not in the key: a star's age rises
    once, after its planets are made (`Star.adjust_age_for_planets`), and the star's makeup and every planet's host ratios must read the
    same draws before and after (GEN.984).

    Raises:
        StarMakeupError: Nothing identifies the star.
    """
    if group is not None:
        return f"{group}"
    if bright_seed is not None:
        return f"bright:{bright_seed}"
    if uid is not None:
        return f"{uid}"
    if mass is None or distance_ly is None:
        raise StarMakeupError("a star's key needs a group, a bright-star seed, a uid, or its mass and distance")
    return f"{mass!r}:{distance_ly!r}"


def key_draws(key):
    """The draws of `key`: a stream of its own, `abundance:{key}`, drawn in the order of `StarDraws`."""
    stream = draw.Stream(f"abundance:{key}")
    stream.random()    # unread: the branch left the key (`stellar_abundances.branch_of`) and every later draw keeps its place
    normals = tuple(stream.gauss(0.0, 1.0) for _ in range(4))
    rest = (stream.random(), stream.random(), stream.gauss(0.0, 1.0), stream.gauss(0.0, 1.0),
            stream.random(), stream.random(), stream.random(), stream.random())
    return StarDraws(normals + (stream.gauss(0.0, 1.0),), *rest)


def white_dwarf_polluted(teff_k, u_pollution, share=POLLUTED_SHARE):
    """
    Whether a white dwarf is polluted (GEN.253): `u_pollution` under `share` for a cool one (cooler than
    `constants.COOL_WHITE_DWARF_BELOW_K`, the floor of the hot white dwarfs), never for a hot one, whose pollution the
    observations do not give. The share is of the cool white dwarfs, so the same star gets the same answer whatever it is
    drawn for.
    """
    return teff_k < constants.COOL_WHITE_DWARF_BELOW_K and u_pollution < share


@dataclass(frozen=True)
class StarMakeup:
    """
    The composition of one star or remnant.

    Attributes:
        kind (str): One of `KINDS`.
        branch (str or None): The abundance sequence, `None` for a black hole or neutron star.
        birth (StellarComposition or None): The composition the star formed with.
        surface (StellarComposition or None): The composition at the surface, after the changes of its phase; the birth
            composition for kinds that change none.
        draws (StarDraws): The star's draws (for the radiogenic birth abundances and the planets).
        detail (dict): What is particular to the kind (see the module docstring); every value is plain data.
        open_items (tuple): What the design leaves open for this kind, so a description does not claim it.
    """

    kind: str
    branch: str | None
    birth: object
    surface: object
    draws: StarDraws
    detail: dict = field(default_factory=dict)
    open_items: tuple = ()


def kind_of(yerkes_class, mass_sol=None, phase=None, z=None):
    """
    The kind of a star from its stored Yerkes class (as `generation.star` and `stellar_evolution.evolve_star` write it), its
    `phase` ("horizontal-branch" or "asymptotic-giant" for a class III star) and its mass. A star below the hydrogen-burning
    limit at its Z is a brown dwarf whatever its class says.

    Raises:
        StarMakeupError: The class is none this module knows.
    """
    if mass_sol is not None and mass_sol < sa.hydrogen_burning_limit_sol(sa.PROTO_Z if z is None else z) \
            and yerkes_class not in ("VII", "D"):
        return BROWN_DWARF
    if yerkes_class == "protostar":
        return PROTOSTAR
    if yerkes_class == "PMS":
        return PRE_MAIN_SEQUENCE
    if yerkes_class in ("V", "VI"):
        return MAIN_SEQUENCE
    if yerkes_class == "IV":
        return SUBGIANT
    if yerkes_class in ("II", "III"):
        return ASYMPTOTIC_GIANT if phase == "asymptotic-giant" else RED_GIANT
    if yerkes_class in SUPERGIANT_CLASSES:
        return SUPERGIANT
    if yerkes_class == wolf_rayet.CLASS:
        return WOLF_RAYET
    if yerkes_class in ("VII", "D"):
        return WHITE_DWARF
    raise StarMakeupError(f"no kind for the Yerkes class {yerkes_class!r}")


def carbon_star_minimum_mass_sol(z):
    """The least initial mass at which an AGB star of metal fraction `z` can become a carbon star, linear in Z between the two
    tabulated metallicities and held outside them."""
    (z_low, m_low), (z_high, m_high) = CARBON_STAR_MINIMUM_MASS_SOL
    share = min(max((z - z_low) / (z_high - z_low), 0.0), 1.0)
    return m_low + share * (m_high - m_low)


def wolf_rayet_row(subtype, teff_k=None):
    """
    The row of `tuning.STELLAR_WOLF_RAYET` a Wolf-Rayet star takes: a row's own name ("WN-early", "WN-late", "WN/WC", "WC", "WO"),
    or the generator's subtype with its subclass number ("WN6", "WC4", "WO2", as the spectral text writes it), or the subtype
    alone ("WN", "WC", "WO") with the star's effective temperature, from which the subclass follows (`wolf_rayet.subclass`). A WN
    star is early up to `tuning.WR_WN_EARLY_MAX_SUBCLASS` and late above it.

    Raises:
        StarMakeupError: The name is none of these, or a WN subtype has neither a subclass number nor a temperature.
    """
    if subtype in tuning.STELLAR_WOLF_RAYET:
        return subtype
    match = re.fullmatch(r"(WN|WC|WO)(\d+)?", str(subtype))
    if match is None:
        raise StarMakeupError(f"no Wolf-Rayet subtype {subtype!r}; expected one of {sorted(tuning.STELLAR_WOLF_RAYET)}, "
                              "or WN, WC or WO with its subclass number")
    name, number = match.group(1), match.group(2)
    if name != "WN":
        return name
    if number is None:
        if teff_k is None:
            raise StarMakeupError("a WN star needs its subclass number or its effective temperature to be early or late")
        number = wolf_rayet.subclass(name, teff_k)
    return "WN-early" if int(number) <= tuning.WR_WN_EARLY_MAX_SUBCLASS else "WN-late"


def wolf_rayet_surface(subtype, u, teff_k=None):
    """The surface mass fractions of a Wolf-Rayet star of `subtype` (`wolf_rayet_row`): every range of the subtype's row taken at the
    position `u` (0 to 1) inside it. The remainder is helium (and some oxygen in a WN or WC star, the split open); a WO star's oxygen is its own range."""
    row = tuning.STELLAR_WOLF_RAYET[wolf_rayet_row(subtype, teff_k)]
    fractions = {name: (value[0] + u * (value[1] - value[0]) if isinstance(value, tuple) else value) for name, value in row.items()}
    fractions["He" if "O" in fractions else "He and O"] = 1.0 - sum(fractions.values())
    return fractions


def _neutron_star_crust(draws):
    """The outer crust of a neutron star: the nuclei, surface inward, of the nuclear mass model `draws.u_crust_model` picks among
    the downloaded ones (the 20 models are the spread of the theory, so one is drawn per star), and the density at which
    neutron drip ends it. A run that disregards missing data takes the one sequence of [R0722] instead."""
    if fetch.locate(neutron_star_matter.DATASET).use_stand_in("the [R0722] outer crust"):
        return {"crust_model": None, "outer_crust": tuple(nucleus for nucleus, _density in sa.NEUTRON_STAR_OUTER_CRUST),
                "neutron_drip_g_cm3": sa.NEUTRON_DRIP_DENSITY_G_CM3, "radius_band_km": tuning.STELLAR_NS_RADIUS_BAND_KM}
    models = neutron_star_matter.crust_models()
    model = models[min(int(draws.u_crust_model * len(models)), len(models) - 1)]
    layers = neutron_star_matter.outer_crust(model)
    return {"crust_model": model, "outer_crust": tuple(f"{elements.element(layer.z).symbol}-{layer.a}" for layer in layers),
            "neutron_drip_g_cm3": neutron_star_matter.mass_density_g_cm3(layers[-1].upper_density_fm3),
            "radius_band_km": tuning.STELLAR_NS_RADIUS_BAND_KM}


def _grid_radius_solar(mass_sol, teff_k, envelope):
    """
    A white dwarf's radius from the cooling sequences, solar radii. The initial-final mass relation reaches 1.33 solar masses and the
    sequences stop at `white_dwarf_models.CRYSTAL_MASS_SOL` (1.30): a heavier star reads the 1.30 sequence at its temperature, scaled by the
    zero-temperature mass-radius relation's ratio of the two masses (the shrinking toward the Chandrasekhar mass).
    """
    top = white_dwarf_models.CRYSTAL_MASS_SOL
    if mass_sol <= top:
        return white_dwarf_models.radius_solar(mass_sol, teff_k, envelope)
    try:
        return (white_dwarf_models.radius_solar(top, teff_k, envelope)
                * sa.white_dwarf_radius_solar(mass_sol) / sa.white_dwarf_radius_solar(top))
    except sa.StellarAbundanceError:
        return None  # heavier than any white dwarf the relation makes


def star_makeup(kind, fe_h, age_gy, mass_sol, key, wolf_rayet_subtype=None, spin=0.0, teff_k=None):
    """
    The composition of a star or remnant.

    Args:
        kind (str): One of `KINDS` (`kind_of` finds it from a stored Yerkes class).
        fe_h (float or None): The stored [Fe/H]; `None` for a neutron star or black hole.
        age_gy (float): The star's age, Gy.
        mass_sol (float): Its present mass, solar masses (a white dwarf's final mass, a giant's near-initial one).
        key (str): `abundance_key(...)`: the draws are `key_draws(key)`.
        wolf_rayet_subtype (str or None): The subtype, for a Wolf-Rayet star (`wolf_rayet_row`: the spectral text's "WN6", or "WN" with `teff_k`).
        spin (float): A black hole's dimensionless spin.
        teff_k (float or None): A white dwarf's effective temperature, K (its stored one); a Wolf-Rayet star's, when its subtype has no subclass number.

    Raises:
        StarMakeupError: An unknown kind, or a kind that needs a [Fe/H], a subtype or a temperature and has none.
        planetgen.datasets.fetch.DataMissing: The kind reads a downloaded file (main sequence, white dwarf, neutron star), the
            file was not downloaded and the run was not told to disregard missing data.
    """
    if kind not in KINDS:
        raise StarMakeupError(f"unknown kind {kind!r}; expected one of {', '.join(KINDS)}")
    draws = key_draws(key)
    if kind == BLACK_HOLE:
        return StarMakeup(kind, None, None, None, draws, {"hole": sa.black_hole(mass_sol, spin)})
    if kind == NEUTRON_STAR:
        crust = _neutron_star_crust(draws)
        open_items = ("inner crust and core need the equation of state integrated (nothing does yet)",)
        if crust["crust_model"] is None:
            open_items += (f"the outer crust is the stand-in sequence: {neutron_star_matter.DATASET.fallback}",)
        return StarMakeup(kind, None, None, None, draws, crust, open_items=open_items)
    if fe_h is None:
        raise StarMakeupError(f"a {kind} needs a [Fe/H]")
    branch = sa.branch_of(fe_h, age_gy)
    offsets = sa.abundance_offsets(fe_h, branch, draws.normals, draws.u_enhanced, draws.u_carbon)
    birth = sa.composition(fe_h, offsets)
    surface, detail, open_items = birth, {}, ()
    if kind == RED_GIANT:
        surface = sa.after_first_dredge_up(birth, mass_sol)
        detail["dredge_up"] = sa.first_dredge_up(mass_sol, birth.z)
    elif kind in (ASYMPTOTIC_GIANT, SUPERGIANT):
        detail["can_be_carbon_star"] = kind == ASYMPTOTIC_GIANT and mass_sol >= carbon_star_minimum_mass_sol(birth.z)
        open_items = ("surface shifts of the later phases (the dredge-up grid stops at 2.5 solar masses)",
                      "the carbon-to-oxygen ratio of a carbon star")
    elif kind == MAIN_SEQUENCE:
        if fetch.locate(isochrones.DATASET).use_stand_in("no main-sequence values"):
            point, detail["main_sequence_lifetime_gy"] = None, None
            open_items = (f"the isochrones are missing: {isochrones.DATASET.fallback}",)
        else:
            track = alpha_tracks.track_fe_h(fe_h, age_gy)    # the track the star's own record reads (GEN.1024)
            point = isochrones.main_sequence_point(mass_sol, age_gy, track)
            detail["main_sequence_lifetime_gy"] = isochrones.main_sequence_lifetime_gy(mass_sol, track)
            if point is None and tuning.IMF_BREAKS_SOL[0] <= mass_sol < tuning.MAIN_SEQUENCE_GRID_MASS_RANGE_SOL[0]:
                # GEN.652: a mass under the grid reads the extension the star's own record reads (`main_sequence.state`).
                row = main_sequence.state(mass_sol, age_gy, track)
                point = (math.log10(row.luminosity_sol), math.log10(row.temperature_k), math.log10(row.radius_sol))
            elif point is None:
                open_items = ("the track does not cover this mass or age (0.1 solar masses and up, 0.1 Myr to about 20 Gy)",)
        detail["main_sequence"] = None if point is None else dict(zip(("log_luminosity_sol", "log_teff_k", "log_radius_sol"), point))
    elif kind == WOLF_RAYET:
        if wolf_rayet_subtype is None:
            raise StarMakeupError("a Wolf-Rayet star needs its subtype")
        detail["surface_fractions"] = wolf_rayet_surface(wolf_rayet_subtype, draws.u_wolf_rayet, teff_k)
        open_items = ("the helium and oxygen split",)
    elif kind == BROWN_DWARF:
        mass_jupiter = mass_sol / sa.JUPITER_MASS_SOL
        deuterium = sa.deuterium_burning_limit_mj(birth.z, birth.y)
        detail.update(
            mass_jupiter=mass_jupiter, deuterium_limit_mj=deuterium, burned_deuterium=mass_jupiter > deuterium,
            burned_lithium=mass_sol > sa.lithium_burning_limit_sol(),
            burns_hydrogen=mass_sol > sa.hydrogen_burning_limit_sol(birth.z))
        if fetch.locate(baraffe_isochrones.DATASET).use_stand_in("no brown dwarf radius"):
            detail["radius_solar"] = None
            open_items = (f"the Baraffe table is missing: {baraffe_isochrones.DATASET.fallback}", "radius at other metallicities")
        else:
            try:
                detail["radius_solar"] = sa.brown_dwarf_radius_solar(mass_sol, age_gy)
            except sa.StellarAbundanceError:
                detail["radius_solar"] = None
            open_items = ("radius at other metallicities",)
    elif kind == WHITE_DWARF:
        if teff_k is None:
            raise StarMakeupError("a white dwarf needs its effective temperature")
        detail["core"] = sa.white_dwarf_core(mass_sol)
        detail["atmosphere"] = "DA" if draws.u_atmosphere < DA_SHARE else "non-DA"
        envelope = "thick" if detail["atmosphere"] == "DA" else "thin"
        detail["envelope"] = envelope
        detail["polluted"] = white_dwarf_polluted(teff_k, draws.u_pollution)
        if fetch.locate(white_dwarf_models.DATASET).use_stand_in("the Chandrasekhar mass-radius table"):
            try:
                detail["radius_solar"] = sa.white_dwarf_radius_solar(mass_sol)
            except sa.StellarAbundanceError:
                detail["radius_solar"] = None
            detail["layers"] = None
            open_items = (f"the white-dwarf sequences are missing: {white_dwarf_models.DATASET.fallback}",)
        else:
            detail["radius_solar"] = _grid_radius_solar(mass_sol, teff_k, envelope)
            detail["layers"] = (white_dwarf_models.layer_mass_fractions(mass_sol, envelope)
                                if detail["core"] == "carbon-oxygen" else None)
            open_items = ("the carbon-to-oxygen split of the core (the sequences assume an even one)",
                          "the layers and the radius of a helium-core or oxygen-neon white dwarf (the sequences are carbon-oxygen)")
    return StarMakeup(kind, branch, birth, surface, draws, detail, open_items)
