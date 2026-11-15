# planetgen/physics/ore_deposits.py

"""
The magmatic ore deposits of a rocky body's crust (GEN.1380, the first slice of GEN.1379's design)
=================================================================================================

`docs/design/composition/ore-processes.md` lists 19 deposit types and gives each a gate on the state the generator already stores. This module builds the first slice: the three **magmatic**
types (section 5, rows 6, 7 and 8), the **three derived fields** every later slice reads (section 4), and the table of what a body stores for them (section 7). Later slices (GEN.1381 to
GEN.1383, GEN.1378) add rows to the same table.

**The derived fields** (all DESIGN stand-ins, ORE1 to ORE3):

* `crust_thickness_km`: the mean thickness of the crust, `D_lid / g * (T_ref / T_m)`, with the surface gravity `g` in Earth gravities, the mantle temperature `T_m` the interior model
  gives at the star's age now (`interior_heat.Heat.mantle_temperature_k`) and `D_lid` 17 km for a mobile lid (thin oceanic crust, thicker continents) and 22 km for a stagnant one. A small body
  has a thick crust, a hot mantle a thin one; Mars (g 0.38) comes out at about 58 km against the 50 km recalled and Venus (g 0.90) at 24 km against 30.
* `craton_fraction`: 0 unless the lid is mobile, the body is at least 2 Gyr old and the ground holds a granitoid share; otherwise 0.3 of the crust, scaled by how far past 2 Gyr the age is
  (full at 4.5 Gyr) and by one draw of 0.5 to 1.5 on the body's own ore stream. Nothing reads it before GEN.1383 (kimberlite and unconformity uranium).
* `crust_sulfur_index`: how much sulfur the crust offers a rising mafic magma, `0.15 r + 0.65 w a`: `r` the star's S/Fe over the Sun's, `w` 1 where liquid water stands (sulfate
  reduction builds sedimentary sulfide; evaporites carry sulfate) and `a` the share of 2 Gyr the age has covered (sedimentary sulfide needs time). Dry bodies keep the primordial 0.15 `r`.

**The opportunity** of a magmatic type is the mantle the body has processed through melting, `degassed_share(released heat, mantle mass) * mantle mass`
(`surface_pressure.degassed_share`, the same measure the degassing of the air uses), over `MAGMA_REFERENCE_KG`, the processed mass of an Earth-like body at 4.5 Gyr. Below `MAGMA_MIN_RATIO`
of it a body has no large magma volume and no magmatic ore; an icy body (a mantle under an ice shell) has none.

**The gates** (rule 1 of the design: magma volume, not the lid): Ni-Cu-PGE sulfide needs the magma volume, a sulfur index of at least `SULFUR_GATE` and a mantle that is not oxidized
(sulfate, not sulfide, is the sulfur of an oxidized magma); chromitite needs the magma volume and a reduced or intermediate mantle; Fe-Ti-V oxide needs only the magma volume.

**Endowment and grade** (section 6). The number of deposits is drawn from a Poisson law with mean `deposits * opportunity * window` (`window` is the share of the body's life in which the
setting existed); each deposit's ore mass is a lognormal of the type's median and spread (the design file's per-deposit tonnage range); the endowment is their sum. A gate that holds and a draw
of no deposit is `absent_by_chance`. The grade is drawn uniformly inside the type's range and capped by what the host supplies: `concentration * crust_mean`, with the crustal mean of the key
element the body's stored bulk-silicate oxide times `crust_over_bulk` (the crust holds a fraction of the mantle's chromium, nickel and titanium). A cap under the bottom of the range is a body with
the dispersed element and no ore (`absent_by_state`), which is how a metal-poor star gives poorer ore. One stream per body and type (`stream_of`) keeps a type's draws from moving any other.

What is left out, and where it goes: the platinum-group elements of the sulfide type (the star set has none; GEN.1371 gives them, and the grade here is the nickel's); vanadium (the key
element of the oxide type is TiO2); the impact-melt route to sulfide that the design names (the collision record is not stored); a body's age windows are shares of its age, not the absolute
ages of the Earth's deposits.
"""

import functools
import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import bulk_draw, interior_heat, radiogenic, stellar_abundances, stellar_composition, surface_pressure, tidal_heating
from planetgen.util import draw

PRESENT, ABSENT_BY_STATE, ABSENT_BY_CHANCE = "present", "absent_by_state", "absent_by_chance"
STATUSES = (PRESENT, ABSENT_BY_STATE, ABSENT_BY_CHANCE)
"""tuple: The `status` flag of a deposit row (design section 6), so a page and the census can tell a dry planet from an unlucky one."""

NICKEL_COPPER_PGE, CHROMITITE, IRON_TITANIUM_VANADIUM = "ni_cu_pge_sulfide", "chromitite", "fe_ti_v_oxide"

COLUMN_FIELDS = ("crust_thickness_km", "craton_fraction", "crust_sulfur_index")
"""tuple: The three derived fields a rocky planet or moon stores (columns of `planets`, `moons` and `rogue_planets`)."""

ATTRIBUTE = "deposit_types"
"""str: The body attribute that holds the assessed types, a list of `ROW_FIELDS` dicts (the child tables `planet_deposit_types`, `moon_deposit_types` and `rogue_planet_deposit_types`)."""

ROW_FIELDS = ("deposit_type", "status", "endowment_kg", "grade", "key_element", "age_window_gyr_start", "age_window_gyr_end", "capped_by_host")
"""tuple: The columns of a deposit row after the owner id."""

# --- the derived fields (ORE1 to ORE3) ------------------------------------------------------------------------------
CRUST_KM_MOBILE, CRUST_KM_STAGNANT = 17.0, 22.0
"""float: `D_lid`, km: the mean crust of an Earth-gravity body at the reference mantle temperature (a mobile lid's thin oceanic crust and thick continents, 17 km; a stagnant lid's one thick plate, 22 km). Stand-in ORE1, tier T6, L3: two values fitted by eye to the Earth (about 17 km area-weighted), Mars (50 km) and Venus (30 km) recalled."""
REFERENCE_MANTLE_K = 2200.0
"""float: The mantle temperature (K) of the interior model for an Earth-mass body at 4.5 Gyr (the middle of the draws, 1960 to 2760 K). Stand-in ORE1."""
MIN_GRAVITY_G = 0.05
CRUST_KM_RANGE = (3.0, 150.0)
"""tuple: The thinnest and thickest crust (km) the law returns (the law stays monotone inside it)."""

CRATON_MIN_AGE_GY = 2.0
"""float: A craton needs a body at least this old (design section 4.2)."""
CRATON_AGE_FULL_GY = 4.5
CRATON_SHARE = 0.3
CRATON_SCATTER = (0.5, 1.5)
"""Stand-in ORE2, tier T7, L3: the design's 'a DESIGN fraction around 0.3 for an Earth-like body', scaled by age, with the draw's scatter."""

SULFUR_PRIMORDIAL, SULFUR_SEDIMENTARY = 0.15, 0.65
SULFUR_AGE_GY = 2.0
"""Stand-in ORE3, tier T7, L3: the primordial and sedimentary parts of the sulfur index (they sum to 0.8 for an Earth-like body) and the time sedimentary sulfide takes to build."""
LIQUID_WATER_STATES = ("surface ocean", "ice-covered ocean", "hycean", "brine seep")
"""tuple: The stored hydrosphere states (`hydrosphere.HYDROSPHERE_STATES`) in which liquid water stands, for the index's wetness."""

# --- the opportunity and the gates (ORE4) ---------------------------------------------------------------------------
MAGMA_REFERENCE_KG = 1.5e24
"""float: The mantle mass an Earth-like body has processed through melting by 4.5 Gyr (about 0.28 of its 5e24 kg mantle). Stand-in ORE4, tier T7, L2."""
MAGMA_MIN_RATIO = 0.05
"""float: The least processed mantle (over `MAGMA_REFERENCE_KG`) that gives a large magma volume; a body under about 0.07 Earth masses is below it, a Mars-mass one above. Stand-in ORE4, tier T7, L3."""
SULFUR_GATE = 0.3
"""float: The least sulfur index of the Ni-Cu-PGE type (design: 'sulfur index high'). Stand-in ORE4, L3."""
CONCENTRATING_REDOX = {NICKEL_COPPER_PGE: ("reduced", "intermediate"), CHROMITITE: ("reduced", "intermediate"), IRON_TITANIUM_VANADIUM: ("reduced", "intermediate", "oxidized")}
"""dict: The mantle redox classes (`mantle_redox`) that admit each type: sulfide and chromite in a reduced or moderate mantle, the Fe-Ti oxide in any."""


@dataclass(frozen=True)
class DepositType:
    """
    One magmatic type of the design's section 5.

    Attributes:
        id (str): The type id stored.
        key_element (str): The element or oxide the grade is of.
        grade_range (tuple): `(low, high)` mass fraction of the key element in the ore (the data file's range, RECALLED).
        window (tuple): `(start, end)` as shares of the body's age in which the setting existed.
        deposits (float): Expected count of deposits on an Earth-like body (stand-in).
        median_ore_kg, sigma (float): Median ore mass of one deposit and the lognormal spread (natural log).
        oxide (str): The stored bulk-silicate oxide column the key element comes from.
        element_per_oxide (float): Mass of the key element in a unit mass of that oxide.
        crust_over_bulk (float): The crust's mean of the element over the bulk silicate's.
        concentration (float): The factor by which the process concentrates the crustal mean into ore.
    """

    id: str
    key_element: str
    grade_range: tuple
    window: tuple
    deposits: float
    median_ore_kg: float
    sigma: float
    oxide: str
    element_per_oxide: float
    crust_over_bulk: float
    concentration: float


TYPES = (
    DepositType(NICKEL_COPPER_PGE, "Ni", (0.005, 0.03), (0.30, 0.90), 60.0, 3.0e10, 1.5, "rock_nio_fraction", 0.7858, 0.024, 320.0),
    DepositType(CHROMITITE, "Cr2O3", (0.30, 0.60), (0.15, 0.85), 40.0, 1.0e11, 2.0, "rock_cr2o3_fraction", 1.0, 0.044, 3300.0),
    DepositType(IRON_TITANIUM_VANADIUM, "TiO2", (0.10, 0.35), (0.35, 0.75), 50.0, 3.0e11, 1.0, "rock_tio2_fraction", 1.0, 3.5, 45.0),
)
"""tuple: The three magmatic types. Grade ranges, concentration factors (Ni 320, Cr 3300) and per-deposit tonnage (1e6 to 1e9 t, 1e7 to 1e10 t, 1e8 to 1e9 t) are the design file's, RECALLED; the
counts, the windows and the crust-over-bulk shares are stand-ins ORE5 to ORE7 (T6, L2 to L3): the crust-over-bulk of Ni (0.024, the crust's 47 ppm over the mantle's 1960 ppm) and Cr (0.044,
set so the Earth's bulk-silicate Cr2O3 of 0.38 percent caps the grade at the middle of the range) and TiO2 (3.5, the crust's 0.7 percent over the mantle's 0.2) are RECALLED shares, the factor of Ti (45,
ilmenite's 32 percent over the crust's 0.7) is derived from them."""
TYPE_BY_ID = {kind.id: kind for kind in TYPES}


# --- the state a body offers ----------------------------------------------------------------------------------------
@dataclass(frozen=True)
class OreState:
    """
    What the ore types read of a rocky body.

    Attributes:
        gravity_g (float): Surface gravity in Earth gravities.
        lid (str): `interior_heat.MOBILE` or `interior_heat.STAGNANT`.
        mantle_k (float): The mantle temperature now.
        age_gy (float): The body's age.
        magma_ratio (float): The processed mantle over `MAGMA_REFERENCE_KG`.
        redox (str or None): The mantle redox class.
        wet (bool): Liquid water stands on the surface or in the crust.
        sulfur_ratio (float): The star's S/Fe over the Sun's.
        felsic_share (float): The areal share of granitoid ground.
        oxides (dict): The stored bulk-silicate oxide fractions by column name.
        icy (bool): The body is icy.
    """

    gravity_g: float
    lid: str
    mantle_k: float
    age_gy: float
    magma_ratio: float
    redox: str
    wet: bool
    sulfur_ratio: float
    felsic_share: float
    oxides: dict
    icy: bool


@functools.lru_cache(maxsize=4096)
def _over_iron_of(symbol, fe_h, age_gy, key):
    draws = stellar_composition.key_draws(key)
    branch = stellar_abundances.branch_of(fe_h, age_gy)
    offsets = stellar_abundances.abundance_offsets(fe_h, branch, draws.normals, draws.u_enhanced, draws.u_carbon)
    ratio = stellar_abundances.number_ratio
    return ratio(offsets, symbol, "Fe") / ratio({}, symbol, "Fe")


def over_iron(star, symbol):
    """
    The star's `symbol`/Fe by number over the reference Sun-like mix's (1 for a star with no abundances: a compact remnant), from the same draws as its other elements. Only the elements of
    `stellar_abundances` have a ratio; a metal the set does not hold (gold, lead) follows the iron: 1.
    """
    fe_h = getattr(star, "metallicity_fe_h", None)
    if fe_h is None or symbol not in tuning.STELLAR_SOLAR_LOG_EPS:
        return 1.0
    return _over_iron_of(symbol, fe_h, star.age, bulk_draw.star_key(star))


def _sulfur_over_iron(star):
    """The star's S/Fe by number over the reference Sun-like mix's."""
    return over_iron(star, "S")


def state_of(body):
    """The `OreState` of a rocky planet or moon, or `None` for a body with no solid interior stored (a giant, a belt, a row from before GEN.699). No draws."""
    if getattr(body, "body_type", None) in ("g", "a") or radiogenic.content_of(body) is None:
        return None
    found = interior_heat.heat_of(body)
    mantle = interior_heat.mantle_of(body.mass, body.radius, radiogenic.rock_of(body)) if found is not None else None
    if found is None or mantle is None:
        return None
    processed = mantle.dry_kg * surface_pressure.degassed_share(found.released_j, mantle.dry_kg)
    wet = getattr(body, "hydrosphere", None) in LIQUID_WATER_STATES
    oxides = {name: (getattr(body, name, None) or 0.0) for name in {kind.oxide for kind in TYPES}}
    return OreState(gravity_g=body.gravity, lid=interior_heat.lid_of(body), mantle_k=found.mantle_temperature_k, age_gy=getattr(body.star, "age", None) or tuning.ROGUE_EARTH_REFERENCE_AGE_GY,
                    magma_ratio=processed / MAGMA_REFERENCE_KG, redox=getattr(body, "mantle_redox", None), wet=wet, sulfur_ratio=_sulfur_over_iron(body.star),
                    felsic_share=getattr(body, "ground_felsic_fraction", None) or 0.0, oxides=oxides, icy=tidal_heating.is_icy(body))


def stream_of(body, kind):
    """The body's own stream for one type (or `"craton"`): seeded from what is final about the body and the type, so no draw of the run's stream moves and adding a type moves no other."""
    return body_stream("ore", body, kind)


def body_stream(prefix, body, kind):
    """
    The stream of one deposit pass (`prefix`: `ore`, `sediment`, `carbide`, `hydrothermal`) for one type on one body, keyed to the body's birth (`draw.body_stream_of`, GEN.1547:
    not its radius, orbit or the star's age) and the type, so no draw of the run's stream moves and adding a type moves no other. Shared by the four deposit modules (GEN.1667).
    """
    return draw.body_stream_of(f"{prefix}:{kind}", body, slot=str(bool(getattr(body, 'is_moon', False))))      # GEN.1547: keyed to the body's birth, not its radius, orbit or the star's age


# --- the derived fields ---------------------------------------------------------------------------------------------
def crust_thickness_km(state):
    """The mean crust thickness (km): `D_lid / g * (T_ref / T_m)`, clamped to `CRUST_KM_RANGE` (ORE1)."""
    base = CRUST_KM_MOBILE if state.lid == interior_heat.MOBILE else CRUST_KM_STAGNANT
    km = base / max(state.gravity_g, MIN_GRAVITY_G) * REFERENCE_MANTLE_K / max(state.mantle_k, 1.0)
    return min(max(km, CRUST_KM_RANGE[0]), CRUST_KM_RANGE[1])


def craton_fraction(state, scatter):
    """The craton share of the crust (ORE2): 0 unless the lid is mobile, the body at least `CRATON_MIN_AGE_GY` old and the ground holds a granitoid share; `scatter` is the body's draw in `CRATON_SCATTER`."""
    if state.lid != interior_heat.MOBILE or state.age_gy < CRATON_MIN_AGE_GY or state.felsic_share <= 0.0:
        return 0.0
    age_factor = min(1.0, (state.age_gy - CRATON_MIN_AGE_GY) / (CRATON_AGE_FULL_GY - CRATON_MIN_AGE_GY))
    return min(1.0, CRATON_SHARE * age_factor * scatter)


def crust_sulfur_index(state):
    """The sulfur the crust offers a mafic magma, 0 to 1 (ORE3): the primordial part scaled by the star's S/Fe plus the sedimentary part where liquid water has stood for time."""
    sedimentary = SULFUR_SEDIMENTARY * (1.0 if state.wet else 0.0) * min(1.0, state.age_gy / SULFUR_AGE_GY)
    return min(1.0, max(0.0, SULFUR_PRIMORDIAL * state.sulfur_ratio + sedimentary))


# --- the types -------------------------------------------------------------------------------------------------------
def admits(kind, state, sulfur_index):
    """Whether the gate of `kind` holds for `state` (design section 3, rule 1): the magma volume, the redox class and, for the sulfide, the sulfur index; never on an icy body."""
    if state.icy or state.magma_ratio < MAGMA_MIN_RATIO or state.redox not in CONCENTRATING_REDOX[kind.id]:
        return False
    return kind.id != NICKEL_COPPER_PGE or sulfur_index >= SULFUR_GATE


def grade_cap(kind, state):
    """The highest grade the host supplies: the concentration factor times the crustal mean of the key element (the stored oxide times its element share and the crust's share of the bulk)."""
    return min(1.0, kind.concentration * state.oxides[kind.oxide] * kind.element_per_oxide * kind.crust_over_bulk)


def poisson(stream, mean):
    """A Poisson draw: inversion for a small mean, the normal approximation above 30."""
    if mean <= 0.0:
        return 0
    if mean > 30.0:
        return max(0, round(stream.gauss(mean, math.sqrt(mean))))
    limit, count, product = math.exp(-mean), 0, stream.random()
    while product > limit:
        count += 1
        product *= stream.random()
    return count


def lognormal_sum(median_kg, sigma, count, stream):
    """The summed ore mass of `count` deposits of lognormal size, the median and the natural-log spread given (the mean of the law past 200 deposits, which the sum would only approach)."""
    if count > 200:
        return count * median_kg * math.exp(0.5 * sigma ** 2)
    return math.fsum(median_kg * math.exp(sigma * stream.gauss()) for _ in range(count))


def assess(kind, state, sulfur_index, stream):
    """
    The stored row of one type for a body: the gate, the cap on the grade, the count of deposits and their sum.

    Returns:
        dict: The `ROW_FIELDS`; `endowment_kg`, `grade` and the window are `None` where the type is absent.
    """
    row = {"deposit_type": kind.id, "status": ABSENT_BY_STATE, "endowment_kg": None, "grade": None, "key_element": kind.key_element,
           "age_window_gyr_start": None, "age_window_gyr_end": None, "capped_by_host": 0}
    if not admits(kind, state, sulfur_index):
        return row
    cap = grade_cap(kind, state)
    low, high = kind.grade_range
    if cap < low:
        return row
    window = kind.window[1] - kind.window[0]
    count = poisson(stream, kind.deposits * state.magma_ratio * (sulfur_index if kind.id == NICKEL_COPPER_PGE else 1.0) * window)
    drawn = stream.uniform(low, high)
    if count == 0:
        row["status"] = ABSENT_BY_CHANCE
        return row
    row.update(status=PRESENT, grade=min(drawn, cap), endowment_kg=lognormal_sum(kind.median_ore_kg, kind.sigma, count, stream), capped_by_host=int(drawn > cap),
               age_window_gyr_start=kind.window[0] * state.age_gy, age_window_gyr_end=kind.window[1] * state.age_gy)
    return row


def generate(body):
    """
    Sets `COLUMN_FIELDS` and `ATTRIBUTE` on a rocky planet or moon from its stored state (call after the ground; no draw of the run's stream). A giant, a belt or a row from before GEN.699
    gets `None` for each, so a stored row tells a body that was not assessed (no rows) from one with no ore (rows with `absent_by_state`).
    """
    for name in COLUMN_FIELDS:
        setattr(body, name, None)
    setattr(body, ATTRIBUTE, None)
    state = state_of(body)
    if state is None:
        return
    sulfur_index = crust_sulfur_index(state)
    body.crust_thickness_km = crust_thickness_km(state)
    body.craton_fraction = craton_fraction(state, stream_of(body, "craton").uniform(*CRATON_SCATTER))
    body.crust_sulfur_index = sulfur_index
    setattr(body, ATTRIBUTE, [assess(kind, state, sulfur_index, stream_of(body, kind.id)) for kind in TYPES])
