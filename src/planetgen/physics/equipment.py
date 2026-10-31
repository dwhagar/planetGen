# planetgen/physics/equipment.py

"""
Survival equipment (GEN.655): the gear a person needs on a world, as items that add up instead of one number
(docs/design/survival-equipment-fields.md, Boss's decisions of 2026-10-23).

Four slots, each a ladder of levels, and every item is **required** (past the physical limit) or **recommended**
(the fuzzy band before it, where a person survives but is impaired or at risk):

* **Body Protection**, the lowest suit that fits: 0 none, 1 clothing, 2 thermal suit, 3 sealed suit, 4 pressure suit,
  5 armoured suit, 6 vehicle only, with the stress that set it (`BODY_KINDS`) and a powered frame for high gravity.
* **Air Filter or Scrubber**, three parts that combine: dust, a gas filter of three grades and a CO2 scrubber.
* **Air Supply**: 0 ambient, 1 gas-mix mask, 2 own gas, 3 mixed gas, 4 isolated.
* **Radiation Protection**, a five-step ladder (`protection_limits.RADIATION_STEPS`): shielded garment, storm vest,
  hardened suit, then shielded shelter and hardened habitat. Each step is a layer of shielding laid on the world's own
  air column and run through the dose curves of `radiation`, against the stored parts of the dose.

The gear follows the physics wherever it is stricter than the score's bands (Boss, decision B): the PHI-6 scores are
unchanged. Only what is needed is ever listed (decision E); a body's text and API give a Required line and a Recommended line
(`lines`), and a list with nothing in it is left off (Boss 2026-10-24).

Everything here is a pure function of a `habitability.World` (`compute`, from which the stored levels come) or of
the stored levels alone (`items`, `required_names`, `recommended_names`, `code`, `summary`).
"""

from collections import namedtuple
from dataclasses import replace

from planetgen.physics import constants, habitability, radiation, radiation_belts, regional_temperature
from planetgen.physics.protection_limits import RADIATION_STEPS, TOXIC, limit

FIELDS = (
    "eq_body", "eq_body_kind", "eq_body_rec", "eq_body_rec_kind", "eq_rad_suit", "eq_rad_suit_rec", "eq_mobility",
    "eq_dust", "eq_gas_filter", "eq_co2_scrub", "eq_supply", "eq_supply_rec", "eq_shelter", "eq_shelter_rec",
    "shelter_column_g_cm2", "outside_hours",
)
"""tuple: The stored levels (a column of each on `planets` and `moons`). `eq_<x>_rec` is the recommended level, stored
only above the required one (else 0 or `None`). `eq_mobility`, `eq_dust` and `eq_co2_scrub` are 0 none, 1
recommended, 2 required. `eq_rad_suit` is the radiation suit step 0 to 3 (garment, vest, hardened suit) and
`eq_shelter` 0 none, 1 shielded shelter, 2 hardened habitat; `shelter_column_g_cm2` is the shielding the listed
shelter needs and `outside_hours` the hours outside a year the listed suit allows (`None` when it is not a limit)."""

BODY_KINDS = ("none", "cold", "heat", "pressure", "chemical", "diving", "uv", "gravity")
"""tuple: What set a Body Protection level, by `eq_body_kind` number."""

KIND_LETTER = {"none": "", "cold": "c", "heat": "h", "pressure": "p", "chemical": "x", "diving": "d", "uv": "u",
               "gravity": "g"}

KIND_DOMAIN = {"cold": "Temperature", "heat": "Temperature", "pressure": "Pressure", "diving": "Pressure",
               "chemical": "Air", "uv": "Radiation", "gravity": "Pressure"}
"""dict: The PHI-6 domain a body need comes from (gravity has none of its own: Pressure carries it)."""

BODY_NAMES = {
    (1, "cold"): "cold-weather clothing", (1, "heat"): "cooling clothing", (1, "uv"): "sun-cover clothing",
    (2, "cold"): "insulated suit", (2, "heat"): "liquid-cooled suit", (2, "diving"): "heated diving suit",
    (2, "uv"): "full-cover suit",
    (3, "cold"): "heated sealed suit", (3, "heat"): "heat-shielded suit", (3, "chemical"): "chemical-barrier suit",
    (4, "pressure"): "pressure suit", (5, "pressure"): "armoured suit", (5, "heat"): "armoured suit",
    (6, "heat"): "vehicle only", (6, "gravity"): "vehicle only", (6, "pressure"): "vehicle only",
}
"""dict: The body word of `(level, kind)`."""

BODY_DEFAULT = {1: "weather clothing", 2: "thermal suit", 3: "sealed suit", 4: "pressure suit", 5: "armoured suit",
                6: "vehicle only"}
SUPPLY_NAMES = {1: "gas-mix mask", 2: "self-contained air", 3: "mixed-gas rig", 4: "isolated air"}
GAS_FILTER_NAMES = {1: "gas filter", 2: "full-face gas filter", 3: "powered gas filter"}
SHELTER_NAMES = {1: "shielded shelter", 2: "hardened habitat"}
RAD_SUIT_NAMES = {1: "shielded garment", 2: "storm vest", 3: "hardened suit"}
RAD_STEP_KEYS = tuple(step.key for step in RADIATION_STEPS)

REQUIRED, RECOMMENDED = "required", "recommended"
SEVERITIES = (REQUIRED, RECOMMENDED)

_Item = namedtuple("_Item", "slot name severity domains")


class Item(_Item):
    """An item of the equipment line: its `slot`, `name`, `severity` and the PHI-6 `domains` it can come from, the likeliest
    first (GEN.702: a self-contained air supply may be for the air, the pressure or the cold; a hardened pressure suit is for
    the pressure and the radiation). `domain` is the first."""

    __slots__ = ()

    @property
    def domain(self):
        return self.domains[0]

    def rows(self, colours):
        """The PHI-6 rows the item is listed under: those of its domains whose colour (`colours`, domain to Blue, Green, Yellow
        or Red) is not Blue, since a Blue row asks for nothing; the first domain when all of them are Blue."""
        shown = [domain for domain in self.domains if colours.get(domain) != "Blue"]
        return shown or [self.domain]


def _item(slot, name, severity, *domains):
    return Item(slot, name, severity, tuple(domains))
"""An item of the equipment line: its `slot` (body, filter, supply, radiation), `name`, `severity` and the PHI-6
`domain` it comes from."""

R_J_MOL_K = 8.314462618
"""float: The gas constant, J/(mol K)."""

BEST_MIX_PO2_KPA = 40.0
"""float: The oxygen partial pressure of the best breathing mix (heliox), kPa."""

HAZARD_RATIOS = (10.0, 50.0, 1000.0)
"""tuple: The most a gas can exceed its long-term limit (the ratio to it) for a gas filter of grade 1 (half mask,
assigned protection factor 10), 2 (full face, 50) and 3 (powered, 1,000); past 1,000 the person breathes own gas."""

SHELTER_SEARCH_MAX_G_CM2 = 2000.0
"""float: How deep the shielding needed for a limit is searched, g/cm2."""

HOURS_PER_YEAR = 8766.0


def best_mix_density_g_l(pressure_kpa, temperature_c):
    """The density, g/L, of the best breathing mix (heliox with 40 kPa of oxygen, all oxygen when the pressure is
    lower) at the world's pressure and temperature: P M / (R T)."""
    if pressure_kpa <= 0:
        return 0.0
    oxygen = min(1.0, BEST_MIX_PO2_KPA / pressure_kpa)
    molar_g = oxygen * 31.998 + (1.0 - oxygen) * 4.0026
    return pressure_kpa * 1000.0 * molar_g * 1e-3 / (R_J_MOL_K * (temperature_c + constants.CELSIUS_ZERO_K))


# --- Body ---------------------------------------------------------------------------------

_KIND_RANK = ("uv", "heat", "cold", "diving", "chemical", "pressure", "gravity")


def _top(needs):
    """`(level, kind)` of the highest need (ties by the stress that is worse to be without), or `(0, "none")`."""
    if not needs:
        return 0, "none"
    return max(needs, key=lambda need: (need[0], _KIND_RANK.index(need[1])))


def _body(world, limits, cold, hot, density, scale):
    """`(required, recommended)` Body Protection needs, each a list of `(level, kind)`."""
    required, recommended = [], []
    pressure = world.pressure_kpa
    low = limit("suit_required_kpa") * scale
    if pressure < low:
        required.append((4, "pressure"))
    elif pressure < limit("suit_recommended_kpa") * scale:
        recommended.append((4, "pressure"))
    if density > limit("density_max_g_l"):
        required.append((5, "pressure"))
    elif density > limit("density_rec_g_l"):
        recommended.append((5, "pressure"))
    if hot > limit("vehicle_heat_c") or pressure > limit("vehicle_pressure_kpa") or _gravity_g(world) > limit("vehicle_gravity_g"):
        required.append((6, "gravity" if _gravity_g(world) > limit("vehicle_gravity_g") else "heat"))
    if pressure > limits.pressure[5]:
        required.append((2, "diving"))
    cold_suit, cold_clothing, comfort_low = limits.temperature[1], limits.temperature[2], limits.temperature[3]
    if cold < cold_suit:
        required.append((3, "cold"))
    elif cold < cold_clothing:
        required.append((2, "cold"))
        if cold < limit("cold_heated_rec_c"):
            recommended.append((3, "cold"))
    elif cold < comfort_low:
        required.append((1, "cold"))
    heat_garment, comfort_high = limits.temperature[5], limits.temperature[4]
    if hot > limit("armoured_heat_c"):
        required.append((5, "heat"))
    elif hot > limit("heat_shielded_c"):
        required.append((3, "heat"))
    elif hot > heat_garment:
        required.append((2, "heat"))
    elif hot > comfort_high:
        required.append((1, "heat"))
    if -20.0 <= hot <= 50.0:
        bulb = habitability.wet_bulb_c(hot, world.relative_humidity)
        if bulb >= limits.wet_bulb:
            required.append((1, "heat"))
        elif bulb >= limits.wet_bulb - (limit("wet_bulb_c") - limit("wet_bulb_rec_c")):
            recommended.append((1, "heat"))
    if _corrosive(world):
        required.append((3, "chemical"))
    if world.uv_index > limit("uv_suit_index"):
        required.append((2, "uv"))
    elif world.uv_index > limit("uv_cover_index"):
        required.append((1, "uv"))
    return required, recommended


def body_at(world, temperature_c, limits=habitability.HUMAN):
    """`(level, kind)` of the Body Protection the ladder lists for `world` at one surface temperature, `temperature_c` (the cold
    and hot ends of a day outside, GEN.768): the world's own pressure, air, gravity and ultraviolet, a steady temperature."""
    probe = replace(world, temperature_c=temperature_c, flux_variation=0.0)
    density = best_mix_density_g_l(probe.pressure_kpa, temperature_c)
    scale = limits.pressure[2] / limit("mask_min_kpa")
    required, _ = _body(probe, limits, temperature_c, temperature_c, density, scale)
    return _top(required)


def _gravity_g(world):
    return world.gravity_ms2 / constants.EARTH_GRAVITY


def _corrosive(world):
    """A corrosive gas past its IDLH (HF, HCl, NH3, Cl2) or a surface liquid below pH 2 or above 12."""
    if world.ph is not None and not limit("corrosive_ph_low") <= world.ph <= limit("corrosive_ph_high"):
        return True
    return any(row.skin and world.partial_kpa(row.gas) > row.idlh_kpa for row in TOXIC.values())


# --- Filter -------------------------------------------------------------------------------

def _gas_hazard(world):
    """`(worst ratio to a long-term limit, any gas past its IDLH)` over the table's gases."""
    worst, over = 0.0, False
    for row in TOXIC.values():
        partial = world.partial_kpa(row.gas)
        if partial > row.chronic_kpa:
            worst = max(worst, partial / row.chronic_kpa)
        if partial > row.idlh_kpa:
            over = True
    return worst, over


def _dust(world):
    """The dust part, 0 none, 1 recommended, 2 required: a body under `dust_pressure_kpa` has regolith dust that wrecks open
    seals (required); a body with air has the airborne dust of its stored optical depth (`World.dust_mg_m3`, GEN.658): recommended
    from the long-term limit, required from the short-term ceiling."""
    if world.pressure_kpa < limit("dust_pressure_kpa"):
        return 2
    loading = world.dust_mg_m3
    if not loading:
        return 0
    return 2 if loading >= limit("dust_required_mg_m3") else 1 if loading >= limit("dust_recommended_mg_m3") else 0


def _filters(world):
    """`(dust, gas filter grade, CO2 scrubber, over IDLH or past a filter)`; the scrubber part is 0 or 2 (required), the
    dust part 0, 1 (recommended) or 2 (required)."""
    dust = _dust(world)
    chronic, acute = habitability.GAS_LIMITS_KPA["CO2"]
    carbon_dioxide = world.partial_kpa("CO2")
    scrub = 2 if chronic < carbon_dioxide <= acute else 0
    worst, over = _gas_hazard(world)
    grade = 0 if worst <= 0 else 1 if worst <= HAZARD_RATIOS[0] else 2 if worst <= HAZARD_RATIOS[1] else 3
    return dust, grade, scrub, over or worst > HAZARD_RATIOS[2]


# --- Supply -------------------------------------------------------------------------------

def _supply(world, limits, cold, hot, density, scale, sealed_gas, filtering, body_level):
    """`(required, recommended)` Air Supply level."""
    pressure = world.pressure_kpa
    oxygen = habitability.inspired_o2_kpa(world, limits)
    low = limit("suit_required_kpa") * scale
    if density > limit("density_max_g_l") or body_level >= 5:
        return 4, 0
    if pressure > limits.pressure[5]:
        return 3, 0
    level = 0
    if pressure < low or body_level == 4:
        level = 2
    if oxygen > limits.o2_acute or world.partial_kpa("CO2") > habitability.GAS_LIMITS_KPA["CO2"][1] or sealed_gas:
        level = 2
    if cold < limits.temperature[1] or hot > limit("heat_shielded_c"):
        level = 2
    recommended = 0
    if level < 2:
        thin = oxygen < limits.o2_min and pressure >= low
        rich = oxygen > limits.o2_chronic
        dense = pressure > limits.pressure[4]
        if thin or rich or dense:
            level = 2 if filtering else 1
        elif oxygen < limit("inspired_o2_rec_kpa") * (limits.o2_min / limit("inspired_o2_min_kpa")) and pressure >= low:
            recommended = 1
    return level, recommended


# --- Radiation ----------------------------------------------------------------------------

def _dose_parts(world):
    """`(cosmic rays, stellar particles, ground)`, Sv/yr; the whole dose is cosmic rays when the parts are not known."""
    if world.dose_gcr_sv_yr is None or world.dose_sep_sv_yr is None or world.dose_ground_sv_yr is None:
        return world.dose_sv_yr, 0.0, 0.0
    return world.dose_gcr_sv_yr, world.dose_sep_sv_yr, world.dose_ground_sv_yr


class _Shield:
    """The dose under extra shielding: the cosmic-ray and storm parts scale by the dose curves of `radiation` at the
    world's air column plus the shielding; the ground's part does not (burial does not remove it)."""

    def __init__(self, world):
        self.gcr, self.sep, self.ground = _dose_parts(world)
        self.belt = world.dose_belt_sv_yr or 0.0
        pressure_pa = world.pressure_kpa * 1000.0
        self.column = radiation.column_g_cm2(pressure_pa, _gravity_g(world))
        self.cutoff = world.cutoff_gv
        self.gcr_base = radiation.gcr_dose_msv_yr(self.column, self.cutoff)
        self.sep_base = radiation.sep_dose_msv_yr(1.0, self.column, self.cutoff)
        self._attenuated_at = {}

    def dose(self, layers):
        """The dose, Sv/yr, with body shares and g/cm2 of shielding `layers`."""
        total = 0.0
        for share, grams in layers:
            total += share * self._attenuated(grams)
        return total + self.ground

    def _attenuated(self, grams):
        """The cosmic-ray and storm dose behind `grams` g/cm2 of shielding; the plans ask the same few depths again and again
        (every suit step and shelter for each limit), so each is worked out once (TEST.392)."""
        found = self._attenuated_at.get(grams)
        if found is not None:
            return found
        if grams <= 0:
            found = self.gcr + self.sep + self.belt
        else:
            column = self.column + grams
            gcr = self.gcr * radiation.gcr_dose_msv_yr(column, self.cutoff) / self.gcr_base if self.gcr_base > 0 else 0.0
            sep = self.sep * radiation.sep_dose_msv_yr(1.0, column, self.cutoff) / self.sep_base if self.sep_base > 0 else 0.0
            found = gcr + sep + self.belt * radiation_belts.air_cut(grams)
        self._attenuated_at[grams] = found
        return found

    def column_for(self, limit_sv):
        """The g/cm2 of shielding all round that brings the dose to `limit_sv`, or `None` when the ground's own dose
        is already past it or no depth up to `SHELTER_SEARCH_MAX_G_CM2` does."""
        if self.ground >= limit_sv:
            return None
        if self.dose(((1.0, SHELTER_SEARCH_MAX_G_CM2),)) > limit_sv:
            return None
        low, high = 0.0, SHELTER_SEARCH_MAX_G_CM2
        for _ in range(40):
            mid = 0.5 * (low + high)
            if self.dose(((1.0, mid),)) > limit_sv:
                low = mid
            else:
                high = mid
        return high


def _radiation_plan(shield, limit_sv):
    """`(suit step 0 to 3, shelter 0 to 2)` that bring the year's dose to `limit_sv` or less. A suit alone when one
    does; else a shelter (the first that does, the habitat when none does) and the lowest suit step that keeps a
    working year outside (`outside_hours`) under the limit, the hardened suit when none does."""
    steps = RADIATION_STEPS
    bare = shield.dose(((1.0, 0.0),))
    if bare <= limit_sv:
        return 0, 0
    for index in range(3):
        if shield.dose(steps[index].layers) <= limit_sv:
            return index + 1, 0
    shelter = 1 if shield.dose(steps[3].layers) <= limit_sv else 2
    year = limit("work_hours_outside_yr") / HOURS_PER_YEAR
    if bare * year <= limit_sv:
        return 0, shelter
    for index in range(3):
        if shield.dose(steps[index].layers) * year <= limit_sv:
            return index + 1, shelter
    return 3, shelter


def _outside_hours(shield, suit, limit_sv):
    """The hours outside a year the worn `suit` step allows before the dose reaches `limit_sv`, or `None` when a
    working year outside is within it."""
    layers = RADIATION_STEPS[suit - 1].layers if suit else ((1.0, 0.0),)
    yearly = shield.dose(layers)
    if yearly <= 0:
        return None
    hours = limit_sv / yearly * HOURS_PER_YEAR
    return hours if hours < limit("work_hours_outside_yr") else None


def _radiation(world, limits):
    """`(req suit, rec suit, req shelter, rec shelter, shelter column, outside hours)`: the required step is the
    lowest that brings the year's dose under the Yellow edge (out of Red), the recommended the lowest that brings it
    under the Green edge; a recommended step is kept only above the required one."""
    dose = world.dose_sv_yr
    green, yellow = limits.dose[5], limits.dose[6]
    if dose <= green:
        return 0, 0, 0, 0, None, None
    shield = _Shield(world)
    required_suit, required_shelter = _radiation_plan(shield, yellow) if dose > yellow else (0, 0)
    recommended_suit, recommended_shelter = _radiation_plan(shield, green)
    if recommended_shelter and not required_shelter:
        # No ordinary garment brings the year under the Green edge (cosmic rays): in the game a hardened suit is rated
        # against them, so it is required and the habitat is the recommended step above it (Boss, 2026-10-23).
        required_suit = 3
    column = None
    if required_shelter or recommended_shelter:
        column = shield.column_for(yellow if required_shelter else green)
    hours = _outside_hours(shield, required_suit, yellow) if required_shelter else None
    return (required_suit, recommended_suit if recommended_suit > required_suit else 0, required_shelter,
            recommended_shelter if recommended_shelter > required_shelter else 0, column, hours)


# --- Putting it together ------------------------------------------------------------------

def compute(world, limits=habitability.HUMAN):
    """
    The stored levels (`FIELDS`) of a world, for a visitor whose limits are `limits` (`habitability.bands`; a species'
    own for its fit).

    Cross-slot rules, in this order: a pressure suit carries its own loop (Supply 2 at least); an armoured suit or a
    vehicle carries its own isolated air (Supply 4); an open mask that is also wanted with a filter becomes a sealed
    supply; with a sealed supply the gas filter and the CO2 scrubber drop out (the person breathes own gas) and the
    dust part stays as seals.
    """
    cold, hot = habitability.temperature_extremes_c(world)
    density = best_mix_density_g_l(world.pressure_kpa, world.temperature_c)
    scale = limits.pressure[2] / limit("mask_min_kpa")
    required_body, recommended_body = _body(world, limits, cold, hot, density, scale)
    body, kind = _top(required_body)
    body_rec, kind_rec = _top(recommended_body)
    if body_rec <= body:
        body_rec, kind_rec = 0, "none"
    dust, gas_grade, scrub, sealed_gas = _filters(world)
    supply, supply_rec = _supply(world, limits, cold, hot, density, scale, sealed_gas, bool(gas_grade or scrub), body)
    if body >= 5:
        supply = 4
    elif body == 4:
        supply = max(supply, 2)
    if supply >= 2:
        gas_grade, scrub = 0, 0
    gravity = _gravity_g(world)
    mobility = 2 if limit("frame_gravity_g") <= gravity <= limit("vehicle_gravity_g") else (
        1 if limit("frame_gravity_rec_g") <= gravity < limit("frame_gravity_g") else 0)
    suit, suit_rec, shelter, shelter_rec, column, hours = _radiation(world, limits)
    return {
        "eq_body": body, "eq_body_kind": BODY_KINDS.index(kind), "eq_body_rec": body_rec,
        "eq_body_rec_kind": BODY_KINDS.index(kind_rec), "eq_rad_suit": suit, "eq_rad_suit_rec": suit_rec,
        "eq_mobility": mobility, "eq_dust": dust, "eq_gas_filter": gas_grade, "eq_co2_scrub": scrub,
        "eq_supply": supply, "eq_supply_rec": supply_rec if supply_rec > supply else 0, "eq_shelter": shelter,
        "eq_shelter_rec": shelter_rec, "shelter_column_g_cm2": column, "outside_hours": hours,
    }


# --- Search ticks -------------------------------------------------------------------------

TICKS = (
    ("clothing", "Weather clothing", "{a}.eq_body = 1"), ("thermal_suit", "Thermal suit", "{a}.eq_body = 2"),
    ("sealed_suit", "Sealed suit", "{a}.eq_body = 3"), ("pressure_suit", "Pressure suit", "{a}.eq_body = 4"),
    ("armoured_suit", "Armoured suit", "{a}.eq_body = 5"), ("vehicle_only", "Vehicle only", "{a}.eq_body = 6"),
    ("powered_frame", "Powered frame", "{a}.eq_mobility = 2"),
    ("shielded_garment", "Shielded garment", "{a}.eq_rad_suit = 1"), ("storm_vest", "Storm vest", "{a}.eq_rad_suit = 2"),
    ("hardened_suit", "Hardened suit", "{a}.eq_rad_suit = 3"),
    ("dust", "Dust filter or seals", "{a}.eq_dust = 2"), ("gas_filter", "Gas filter", "{a}.eq_gas_filter > 0"),
    ("co2_scrubber", "CO2 scrubber", "{a}.eq_co2_scrub = 2"),
    ("gas_mix_mask", "Gas-mix mask", "{a}.eq_supply = 1"), ("own_air", "Self-contained air", "{a}.eq_supply = 2"),
    ("mixed_gas", "Mixed-gas rig", "{a}.eq_supply = 3"), ("isolated_air", "Isolated air", "{a}.eq_supply = 4"),
    ("shielded_shelter", "Shielded shelter", "{a}.eq_shelter = 1"),
    ("hardened_habitat", "Hardened habitat", "{a}.eq_shelter = 2"),
)
"""tuple: The search page's "needs" ticks: `(value, label, SQL condition)`, the condition on a `planets` or `moons` row
with the table alias in `{a}`. A tick is a required item (UX.153 adds the recommended ones)."""

TICK_VALUES = frozenset(value for value, _label, _sql in TICKS)
"""frozenset: The tick values."""


def tick_condition(alias, values):
    """The SQL condition (a bare `OR` of the ticks in `values`, parenthesised; unknown values are dropped) on the table
    `alias`, or `None` when none of `values` is a tick."""
    conditions = [sql.format(a=alias) for value, _label, sql in TICKS if value in values]
    return "(" + " OR ".join(conditions) + ")" if conditions else None


def tick_counts_select(alias):
    """The select list that counts the bodies of table `alias` each tick needs (`t0`, `t1`, ... in `TICKS` order)."""
    return ", ".join(f"COALESCE(SUM({sql.format(a=alias)}), 0) AS t{i}" for i, (_v, _l, sql) in enumerate(TICKS))


# --- Names, from the stored levels alone ----------------------------------------------------

def _get(levels, key):
    value = levels.get(key)
    return value or 0


def scored(levels):
    """True when `levels` hold a computed line (a gas giant or an unscored body has none)."""
    return levels.get("eq_body") is not None


def _body_word(level, kind_number):
    kind = BODY_KINDS[kind_number] if 0 <= kind_number < len(BODY_KINDS) else "none"
    return BODY_NAMES.get((level, kind), BODY_DEFAULT.get(level, "suit"))


def _list(levels, severity):
    """The `Item`s of one severity, in Boss's order: body (with its frame and radiation suit), filter, supply,
    shelter."""
    required = severity == REQUIRED
    suffix = "" if required else "_rec"
    body = _get(levels, f"eq_body{suffix}")
    kind = _get(levels, f"eq_body{suffix}_kind")
    suit = _get(levels, f"eq_rad_suit{suffix}")
    shelter = _get(levels, f"eq_shelter{suffix}")
    supply = _get(levels, f"eq_supply{suffix}")
    out = []
    worn = body or (0 if required else _get(levels, "eq_body"))  # a recommended hardened suit hardens the required one
    worn_kind = kind if body else _get(levels, "eq_body_kind")
    word = _body_word(worn, worn_kind) if worn else None
    domain = KIND_DOMAIN.get(BODY_KINDS[worn_kind], "Pressure") if worn_kind else "Pressure"
    domains = (domain,)
    if worn >= 4:  # a suit that holds a pressure round the wearer also keeps out the cold, the heat and the air
        domains = tuple(dict.fromkeys((domain, "Pressure", "Temperature", "Air")))
    elif worn == 3:
        domains = (domain, "Air")
    if suit == 3:
        hardened = worn >= 3 and worn < 6
        word = "hardened " + word if hardened else "hardened suit"
        domains = (*domains, "Radiation") if hardened else ("Radiation",)
    if word and (body or suit == 3):
        out.append(_item("body", word, severity, *domains))
    elif suit == 3:
        out.append(_item("body", "hardened suit", severity, "Radiation"))
    if _get(levels, "eq_mobility") == (2 if required else 1):
        out.append(_item("body", "powered frame", severity, "Pressure"))
    if suit in (1, 2):
        out.append(_item("body", RAD_SUIT_NAMES[suit], severity, "Radiation"))
    if required:
        sealed = _get(levels, "eq_supply") >= 2 or _get(levels, "eq_body") >= 4
        if _get(levels, "eq_dust") == 2:
            out.append(_item("filter", "dust-tight seals" if sealed else "dust filter", severity, "Air"))
        grade = _get(levels, "eq_gas_filter")
        if grade:
            out.append(_item("filter", GAS_FILTER_NAMES[grade], severity, "Air"))
        if _get(levels, "eq_co2_scrub"):
            out.append(_item("filter", "CO2 scrubber", severity, "Air"))
    elif _get(levels, "eq_dust") == 1:
        sealed = _get(levels, "eq_supply") >= 2 or _get(levels, "eq_body") >= 4
        out.append(_item("filter", "dust-tight seals" if sealed else "dust filter", severity, "Air"))
    carried = _get(levels, "eq_body")
    if supply and not (carried >= 4 and supply <= 2) and not (carried >= 5 and supply == 4):
        out.append(_item("supply", SUPPLY_NAMES[supply], severity,
                         *(("Air", "Pressure", "Temperature") if supply < 3 else ("Pressure",))))
    if shelter:
        out.append(_item("radiation", SHELTER_NAMES[shelter], severity, "Radiation"))
    return out


def items(levels):
    """All the line's `Item`s, required first (an item of the same name is not listed twice, the required one wins)."""
    if not scored(levels):
        return []
    out = _list(levels, REQUIRED)
    names = {item.name for item in out}
    return out + [item for item in _list(levels, RECOMMENDED) if item.name not in names]


def required_names(levels):
    """The names of the required items, in order."""
    return [item.name for item in items(levels) if item.severity == REQUIRED]


def recommended_names(levels):
    """The names of the recommended items, in order."""
    return [item.name for item in items(levels) if item.severity == RECOMMENDED]


def line(levels):
    """The equipment line as text: "Required: a + b. Recommended: c." A list with nothing in it is left out, and a world
    that needs nothing gives the empty string (never "none", Boss, decision E)."""
    parts = []
    required, recommended = required_names(levels), recommended_names(levels)
    if required:
        parts.append("Required: " + " + ".join(required) + ".")
    if recommended:
        parts.append("Recommended: " + " + ".join(recommended) + ".")
    return " ".join(parts)


def lines(levels):
    """`(required line, recommended line)`: "Required: a + b." and "Recommended: c.". A list with nothing in it is left off
    (the empty string, never "none"; Boss, 2026-10-24), so a world that needs nothing has neither."""
    required, recommended = required_names(levels), recommended_names(levels)
    return (("Required: " + " + ".join(required) + ".") if required else "",
            ("Recommended: " + " + ".join(recommended) + ".") if recommended else "")


def _tokens(levels, suffix):
    body = _get(levels, f"eq_body{suffix}")
    out = []
    if body:
        out.append(f"B{body}{KIND_LETTER[BODY_KINDS[_get(levels, f'eq_body{suffix}_kind')]]}")
    if suffix == "":
        parts = (_get(levels, "eq_dust") == 2, _get(levels, "eq_gas_filter"), _get(levels, "eq_co2_scrub") > 0)
        if any(parts):
            out.append(f"F{int(parts[0])}{parts[1]}{int(parts[2])}")
    elif _get(levels, "eq_dust") == 1:
        out.append("F100")
    supply = _get(levels, f"eq_supply{suffix}")
    if supply:
        out.append(f"S{supply}")
    if _get(levels, "eq_mobility") == (2 if suffix == "" else 1):
        out.append("M")
    suit = _get(levels, f"eq_rad_suit{suffix}")
    if suit:
        out.append(f"H{suit}")
    shelter = _get(levels, f"eq_shelter{suffix}")
    if shelter:
        out.append(f"R{shelter}")
    return out


def code(levels):
    """The short code: `B<level><kind> F<dust><gas><co2> S<level> M H<suit step> R<shelter>` for the required items, then
    the recommended ones in brackets (`[S1 R2]`); empty for a world that needs nothing. Mars: `B4p F100 S2 [R2]`."""
    if not scored(levels):
        return ""
    required, recommended = _tokens(levels, ""), _tokens(levels, "_rec")
    text = " ".join(required)
    if recommended:
        text = (text + " " if text else "") + "[" + " ".join(recommended) + "]"
    return text


def summary(levels):
    """The `equipment` object the API serves (`None` for an unscored body): `required` and `recommended` names, the
    `code`, the `items` (slot, item, severity, domain), the `line` and the shelter's depth and the hours outside (a year under the dose limit, and a day with what limits it, GEN.768)."""
    if not scored(levels):
        return None
    return {
        "required": required_names(levels), "recommended": recommended_names(levels), "code": code(levels),
        "items": [{"slot": item.slot, "item": item.name, "severity": item.severity, "domain": item.domain,
                   "domains": list(item.domains)}
                  for item in items(levels)],
        "line": line(levels), "required_line": lines(levels)[0], "recommended_line": lines(levels)[1],
        "shelter_column_g_cm2": levels.get("shelter_column_g_cm2"),
        "outside_hours": levels.get("outside_hours"), "outside_hours_day": levels.get("outside_hours_day"),
        "outside_hours_limit": levels.get("outside_hours_limit"), "regions": regional_temperature.summary(levels),
    }


def tech_needed(levels, table):
    """The tech level a species needs to make the required gear: the highest of `table`
    (`tuning.SPECIES_EQUIPMENT_TECH`, item name to tech level) over the required items; a "hardened" suit needs the
    larger of the hardened suit's level and its own. 0 for a world that needs nothing."""
    return max([0.0] + [item_tech(name, table) for name in required_names(levels)])


def item_tech(name, table):
    """The tech level `table` (`tuning.SPECIES_EQUIPMENT_TECH`) gives the item `name`: a "hardened" suit needs the larger
    of the hardened suit's level and its own; 0 for a name the table does not know."""
    level = table.get(name)
    if level is None and name.startswith("hardened "):
        level = max(table["hardened suit"], table.get(name[len("hardened "):], 0.0))
    return level or 0.0
