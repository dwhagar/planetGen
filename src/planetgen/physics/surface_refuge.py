# planetgen/physics/surface_refuge.py

"""
Surface refuges for oxygen-making life (GEN.1723, GEN.1722 research)
======================================================================

Boss's rule (8:38 PM EST 11-14): *life will always find a way*. Oxygenic photosynthesis is present on every world that has **light and any liquid water at its
surface**, an open sea, a brine seep, a seasonal melt, a spring or a hot pool, and not only where the hydrosphere is a wide sunlit ocean. Under ice and in an aquifer there is no
light, so those places count for nothing. The sourced limits and the size of the production are in `reports/life-refuges-research-gen1722-2026-11-14.md` (Lane 03).

**The gate (RW1).** A place is a refuge when it has liquid water, lit, `MIN_HOURS_YR` hours a year or more, at a temperature of `LIQUID_MIN_K` to `LIQUID_MAX_K` (the brine's eutectic to
the limit of cyanobacterial mats in hot springs). `lit_liquid_hours` is the most hours any place of the world has. An open sea (`surface ocean`, `hycean`) has the whole year, as it
always did, so a world with one is scored exactly as before. The rest is read from the temperature of every place (`temperature_coverage.series_of_values`, a stored body is
worked once and kept):

* a `brine seep`: its liquid (`ocean_fraction`) where the temperature is in the range;
* an `ice-covered ocean` or an `ice` world: the meltwater where a place rises above 273.15 K (the ice cover of the place for an `ice` world: a melt hole, a cryoconite pond);
* the sea under the lid of an `ice-covered ocean` (GEN.1724 research, LA4, LA5; GEN.1735): each cell of the world wears the lid its own mean temperature, the swing of its year and the starlight its
  surface absorbs give (`ice_growth.lid_of`: the ice a freezing season grows under Stefan's law, the snow on it, and for a year that never thaws the balance of the ice that is worn from the top and grown from below), the lid lets
  `ice_light.transmission` (GEN.1726) of the light through, and the sea (`ocean_fraction`) makes the ocean's own production in proportion to it, at 273.15 K, for the share of the year the cell is
  frozen, when what gets through is over the light floor (`LIGHT_FLOOR_W_M2`, RW2). A lid of a few metres counts (the lake mats under 3 to 5 m of ice, most of the Arctic Ocean's production is under ice), a lid of a kilometre
  does not; a cell whose year never freezes has no lid, and its open water is counted where it is wet and lit (the melting steps above);
* a world of `SURFACE_STATES` with an aquifer (`aquifers`): the springs and seeps where the water table meets the ground, `SEEP_LAND_SHARE` of the land at a table at the surface and
  falling as `exp(-depth / SEEP_DEPTH_KM)` with the depth of its top ([R5271, R5272]: nearly all the groundwater that reaches the surface comes from a table within a metre), in
  the range. A table deeper than a few metres gives none (the springs where relief cuts it are left out: the minimum reasonable case); the water of a deep aquifer reaches no light. Once the world's life has
  reached mosses and the first land plants (`native_life` `macroscopic` or `biosphere`, LA7), the roots reach the table: `ROOT_SEEP_SHARE` of the land (22 to 32 percent is under shallow groundwater influence) and the table's depth is
  measured on the rooting depth of a desert plant, 9.5 m [R5274], in place of the metre. The oxygen of the rise is worked before that stage (the stage needs the oxygen), so the widening is in the scores and the source of a world already there.

**The production (RW3).** `lit_area` is the share of the world's surface that is wet and lit, each place weighed by `production_share(hours)`, the production per area against an
ocean's: `(hours / 8760) ** k * cover`, with `cover` rising from `COVER_MIN` at the 75-hour limit to `COVER_MAX` at `COVER_FULL_HOURS` (a log scale; the hypolith and soil crust
covers of [R5269] run 0.1 to 28 percent, and the minimum reasonable one is taken) and `k` 1 above `K_ONE_HOURS` rising to `K_MAX` at the limit (the five orders of
magnitude of the Warren-Rhodes gradient over a tenfold fall of the water, taken only as its direction). `biogenic_oxygen.productivity` turns that into the source `S` with
Earth's own burial per area of ocean: **no new burial constant**. The result is small on purpose: about 1e-2 of Earth's `S` for a seep world with a year of warm light, 1e-3 or less for a melt hole or a cold edge, and 0.14 at the very most (the cover cap of 0.1 over a whole sea that melts all year), so it never beats the reductant flux `R` (0.3 of Earth's `S`) on a world whose interior is Earth's. If dry worlds are to hold oxygen, the lever is `R` (a
more oxidized mantle), not this number.

Limits. The series is a circular orbit with an evenly spread sea, and "lit" is the mean share of the day the star is up (`day_share_mean`, 0.5 when not stored) applied to every
step, not each step's own light; the aquifer's supply of water is not scaled by its thickness (the research's `G / G_Earth` needs a head the generator does not keep), so every
aquifer with its table at the surface counts as Earth's. Tiers and tweak levels are in the stand-in registry (RW1 to RW3).
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import aquifers, constants, ice_growth, ice_light, temperature_coverage

YEAR_HOURS = 8760.0
"""float: The hours of a year, 365 days of 24 hours (a definition: source tier T1, tweak level L0)."""
LIQUID_MIN_K = 252.0
"""float: The coldest a place's liquid water is counted at, K: the eutectic of a chloride brine (RW1, tweak level L0). The oxygen makers are active to about 253 K in the field (LA1; the production taper below)."""
RATE_FULL_K = 263.0
RATE_ZERO_K = 253.0
"""float: The production factor of a wet place is 1 at `RATE_FULL_K` and above and falls linearly to 0 at `RATE_ZERO_K`, K: Antarctic lichens photosynthesise to -20 C in the field and the rate is small below -10 C ([R5275, R5276] and the field limit of Antarctic lichens (-20 C), read; reports/ice-light-and-lichen-cold-edge-2026-11-14.md). LA1 (the edge, T1, L0) and LA1b (the straight line between the two points, T6, L1)."""
LIQUID_MAX_K = 346.0
"""float: The hottest, K: the upper limit of cyanobacterial mats in hot springs, about 73 C (Castenholz, from memory). RW1, tweak level L0."""
MELT_K = 273.15
"""float: Meltwater needs a place above this temperature, K (ice has no liquid in it below it, a definition)."""
MIN_HOURS_YR = 75.0
"""float: The least hours of lit liquid water in a year a refuge has ([R5269]: the hypolith mats of the Mojave and Namib live on about 75 h a year of liquid water, from memory). RW1, educated guess from physics, tweak level L0."""
SEEP_LAND_SHARE = 0.15
"""float: The share of the land a groundwater-fed surface feature (a spring, seep, oasis, playa or wetland) covers on Earth ([R5271], Science). RW1, tweak level L2."""
SEEP_DEPTH_KM = 0.001
"""float: The depth, km, of a water table over which the share of seeps falls by 1/e: 1 m ([R5272]: nearly all the groundwater that evaporates comes from a table within 0.5 m; the share ends at 30 to 100 cm). RW1, tweak level L2."""
COVER_MIN = 0.001
COVER_MAX = 0.1
"""float: The cover of the life on the wet lit area, 0.001 at the dry limit up to 0.1 on a wet crust (the low end of the 0.1 to 28 percent of [R5269]: the minimum reasonable one). RW3, educated guess from physics, tweak level L2."""
COVER_FULL_HOURS = 2000.0
"""float: The hours a year at which a place has the full cover `COVER_MAX` (RW3, tweak level L2)."""
K_ONE_HOURS = 500.0
K_MAX = 5.0
"""float: The hours a year above which production goes in proportion to them (`k` 1), and the exponent it rises to at the 75-hour limit (RW3, tweak level L2)."""
EARTH_OCEAN_SHARE = 0.71
"""float: The share of its surface Earth's ocean covers (the area `biogenic_oxygen.EARTH_BURIAL_MOL_YR` is the burial of; the same 0.71 as `biogenic_oxygen.EARTH_OCEAN_FRACTION`). Source tier T1, tweak level L0."""
LIGHT_FLOOR_W_M2 = 2.0e-3
"""float: The least usable light, W/m2, under which no photosynthesis is counted: 0.01 umol photons m-2 s-1 ([R5273], RW2), about 2e-3 W/m2 of photosynthetically active light. Tier T1, tweak level L0."""
UNDER_ICE_K = 273.15
"""float: The temperature of the water under a lid of ice, K (the melting point; a definition: source tier T1, tweak level L0)."""
ROOTED = ("macroscopic", "biosphere")
"""tuple: The native-life stages (`habitability_world.native_life`) with mosses and the first land plants (`habitability.LIFE_STAGES`): the roots reach the water table (LA7)."""
ROOT_SEEP_SHARE = 0.27
"""float: The share of the land a rooted plant draws groundwater from once the stage is reached: the middle of the 22 to 32 percent of Earth's land under shallow groundwater influence ([R5271], read). LA7, tier T1 to T6, tweak level L1."""
RANGE_SHARE = 0.95
"""float: The share of the wet lit time of the refuge places that the temperature range of `Refuge.coldest_k` to `Refuge.hottest_k` holds (GEN.1808): the life lives where most of the water is, so the range leaves out the coldest and hottest
`(1 - RANGE_SHARE) / 2` of it each. A design choice of the adaptation record, not a measured one. [Stand-in 6.318 AE2, tier T6, tweak L0.]"""
ROOT_DEPTH_KM = 0.0095
"""float: The depth, km, a water table falls by 1/e in the share of the land a rooted plant reaches: the mean maximum rooting depth of desert plants, 9.5 m (the record is 68 m; [R5274], read at second hand). LA7, tier T1, tweak level L0."""
DEFAULT_ALBEDO = 0.3
"""float: The Bond albedo taken for a body that stores none, when the starlight is worked from its star and distance (Earth's, a stand-in: source tier T6, tweak level L0)."""
DEFAULT_LIT_SHARE = 0.5
"""float: The share of the time the star is up when `day_share_mean` is not stored (a half: stand-in, source tier T6, tweak level L0)."""


@dataclass(frozen=True)
class Refuge:
    """
    The surface refuges of a world.

    Attributes:
        hours (float): The most lit liquid hours a year any place has (`YEAR_HOURS` for an open sea, 0 for none).
        area (float): The share of the world's surface that is wet and lit, weighed by each place's production against an ocean's (`production_share`); 0 for an open sea
            (its source is the ocean's own).
        temperature_k (float or None): The mean temperature of the wet lit steps, K, weighed by area; `None` for none.
        plants (float): The share of the world's surface that rooted land plants hold wet and lit (GEN.1729): the groundwater-fed places once the life has reached the
            rooted stage (`ROOTED`), on the land, each weighed by the share of the year it is wet and lit against Earth's lit half (production in proportion to the time, with no cover factor: a
            plant stand fills its wet land, as an ocean fills its sea; Earth's land is wet the whole year); 0 before that stage and for a world with no water table the roots reach.
        plants_cover (float): The part of `area` that those same places make under the cover law of `production_share` (what `plants` replaces in the oxygen source).
        coldest_k (float or None): The cold end of the temperatures the wet lit steps of the places counted have, K (GEN.1808: what the life's cold limit must reach): the temperature under which `1 - RANGE_SHARE` over two
            of the wet lit time of the places lies, each step weighed by its place's share of the surface; `None` for an open sea, which holds no such place, and for none.
        hottest_k (float or None): The hot end of the same range, K (the hot limit it must reach).
        under_ice_w_m2 (float or None): The least usable light, W/m2, of the cells of the sea under its lid that count (the light floor it must reach); `None` when no such cell counts.
        habitats (tuple): The ids of the habitats of `life_adaptations.HABITATS` the places counted belong to (`open_water`, `brine`, `melt_film`, `seep`, `water_under_ice`).
    """

    hours: float
    area: float
    temperature_k: float | None
    plants: float = 0.0
    plants_cover: float = 0.0
    coldest_k: float | None = None
    hottest_k: float | None = None
    under_ice_w_m2: float | None = None
    habitats: tuple = ()


NONE = Refuge(0.0, 0.0, None)
OPEN_SEA = Refuge(YEAR_HOURS, 0.0, None, habitats=("open_water",))
OPEN_SEAS = ("surface ocean", "hycean")
"""tuple: The states of a wide sunlit sea (`biogenic_oxygen.LIVING_HYDROSPHERES`); their gate and source are the sea's."""


def production_share(hours):
    """The production per area of a place with `hours` of lit liquid water a year, against an ocean's: `(hours / 8760) ** k * cover` (RW3); 0 under `MIN_HOURS_YR`."""
    if hours < MIN_HOURS_YR:
        return 0.0
    hours = min(hours, YEAR_HOURS)
    span = math.log(COVER_FULL_HOURS / MIN_HOURS_YR)
    cover = COVER_MIN * (COVER_MAX / COVER_MIN) ** min(1.0, math.log(hours / MIN_HOURS_YR) / span)
    exponent = 1.0 if hours >= K_ONE_HOURS else 1.0 + (K_MAX - 1.0) * math.log(K_ONE_HOURS / hours) / math.log(K_ONE_HOURS / MIN_HOURS_YR)
    return (hours / YEAR_HOURS) ** exponent * cover


def _hours_array(hours):
    """`production_share` of an array of hours."""
    return np.vectorize(production_share, otypes=[float])(hours) if hours.size else hours


def refuge_of(values):
    """
    The `Refuge` of a rocky body from its column-named `values` (`habitability_world.values_from_body`, or a stored row); `NONE` for a world with no liquid water lit at its surface and
    `OPEN_SEA` for one with a wide sunlit sea.
    """
    state = values.get("hydrosphere")
    if state in OPEN_SEAS:
        return OPEN_SEA
    mean = values.get("surface_temperature_k")
    if not mean:
        return NONE
    aquifer = bool(values.get("aquifer_thickness_km")) and state in aquifers.SURFACE_STATES
    if state not in ("brine seep", "ice-covered ocean", "ice") and not aquifer:
        return NONE
    series = temperature_coverage.series_of_values(values) or temperature_coverage.still_series(mean)
    temperature = series.temperature_k
    cover = getattr(series, "ice_cover", None)
    ice = np.zeros(temperature.shape[:2]) if cover is None else np.asarray(cover, dtype=float)
    steps = temperature.shape[2]
    lit = values.get("day_share_mean")
    lit = DEFAULT_LIT_SHARE if lit is None else min(max(float(lit), 0.0), 1.0)
    liquid = (temperature >= LIQUID_MIN_K) & (temperature <= LIQUID_MAX_K)
    melting = (temperature >= MELT_K) & (temperature <= LIQUID_MAX_K)
    sources = []        # (share of each cell that is wet, the steps it is wet at, whether rooted plants stand on it, the habitat it is)
    if state == "brine seep":
        sources.append((np.full(ice.shape, float(values.get("ocean_fraction") or values.get("ground_ocean_fraction") or 0.0)), liquid, False, "brine"))
    elif state == "ice-covered ocean":
        sources.append((np.ones(ice.shape), melting, False, "melt_film"))
    elif state == "ice":
        sources.append((ice, melting, False, "melt_film"))
    rooted = values.get("native_life") in ROOTED
    land = 1.0 - min(max(float(values.get("ocean_fraction") or 0.0), 0.0), 1.0)         # the share of each cell that is land, where the roots stand (GEN.1729)
    if aquifer:
        share, scale = (ROOT_SEEP_SHARE, ROOT_DEPTH_KM) if rooted else (SEEP_LAND_SHARE, SEEP_DEPTH_KM)
        seeps = share * math.exp(-(values.get("aquifer_depth_km") or 0.0) / scale)
        if seeps > 1.0e-6:
            sources.append((seeps * (1.0 - ice), liquid, rooted, "seep"))
    cells = temperature_coverage._cell_weights(series)
    rate = np.clip((temperature - RATE_ZERO_K) / (RATE_FULL_K - RATE_ZERO_K), 0.0, 1.0)         # LA1b: the oxygen makers slow from 263 K to none at 253 K
    most, area, warm, weight, plants, plants_cover = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    ends, habitats = [], []
    for share, wet, standing, habitat_id in sources:
        wet_steps = wet.sum(axis=2)
        hours = wet_steps / steps * YEAR_HOURS * lit
        has = (share > 0.0) & (hours >= MIN_HOURS_YR)
        if not has.any():
            continue
        most = max(most, float(hours[has].max()))
        counted = wet & has[..., None]
        if counted.any():       # GEN.1808: the temperatures the life of these places works at, each step weighed by its place's share of the surface
            ends.append((temperature[counted], np.broadcast_to((cells * share)[..., None], counted.shape)[counted]))
            habitats.append(habitat_id)
        effective = (wet * rate).sum(axis=2) / steps * YEAR_HOURS * lit
        made = float((cells * np.where(has, _hours_array(effective) * share, 0.0)).sum())
        area += made
        if standing:        # rooted plants stand on the places the roots reach (GEN.1729): production goes with the time they are wet and lit, not with the cover law
            plants += float((cells * np.where(has, np.minimum(effective / YEAR_HOURS / DEFAULT_LIT_SHARE, 1.0) * share * land, 0.0)).sum())
            plants_cover += made
        mean_wet = np.where(wet_steps > 0, (temperature * wet).sum(axis=2) / np.maximum(wet_steps, 1), 0.0)
        place = cells * np.where(has, share * wet_steps / steps, 0.0)
        warm += float((place * mean_wet).sum())
        weight += float(place.sum())
    under = _under_ice(values, state, lit, series)
    least_light = None
    if under is not None:           # the sea under a lid that lets light through: the ocean's own production, in proportion to the light that gets through (LA4)
        hours, produce, least_light = under
        habitats.append("water_under_ice")
        most = max(most, hours)
        area += produce
        warm += UNDER_ICE_K * produce
        weight += produce
    if most < MIN_HOURS_YR:
        return NONE
    coldest, hottest = _range_of(ends)
    return Refuge(most, area, warm / weight if weight > 0.0 else None, plants, plants_cover, coldest, hottest, least_light, tuple(dict.fromkeys(habitats)))


def _range_of(ends):
    """`(cold end, hot end)` K of the wet lit steps `ends` (`[(temperatures, weights)]`) that hold `RANGE_SHARE` of the weight, or `(None, None)` for none."""
    if not ends:
        return None, None
    temperature = np.concatenate([t for t, _ in ends])
    weight = np.concatenate([w for _, w in ends])
    if weight.sum() <= 0.0:
        return None, None
    order = np.argsort(temperature, kind="stable")
    temperature, weight = temperature[order], weight[order]
    reached = np.cumsum(weight) / weight.sum()
    tail = (1.0 - RANGE_SHARE) / 2.0
    return float(temperature[np.searchsorted(reached, tail)]), float(temperature[min(np.searchsorted(reached, 1.0 - tail), len(temperature) - 1)])


def _daily_mean(series):
    """The series with the turn of the day averaged out, `[row, column, step]` (a sea's lid answers the year, not the night): each step takes the mean of the steps of one solar day round it; a body whose day is
    not short against its year (a slow turner, a locked one) keeps its steps."""
    temperature = series.temperature_k
    steps = temperature.shape[2]
    day = getattr(series, "solar_days", None)
    window = int(round(day / series.step_days)) if day and series.step_days and math.isfinite(day) else 1
    if window < 2 or window * 4 > steps:
        return temperature
    return sum(np.roll(temperature, -k, axis=2) for k in range(window)) / window


def _absorbed_w_m2(values):
    """The mean starlight the surface absorbs, W/m2: a quarter of the starlight (`temperature_coverage.starlight_w_m2`) over the part the albedo does not return; the star's flux at the body's
    distance with a Bond albedo of 0.3 for a body with no stored albedo and greenhouse; 0 for none."""
    light, albedo = temperature_coverage.starlight_w_m2(values)
    if not light and values.get("star_luminosity_w") and values.get("star_distance_au"):
        light = values["star_luminosity_w"] / (4.0 * math.pi * (values["star_distance_au"] * constants.AU_TO_M) ** 2)
        albedo = DEFAULT_ALBEDO
    return 0.25 * light * (1.0 - albedo) if light else 0.0


def lid_of_cells(values, temperature):
    """
    The lid each cell of a sea wears and the light it lets through (GEN.1735, `ice_growth.lid_of`): `(transmission, share of the year under the lid)`, two arrays of the shape of the cells.
    Each cell's lid is the one its own mean temperature and the swing of its year (the root two of the standard deviation of the series with the day averaged out) give, with the starlight
    the surface absorbs and the heat the body's interior gives its water (`internal_heat_flux_w_m2`, the Arctic's `ice_growth.OCEAN_FLUX_W_M2` where none is stored); the light is
    `ice_light.transmission` of the ice and the snow on it. A cell whose year never freezes has no lid (`(0, 0)`: its open water is counted where it is wet and lit). The lid is the one the
    freezing season leaves, so a summer's thin ice counts as the thick one: less light, the cautious side.
    """
    mean = temperature.mean(axis=2)
    swing = math.sqrt(2.0) * temperature.std(axis=2)
    flux = values.get("internal_heat_flux_w_m2") or ice_growth.OCEAN_FLUX_W_M2
    absorbed = _absorbed_w_m2(values)
    key = np.stack([np.round(mean * 4.0) / 4.0, np.round(swing * 4.0) / 4.0], axis=-1).reshape(-1, 2)
    unique, inverse = np.unique(key, axis=0, return_inverse=True)
    through, frozen = np.zeros(len(unique)), np.zeros(len(unique))
    for i, (mean_k, swing_k) in enumerate(unique):
        share = ice_growth.frozen_share(mean_k, swing_k)
        if share > 0.0:
            lid = ice_growth.lid_of(mean_k, swing_k, absorbed, flux)
            through[i] = ice_light.transmission(lid.ice_m, lid.snow_m, kind=ice_light.SEA)
            frozen[i] = share
    inverse = np.asarray(inverse).reshape(mean.shape)
    return through[inverse], frozen[inverse]


def _under_ice(values, state, lit, series):
    """
    `(hours, production share, least light W/m2)` of the sea under the lids of its cells, or `None` (GEN.1735, LA4; the least light of a counted cell is GEN.1808's): each cell's lid is `lid_of_cells`'s; a cell counts when the light that gets through is over
    `LIGHT_FLOOR_W_M2` (the usable light at the surface times the transmission), for the share of the year it is frozen, and its production is the ocean's share of the surface times the transmitted share,
    linear below the full light. The hours are the longest any such cell is under its lid.
    """
    if state != "ice-covered ocean":
        return None
    sea = float(values.get("ocean_fraction") or values.get("ground_ocean_fraction") or 0.0)
    light = surface_light_w_m2(values)
    if sea <= 0.0 or not light:
        return None
    through, frozen = lid_of_cells(values, _daily_mean(series))
    has = (through * light >= LIGHT_FLOOR_W_M2) & (frozen > 0.0)
    if not has.any():
        return None
    cells = temperature_coverage._cell_weights(series)
    return float(frozen[has].max()) * YEAR_HOURS * lit, float((cells * np.where(has, sea * through * frozen, 0.0)).sum()), float((through * light)[has].min())


def lit_liquid_hours(values):
    """The most hours a year any place of the body has lit liquid water (the quantity the oxygenic classes' gate reads): `YEAR_HOURS` for an open sea, 0 for none."""
    return refuge_of(values).hours


def surface_light_w_m2(values):
    """The usable light a refuge gets, W/m2: the starlight of the stored albedo and greenhouse (`temperature_coverage.starlight_w_m2`) times
    `habitability_world.LIGHT_SHARE` and the body's `light_factor`; `None` where the body stores none (an open-sea or chemical-only flux is the stored one)."""
    from planetgen.physics import habitability_world

    light, _ = temperature_coverage.starlight_w_m2(values)
    if not light and values.get("star_luminosity_w") and values.get("star_distance_au"):
        light = values["star_luminosity_w"] / (4.0 * math.pi * (values["star_distance_au"] * constants.AU_TO_M) ** 2)         # a body with no stored albedo and greenhouse: the star's flux at its distance
    if not light:
        return None
    factor = values.get("light_factor")
    return habitability_world.LIGHT_SHARE * light * (1.0 if factor is None else factor)
