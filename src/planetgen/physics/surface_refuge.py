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
* the sea under the lid of an `ice-covered ocean` (GEN.1724 research, LA4, LA5): the lid lets `ice_light.transmission` (GEN.1726) of the light through, and the sea (`ocean_fraction`) makes the
  ocean's own production in proportion to it, at 273.15 K, when what gets through is over the light floor (`LIGHT_FLOOR_W_M2`, RW2). A lid of a few metres counts (the lake mats under 3 to 5 m of ice, most of
  the Arctic Ocean's production is under ice), a lid of a kilometre does not. Snow on the ice is left out;
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

from planetgen.physics import aquifers, constants, ice_light, temperature_coverage

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
ROOT_DEPTH_KM = 0.0095
"""float: The depth, km, a water table falls by 1/e in the share of the land a rooted plant reaches: the mean maximum rooting depth of desert plants, 9.5 m (the record is 68 m; [R5274], read at second hand). LA7, tier T1, tweak level L0."""
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
    """

    hours: float
    area: float
    temperature_k: float | None
    plants: float = 0.0
    plants_cover: float = 0.0


NONE = Refuge(0.0, 0.0, None)
OPEN_SEA = Refuge(YEAR_HOURS, 0.0, None)
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
    sources = []        # (share of each cell that is wet, the steps it is wet at, whether rooted plants stand on it)
    if state == "brine seep":
        sources.append((np.full(ice.shape, float(values.get("ocean_fraction") or values.get("ground_ocean_fraction") or 0.0)), liquid, False))
    elif state == "ice-covered ocean":
        sources.append((np.ones(ice.shape), melting, False))
    elif state == "ice":
        sources.append((ice, melting, False))
    rooted = values.get("native_life") in ROOTED
    land = 1.0 - min(max(float(values.get("ocean_fraction") or 0.0), 0.0), 1.0)         # the share of each cell that is land, where the roots stand (GEN.1729)
    if aquifer:
        share, scale = (ROOT_SEEP_SHARE, ROOT_DEPTH_KM) if rooted else (SEEP_LAND_SHARE, SEEP_DEPTH_KM)
        seeps = share * math.exp(-(values.get("aquifer_depth_km") or 0.0) / scale)
        if seeps > 1.0e-6:
            sources.append((seeps * (1.0 - ice), liquid, rooted))
    cells = temperature_coverage._cell_weights(series)
    rate = np.clip((temperature - RATE_ZERO_K) / (RATE_FULL_K - RATE_ZERO_K), 0.0, 1.0)         # LA1b: the oxygen makers slow from 263 K to none at 253 K
    most, area, warm, weight, plants, plants_cover = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    for share, wet, standing in sources:
        wet_steps = wet.sum(axis=2)
        hours = wet_steps / steps * YEAR_HOURS * lit
        has = (share > 0.0) & (hours >= MIN_HOURS_YR)
        if not has.any():
            continue
        most = max(most, float(hours[has].max()))
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
    under = _under_ice(values, state, lit)
    if under is not None:           # the sea under a lid that lets light through: the ocean's own production, in proportion to the light that gets through (LA4)
        hours, produce = under
        most = max(most, hours)
        area += produce
        warm += UNDER_ICE_K * produce
        weight += produce
    if most < MIN_HOURS_YR:
        return NONE
    return Refuge(most, area, warm / weight if weight > 0.0 else None, plants, plants_cover)


def under_ice_light_share(values):
    """The share of the usable light that gets through the ice lid of an `ice-covered ocean` (`ice_light.transmission` of sea ice, GEN.1726, rows IL1 to IL7) from its stored `ice_shell_km`; 0 for any other
    state or a lid of unknown thickness."""
    if values.get("hydrosphere") != "ice-covered ocean" or values.get("ice_shell_km") is None:
        return 0.0
    return ice_light.transmission(1000.0 * max(float(values["ice_shell_km"]), 0.0), kind=ice_light.SEA)


def _under_ice(values, state, lit):
    """`(hours, production share)` of the sea under a thin lid, or `None`: the light that gets through (`under_ice_light_share`) must be over `LIGHT_FLOOR_W_M2`, and production is the ocean's
    share of the surface times the transmitted share, linear below the full light (LA4)."""
    if state != "ice-covered ocean":
        return None
    through = under_ice_light_share(values)
    light = surface_light_w_m2(values)
    sea = float(values.get("ocean_fraction") or values.get("ground_ocean_fraction") or 0.0)
    if through <= 0.0 or not light or through * light < LIGHT_FLOOR_W_M2 or sea <= 0.0:
        return None
    return YEAR_HOURS * lit, sea * through


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
