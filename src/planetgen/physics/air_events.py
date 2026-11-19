# planetgen/physics/air_events.py

"""
The shocks of a world's history: flood basalts and large impacts on a stream of the body's own, and the share of life each kills (GEN.1746, slice 1; design:
reports/air-history-design-gen1739-2026-11-15.md, sections 3.4, 5.5 and 5.6)
=========================================================================================================================================

The air history of GEN.1745 follows a world's air through smooth laws of age. The record of the Earth is also a record of shocks: the large igneous provinces that coincide with four of the
largest extinctions of the last 260 Myr, the impacts, and the swings in the air they leave. This module is the first part of that layer: **what happens when**, as a Poisson process for each kind of
event whose rate follows the world's own state at that age, and **how much of life each event kills**. It reads the clocks of `air_history` and the heat of `interior_heat`; nothing in the generator reads
it yet (the later slices of GEN.1746 feed the air and the life ages from it), so the output of a run is the one it was.

* *Flood basalts.* One about every 25 Myr on the Earth (a dozen in the last 300 Myr), each about 2.5 million km3 of basalt erupted in about 1 Myr ([R5319], read for the numbers). The volcanism goes with the heat the interior gives
  off, so the rate is the Earth's times the body's heat flux at that age over the model Earth's (`surface_pressure.HEAT_FLUX_EARTH_W_M2`, the one Earth reference of the heat laws): a young body has more, a cold or tidally
  quiet one fewer, a body with no heat none. A lid that does not move still erupts from its plumes, so the lid is not asked.
* *Large impacts.* The cumulative rate of the Earth for bodies over a diameter is a power law through the two figures a secondary source gives: 1 km about every 0.5 Myr, 5 km about every 15 Myr. The rate of a body
  goes with its cross section (its radius squared over the Earth's, gravitational focusing left out) and with the flux over its history relative to now: the lunar cratering chronology of Neukum et al., as printed in
  [R5322], `N(1) = 5.44e-14 (exp(6.93 T) - 1) + 8.38e-4 T` for T in Ga before now, whose rate is `5.44e-14 x 6.93 exp(6.93 T) + 8.38e-4`. Only the impacts over the size the caller names are listed (the small ones
  are too many and kill nothing the record keeps).
* *The kill of an event.* The kill curve of the marine record gives the share of species killed by the events that recur at a given interval ([R5320]: about 10 percent once per Myr, 30 percent once per 10 Myr, 65
  percent once per 100 Myr; "does not imply any cause"). An event kills the share the curve gives for how often events as large as it recur on the Earth: an impact by the recurrence of its size, a flood basalt by its
  25 Myr. A body with more events than the Earth's has more kills, not larger ones.
* *Recovery.* The origination that follows an extinction lags it by about 10 Myr, whatever its size ([R5321]).

The events of a body are drawn on a stream keyed to the body (`flavor_tree.seed_of`, the one the flavors use: what is final about the body and its star, nothing about the order the system made its bodies in), so the
run's own stream does not move and the same body has the same events whenever it is made again. They are a draw of the history, not of the run.

Every number that is not read from a source is a design constant with its tier and tweak level in `docs/design/research/gap-filling-plan.md` (GEN.725 rule, section 6.313).
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import air_history, constants, flavor_tree, interior_heat, surface_pressure
from planetgen.util import draw

FLOOD_BASALT = "flood basalt"
IMPACT = "impact"
KINDS = (FLOOD_BASALT, IMPACT)
"""tuple: The kinds of event."""

FLOOD_BASALT_PER_GYR = 40.0
"""float: The flood basalts of the Earth a gigayear, 40 (one about every 25 Myr; about a dozen in the last 300 Myr, [R5319]). [Stand-in 6.313 EV1, tier T1 verified, tweak L0.]"""
FLOOD_BASALT_VOLUME_KM3 = 2.5e6
"""float: The modal volume of one flood basalt, km3 ([R5319]), erupted in `FLOOD_BASALT_YEARS` (its pulse of CO2 and SO2 is the next slice). [Stand-in 6.313 EV1, tier T1 verified, tweak L0.]"""
FLOOD_BASALT_YEARS = 1.0e6
"""float: The time the bulk of a flood basalt takes, years ([R5319]). [Stand-in 6.313 EV1, tier T1 verified, tweak L0.]"""
FLOOD_BASALT_RECURRENCE_MYR = 1000.0 / FLOOD_BASALT_PER_GYR
"""float: How often the Earth has a flood basalt, Myr: the interval at which the kill curve is read for one. [Stand-in 6.313 EV1.]"""

IMPACT_ONE_KM_PER_GYR = 2000.0
"""float: The Earth's impacts a gigayear of bodies over 1 km across: one about every 0.5 Myr. [Stand-in 6.313 EV3, tier T3 marginal unverified (two figures on one secondary page: the 5 km one is "10 to 20 Myr"), tweak L0.]"""
IMPACT_FIVE_KM_PER_GYR = 1000.0 / 15.0
"""float: The Earth's impacts a gigayear of bodies over 5 km across: one about every 15 Myr (the middle of 10 to 20). [Stand-in 6.313 EV3, tier T3 marginal unverified, tweak L0.]"""
IMPACT_SIZE_EXPONENT = math.log(IMPACT_ONE_KM_PER_GYR / IMPACT_FIVE_KM_PER_GYR) / math.log(5.0)
"""float: The slope of the cumulative size law, 2.1: the number of impacts over `D` goes as `D^-2.1` through the two figures above. [Stand-in 6.313 EV4, tier T7 educated guess based on inferred physics, tweak L0.]"""
IMPACT_LARGEST_KM = 100.0
"""float: The largest body the size law draws, km: over it an event is a giant impact, which is the formation's and not the history's. [Stand-in 6.313 EV8, tier T6 marginal educated guess based on physics, tweak L3.]"""

NEUKUM_A = 5.44e-14
NEUKUM_LAMBDA_PER_GYR = 6.93
NEUKUM_B_PER_GYR = 8.38e-4
"""float: The three constants of the lunar cratering chronology function `N(1) = a (exp(lambda T) - 1) + b T` for craters over 1 km a square kilometre at an age T Ga before now ([R5322], printed from Neukum et al.). [Stand-in 6.313 EV5, tier T2 confident unverified, tweak L0.]"""

KILL_CURVE = ((1.0, 0.10), (10.0, 0.30), (100.0, 0.65))
"""tuple: The kill curve of the marine record: `(recurrence interval Myr, share of species killed)` ([R5320], through a secondary article; the table itself is not read). [Stand-in 6.313 EV6, tier T2 confident unverified, tweak L0.]"""

RECOVERY_GYR = 0.010
"""float: The lag of the origination behind an extinction, Gyr, the same after the large ones as after the small ([R5321]). [Stand-in 6.313 EV7, tier T1 verified, tweak L0.]"""

EVENT_STEP_GYR = 0.05
"""float: The longest span over which an event rate is held constant when the events are drawn, Gyr (a numerical setting, defined by the integration: GEN.1739)."""


@dataclass(frozen=True)
class Event:
    """
    One shock of a world's history.

    Attributes:
        age_gyr (float): The age of the body when it happened, Gyr.
        kind (str): `FLOOD_BASALT` or `IMPACT`.
        size (float): The volume of the basalt, km3, or the diameter of the impactor, km.
        kill_share (float): The share of the world's species it kills (`kill_share`).
    """

    age_gyr: float
    kind: str
    size: float
    kill_share: float


def impact_flux_ratio(age_gyr):
    """The impact rate of a body `age_gyr` Gyr after it formed over the Earth's now: the derivative of the lunar cratering chronology function (`NEUKUM_A`, ...) at the time before now of that age on
    the Earth's own clock (`tuning.STELLAR_SOLAR_SYSTEM_AGE_GY` less `age_gyr`), so the flux falls as the leftover bodies of the disk are used up and is the Earth's present one, 1, from the Earth's age on."""
    before_now = max(tuning.STELLAR_SOLAR_SYSTEM_AGE_GY - max(age_gyr, 0.0), 0.0)
    rate = NEUKUM_A * NEUKUM_LAMBDA_PER_GYR * math.exp(NEUKUM_LAMBDA_PER_GYR * before_now) + NEUKUM_B_PER_GYR
    return rate / (NEUKUM_A * NEUKUM_LAMBDA_PER_GYR + NEUKUM_B_PER_GYR)


def flood_basalt_rate_per_gyr(heat_flux_w_m2):
    """The flood basalts a gigayear of a body whose interior gives off `heat_flux_w_m2` at that age: the Earth's rate times the flux over the model Earth's; 0 for no heat."""
    if not heat_flux_w_m2 or heat_flux_w_m2 <= 0.0:
        return 0.0
    return FLOOD_BASALT_PER_GYR * heat_flux_w_m2 / surface_pressure.HEAT_FLUX_EARTH_W_M2


def impact_rate_per_gyr(diameter_km, radius_m, flux_ratio=1.0):
    """The impacts a gigayear of bodies over `diameter_km` on a body of radius `radius_m` when the flux is `flux_ratio` times the present: the Earth's power law, times the cross section over the Earth's."""
    cross_section = (radius_m / constants.EARTH_RADIUS_M) ** 2
    return IMPACT_FIVE_KM_PER_GYR * (diameter_km / 5.0) ** (-IMPACT_SIZE_EXPONENT) * cross_section * flux_ratio


def kill_share(recurrence_myr):
    """The share of species killed by the events that recur every `recurrence_myr` Myr: the kill curve, straight in the log of the interval between its points and held at the ends (0.10 under 1 Myr, 0.65 over 100)."""
    first, last = KILL_CURVE[0], KILL_CURVE[-1]
    if recurrence_myr <= first[0]:
        return first[1]
    if recurrence_myr >= last[0]:
        return last[1]
    for (low_myr, low_share), (high_myr, high_share) in zip(KILL_CURVE, KILL_CURVE[1:]):
        if recurrence_myr <= high_myr:
            weight = math.log(recurrence_myr / low_myr) / math.log(high_myr / low_myr)
            return low_share + weight * (high_share - low_share)
    return last[1]


def impact_kill_share(diameter_km):
    """The share of species an impact of `diameter_km` kills: the kill curve at how often the Earth has impacts as large (`impact_rate_per_gyr` at the Earth's present flux)."""
    per_myr = impact_rate_per_gyr(diameter_km, constants.EARTH_RADIUS_M) / 1000.0
    return kill_share(1.0 / per_myr)


def flood_basalt_kill_share():
    """The share of species a flood basalt kills: the kill curve at the interval of the Earth's flood basalts."""
    return kill_share(FLOOD_BASALT_RECURRENCE_MYR)


def _waits(stream, rates_per_gyr, spans_gyr):
    """The ages, Gyr from the start of the first span, of a Poisson process whose rate is `rates_per_gyr[i]` over the `i`th of `spans_gyr`: an exponential wait of mean 1 is spent against the rate of each span in turn."""
    ages, start = [], 0.0
    budget = -math.log(1.0 - stream.random())
    for rate, span in zip(rates_per_gyr, spans_gyr):
        used = 0.0
        while rate > 0.0:
            wait = budget / rate
            if used + wait > span:
                budget -= rate * (span - used)
                break
            used += wait
            ages.append(start + used)
            budget = -math.log(1.0 - stream.random())
        start += span
    return ages


def spans(start_gyr, end_gyr):
    """The `(age, span)` pairs, Gyr, that cover `start_gyr` to `end_gyr` with no span over `EVENT_STEP_GYR` and none crossing a step of `air_history.schedule`: the age of each is the middle of its span."""
    if end_gyr <= start_gyr:
        return []
    edges = [start_gyr] + [age for age in air_history.schedule(end_gyr) if start_gyr < age < end_gyr] + [end_gyr]
    pairs = []
    for low, high in zip(edges, edges[1:]):
        count = max(1, math.ceil((high - low) / EVENT_STEP_GYR))
        pairs.extend((low + (index + 0.5) * (high - low) / count, (high - low) / count) for index in range(count))
    return pairs


def events_of(body, start_gyr=0.0, end_gyr=None, min_impact_km=5.0):
    """
    The events of `body`'s history from `start_gyr` to `end_gyr` (the star's age by default) in the order they happen, drawn on the body's own stream (`flavor_tree.seed_of`), none from the run's.

    The flood basalts follow the heat the body's interior gave off at each age (`interior_heat.heat_at`); the impacts of bodies over `min_impact_km` follow its cross section and the lunar chronology of
    the flux (`impact_flux_ratio`), each of a size drawn from the power law between `min_impact_km` and `IMPACT_LARGEST_KM`. A body with no star age, or no stored interior, has no flood basalts.
    """
    star_age = getattr(getattr(body, "star", None), "age", None)
    end_gyr = end_gyr if end_gyr is not None else star_age
    if not end_gyr or end_gyr <= start_gyr:
        return []
    radius_m = (getattr(body, "radius", None) or 0.0) * 1000.0      # the body's radius is in km
    pairs = spans(start_gyr, end_gyr)
    ages = [age for age, _ in pairs]
    lengths = [span for _, span in pairs]
    flood_rates, impact_rates = [], []
    for age in ages:
        heat = interior_heat.heat_at(body, age)
        flood_rates.append(flood_basalt_rate_per_gyr(heat.flux_w_m2 if heat is not None else 0.0))
        impact_rates.append(impact_rate_per_gyr(min_impact_km, radius_m, impact_flux_ratio(age)) if radius_m > 0.0 else 0.0)
    stream = draw.Stream(flavor_tree.seed_of("air-events", body))
    events = [Event(start_gyr + offset, FLOOD_BASALT, FLOOD_BASALT_VOLUME_KM3, flood_basalt_kill_share()) for offset in _waits(stream, flood_rates, lengths)]
    exponent = IMPACT_SIZE_EXPONENT
    for offset in _waits(stream, impact_rates, lengths):
        diameter = min(min_impact_km * (1.0 - stream.random()) ** (-1.0 / exponent), IMPACT_LARGEST_KM)
        events.append(Event(start_gyr + offset, IMPACT, diameter, impact_kill_share(diameter)))
    return sorted(events, key=lambda event: (event.age_gyr, event.kind))
