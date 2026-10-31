# planetgen/physics/field_rates.py

"""
How often an interstellar asteroid field is born (GEN.1101; docs/design/systems/interstellar-asteroid-fields.md sections 2, 4 and 8, docs/design/research/interstellar-field-sources.md)
==============================================================================================================================================

The physical rates and the host rule, nothing else. Boss decided (4:54 PM EDT 2026-10-30) that fields come only from physical sources: no game rate, no floor dial, no field drawn at a
made-up per-sector rate. So a sector holds a field when a source of the Galaxy puts one there, and in most sectors that is never.

* **S1, a collision of two rogue planets** (`collision_field_rate_per_sector_per_year`): `0.5 n V Gamma`, the pairs in a sector volume `V`, with the share of the collisions whose speed lies
  in the rubble window of `rogue_collision.rubble_speed_window_ms`. About 2e-16 standing fields per sector (a field for every 5e15 sectors); the module only does the arithmetic.
* **S2, a belt unbound by a supernova** (`unbound_belt_birth_rate_per_year`): the core-collapse rate of the Galaxy (`GALAXY_CORE_COLLAPSE_RATE_PER_YEAR`) times the share of the progenitors that
  hold a debris belt at their age (`belt_share`, Su et al. 2006 for A stars) times the rule that more than half of the mass must leave (`belt_is_unbound`, IF10). With the belt living the
  time `field_dispersal.retirement_age_yr` it needs to leave a light year (27 kyr at 11 km/s), about 200 fields stand in the whole Galaxy (`standing_number`), 3e-9 to 8e-9 per sector.
  **Hosted**: the field does not exist apart from its remnant. `hosted_belt` is the rule the supernova remnant draws from: the remnant of age `t` hosts an unbound belt of age `t` when
  the progenitor kept a belt (`belt_unit` against `belt_share`), the mass loss unbinds it and `t` is under the retirement age.
* **S3, a young cluster stripping the discs of its stars** (`flyby_stripping_rate_per_star_per_year`, `cluster_swarm_lifetime_yr`, `cluster_swarm_standing`): a flyby closer than the belt
  radius over `TRUNCATION_FRACTION` takes the belt away from the star, with the gravitational focusing of the pair; the swarm spreads at the orbital speed it had and so leaves a light
  year in about 70 thousand years. Standing swarms are the stripping rate times that life.
* **Not sources.** A rogue disrupted by a neutron star or black hole (S4) feeds the pulsar-planet disc, an ejected belt (S5) joins the interstellar background, a destroyed planet (S6) is
  bound and already a belt, a disrupted giant (S7) leaves gas. They have no function here.

Stand-ins are IF19 to IF22 of `docs/design/research/gap-filling-plan.md` section 6.106. Nothing reads the module yet: the supernova host rule and the columns (GEN.1103, GEN.1104) wire it.
"""

import collections
import math

from planetgen.physics import constants, field_dispersal, rogue_collision

GALAXY_CORE_COLLAPSE_RATE_PER_YEAR = 0.02
"""float: Core-collapse supernovae in the Galaxy per year, 2 per century (Diehl et al. 2006 from 26Al gives 1.9 +/- 1.1 per century, recalled). IF19, tier T2, L0."""
PROGENITOR_BELT_SHARE_BY_AGE = ((3.0e7, 0.41), (1.89e8, 0.25), (4.0e8, 0.15), (math.inf, 0.02))
"""tuple: `(upper age, share)` of A stars of 2 to 3 solar masses with a 24 micron excess by age, years (Su et al. 2006, Table 4, row IF9); the progenitors of a core-collapse supernova
live 10 to 40 Myr, so nearly all of them take the first row."""
TRUNCATION_FRACTION = 0.3
"""float: A flyby of periastron `q` truncates a disc at about `0.3 q` (Holman and Wiegert style cut; Breslau et al. 2014, recalled); a belt of radius `a` is stripped by `q < a / 0.3`. IF20,
tier T2, L1."""
CLUSTER_BELT_SHARE = 0.32
"""float: The share of the stars of a young cluster that hold a debris belt: the overall 32 +/- 5 percent of the Spitzer sample of A stars (IF9), used for solar-type stars too until a sample is
read. IF21, tier T3, L1."""
BELT_UNIT_UPPER = 1.0

SECONDS_PER_YEAR = constants.SECONDS_PER_YEAR
PC_PER_YEAR_AT_ONE_MS = rogue_collision.PC_PER_YEAR_AT_ONE_MS

HostedBelt = collections.namedtuple("HostedBelt", "speed_ms age_yr belt_radius_m")
"""namedtuple: A belt a remnant hosts: the speed at infinity it leaves at, its age (the remnant's) and the radius it started from."""


class FieldRatesError(ValueError):
    """A density, speed, mass, age or draw is out of range."""


def _positive(**values):
    for name, value in values.items():
        if not value > 0.0:
            raise FieldRatesError(f"{name} {value!r} is not positive")


# --- S1: collision of two rogues -------------------------------------------------------------------------------------------

def collision_field_rate_per_sector_per_year(rogue_density_pc3, sector_volume_pc3, mass1_kg, mass2_kg, radius1_m, radius2_m, speed_ms, dispersion_ms):
    """
    Fields born from rogue-rogue collisions in one sector per year: the collisions `0.5 n V Gamma` of the sector (`Gamma` the rate per rogue of `rogue_collision.rate_per_rogue_per_year`
    at the characteristic speed `speed_ms`) times the share whose speed lies in the rubble window of the pair (`rogue_collision.speed_fractions` at the dispersion `dispersion_ms`, with the
    focusing of the pair). Two Earths in the disc (0.78 per pc3, 40 km/s, dispersion 25 km/s) give 1e-23 collisions per rogue per year, about 2e-22 fields per sector per year for a cube of 4 pc.
    """
    _positive(rogue_density_pc3=rogue_density_pc3, sector_volume_pc3=sector_volume_pc3, speed_ms=speed_ms, dispersion_ms=dispersion_ms)
    low, high = rogue_collision.rubble_speed_window_ms(max(mass1_kg, mass2_kg), min(mass1_kg, mass2_kg))
    escape = rogue_collision.escape_speed_touching_ms(mass1_kg, mass2_kg, radius1_m, radius2_m)
    in_window = rogue_collision.speed_fractions(dispersion_ms, escape, (low, high))[1]
    per_rogue = rogue_collision.rate_per_rogue_per_year(rogue_density_pc3, mass1_kg, mass2_kg, radius1_m, radius2_m, speed_ms)
    return 0.5 * rogue_density_pc3 * sector_volume_pc3 * per_rogue * in_window


# --- S2: belt unbound by a supernova ---------------------------------------------------------------------------------------

def belt_share(progenitor_age_yr):
    """The share of the progenitors that hold a debris belt at the age `progenitor_age_yr` (41 percent up to 30 Myr, 25 percent to 189 Myr, 15 percent to 400 Myr, 2 percent after)."""
    if progenitor_age_yr < 0.0:
        raise FieldRatesError(f"the age {progenitor_age_yr!r} is negative")
    for upper, share in PROGENITOR_BELT_SHARE_BY_AGE:
        if progenitor_age_yr <= upper:
            return share
    return PROGENITOR_BELT_SHARE_BY_AGE[-1][1]


def belt_is_unbound(progenitor_mass_kg, remnant_mass_kg):
    """True when the sudden loss of the supernova frees a belt on a circular orbit: the remnant keeps under half of the mass (IF10)."""
    _positive(progenitor_mass_kg=progenitor_mass_kg, remnant_mass_kg=remnant_mass_kg)
    return remnant_mass_kg < progenitor_mass_kg / 2.0


def unbound_belt_birth_rate_per_year(progenitor_age_yr, core_collapse_rate_per_year=GALAXY_CORE_COLLAPSE_RATE_PER_YEAR, unbound_share=1.0):
    """
    Unbound belts born per year in a region that holds the core-collapse rate `core_collapse_rate_per_year` (the Galaxy's by default): that rate times `belt_share` at the progenitor's age
    times `unbound_share`, the share of the remnants that free the belt (1 for neutron stars, which every progenitor of 8 to 20 solar masses leaves).
    """
    _positive(core_collapse_rate_per_year=core_collapse_rate_per_year)
    if not 0.0 <= unbound_share <= 1.0:
        raise FieldRatesError(f"the unbound share {unbound_share!r} is not in 0 to 1")
    return core_collapse_rate_per_year * belt_share(progenitor_age_yr) * unbound_share


def standing_number(birth_rate_per_year, lifetime_yr):
    """The fields standing at any time: the birth rate times the life of a field (a Poisson mean). 0.02 x 0.41 per year and 27 kyr give about 220 in the Galaxy."""
    if birth_rate_per_year < 0.0 or lifetime_yr < 0.0:
        raise FieldRatesError("the rate and the lifetime are not negative")
    return birth_rate_per_year * lifetime_yr


def standing_per_sector(standing_in_galaxy, stars_in_sector, stars_in_galaxy):
    """The standing fields to a sector of `stars_in_sector` stars in a Galaxy of `stars_in_galaxy`: the same share of the stars as of the fields."""
    _positive(stars_in_galaxy=stars_in_galaxy)
    if stars_in_sector < 0.0 or standing_in_galaxy < 0.0:
        raise FieldRatesError("the star count and the standing number are not negative")
    return standing_in_galaxy * stars_in_sector / stars_in_galaxy


def hosted_age_share(remnant_age_range_yr, belt_speed_ms, initial_radius_m=0.0):
    """
    The share of the remnants of an age range `(young, old)`, uniform in age, that are younger than a field lives at the speed `belt_speed_ms` and so still host it: 27 kyr of the
    50 to 100000 years of `tuning.SUPERNOVA_REMNANT_AGE_RANGE_YEARS` at 11 km/s, about a quarter.
    """
    young, old = remnant_age_range_yr
    if not 0.0 <= young < old:
        raise FieldRatesError("the remnant age range must be ascending and not negative")
    limit = field_dispersal.retirement_age_yr(belt_speed_ms, initial_radius_m)
    return min(max(limit - young, 0.0), old - young) / (old - young)


def hosted_belt(belt_unit, progenitor_age_yr, progenitor_mass_kg, remnant_mass_kg, belt_radius_m, remnant_age_yr):
    """
    The belt a supernova remnant hosts, or `None`. `belt_unit` is a uniform draw in `[0, 1)` of the remnant: the progenitor kept a belt when it is under `belt_share(progenitor_age_yr)`. The
    belt is unbound when the remnant keeps under half of the mass (`belt_is_unbound`), it leaves at `field_dispersal.unbound_belt_speed_ms`, and the remnant hosts it while its own age is
    under the time that belt needs to pass the retirement radius (`field_dispersal.is_retired`, from the radius of the belt). Past that the field is retired and the remnant hosts nothing.
    """
    if not 0.0 <= belt_unit < BELT_UNIT_UPPER:
        raise FieldRatesError(f"the draw {belt_unit!r} is not in [0, 1)")
    _positive(belt_radius_m=belt_radius_m)
    if remnant_age_yr < 0.0:
        raise FieldRatesError(f"the age {remnant_age_yr!r} is negative")
    if belt_unit >= belt_share(progenitor_age_yr) or not belt_is_unbound(progenitor_mass_kg, remnant_mass_kg):
        return None
    speed = field_dispersal.unbound_belt_speed_ms(progenitor_mass_kg, remnant_mass_kg, belt_radius_m)
    if field_dispersal.is_retired(remnant_age_yr, speed, belt_radius_m):
        return None
    return HostedBelt(speed, remnant_age_yr, belt_radius_m)


# --- S3: young cluster stripping -------------------------------------------------------------------------------------------

def flyby_stripping_rate_per_star_per_year(star_density_pc3, dispersion_ms, belt_radius_m, star_mass_kg):
    """
    Flybys per star per year that come inside `belt_radius_m / TRUNCATION_FRACTION` and strip its belt: `n pi q^2 (1 + 2 G M / (q v^2)) v` for the relative speed `v` (the dispersion
    `dispersion_ms` times sqrt(2)) and the pair mass `M` of two stars of `star_mass_kg`, `q` the periastron. A belt of 50 AU in a cluster of 1000 stars per pc3 at 2 km/s: about 20 per Myr
    per 1000 stars.
    """
    _positive(star_density_pc3=star_density_pc3, dispersion_ms=dispersion_ms, belt_radius_m=belt_radius_m, star_mass_kg=star_mass_kg)
    speed = math.sqrt(2.0) * dispersion_ms
    periastron_m = belt_radius_m / TRUNCATION_FRACTION
    focusing = 1.0 + 2.0 * constants.G * 2.0 * star_mass_kg / (periastron_m * speed ** 2)
    area_pc2 = math.pi * (periastron_m / constants.PARSEC_M) ** 2 * focusing
    return star_density_pc3 * area_pc2 * speed * PC_PER_YEAR_AT_ONE_MS


def cluster_swarm_lifetime_yr(belt_radius_m, star_mass_kg):
    """
    How long a stripped swarm is a field, years: the time it needs to pass the retirement radius. The bodies keep the orbital speed they had round their star, `sqrt(G M / a)`, so the swarm
    spreads at that speed from the radius `belt_radius_m` (4.2 km/s for 50 AU round the Sun-mass star: 70 kyr), whatever the dispersion of the cluster.
    """
    _positive(belt_radius_m=belt_radius_m, star_mass_kg=star_mass_kg)
    speed = math.sqrt(constants.G * star_mass_kg / belt_radius_m)
    return field_dispersal.retirement_age_yr(speed, belt_radius_m)


def cluster_swarm_standing(star_count, star_density_pc3, dispersion_ms, belt_radius_m, star_mass_kg, belt_share_of_stars=CLUSTER_BELT_SHARE):
    """The stripped swarms standing in a cluster of `star_count` stars: stars with a belt times the stripping rate times the life of a swarm (0.5 for 1000 stars of a 1000-per-pc3 cluster)."""
    if star_count < 0.0 or not 0.0 <= belt_share_of_stars <= 1.0:
        raise FieldRatesError("the star count is not negative and the belt share is in 0 to 1")
    rate = flyby_stripping_rate_per_star_per_year(star_density_pc3, dispersion_ms, belt_radius_m, star_mass_kg)
    return standing_number(star_count * belt_share_of_stars * rate, cluster_swarm_lifetime_yr(belt_radius_m, star_mass_kg))
