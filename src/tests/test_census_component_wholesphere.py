# tests/test_census_component_wholesphere.py

"""
TEST.185: the whole-sphere tests of the component side, on the census harness (`tests.audit`): nothing that needs gas or
young stars in a halo or an old cluster cell; every D or E region holds its members inside its radius at any distance of
the reach; the pre-placed share of every population; and no star of a cluster counted twice between the passes and the
fill. TEST.184 holds the nebula and Wolf-Rayet side; TEST.159 is done when both are.

The fills run with every per-star rate dialed up (`tests.audit.representation`), so a forbidden kind has a chance to
appear and the disk, where it is allowed, shows it does.
"""

import collections
import math

import pytest

from planetgen import tuning
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import component_pass, nebula_contents
from planetgen.generation import star_population as population
from planetgen.galaxy import scatter_limits
from planetgen.physics import hii, stellar_evolution as evolution
from planetgen.physics.units import ly_to_pc
from tests.audit import representation as rep
from tests.test_emission_members import D_REGIONS, E_REGIONS, SEED

GAS_AND_YOUNG_KINDS = ("nebula", "molecular-cloud", "supernova-remnant")
"""The kinds a place with no gas and no young stars never makes: emission and reflection nebulae grown round a hot
star, molecular clouds, and the supernova remnants (a core-collapse one needs young stars, a Type Ia one the gas to
be seen in: GEN.236)."""

OLD_REGIONS = ("halo", "globular-cluster", "intergalactic")
"""The regions of old stars with no cold gas: the halo, a cluster away from the disk, and the floor between galaxies."""


FILLED_SECTORS = 30
"""int: Sectors filled in each region: the rates are dialed up, so thirty show every kind the disk makes, and the fixture takes a third of the
nine minutes a hundred took on four cores (TEST.464)."""


@pytest.fixture(scope="module")
def audit_run():
    """The boosted census of every region, with every hot star hosting a nebula (the rules' chances at 1)."""
    patch = pytest.MonkeyPatch()
    try:
        patch.setattr(tuning, "NEBULA_HOST_RULES", tuple((*rule[:4], 1.0) for rule in tuning.NEBULA_HOST_RULES))
        yield rep.run(patch, sectors=FILLED_SECTORS)
    finally:
        patch.undo()


# --- Nothing hosted or made of gas in a halo or an old cluster cell -----------------------------------------------------


def test_the_regions_are_the_populations_they_are_named_for():
    assert rep.dominant_population("halo") == "halo"
    assert rep.dominant_population("globular-cluster") == "cluster"
    assert rep.dominant_population("satellite") == "satellite"
    assert rep.dominant_population("bulge") == "bulge"


@pytest.mark.parametrize("name", OLD_REGIONS)
def test_a_halo_or_old_cluster_cell_makes_no_emission_reflection_or_core_collapse_remnant(audit_run, name):
    census, _fill, _names = audit_run
    region = census[name]
    assert region.sectors == FILLED_SECTORS
    for kind in GAS_AND_YOUNG_KINDS:
        assert region.observed(kind) == 0, (name, kind)
    factors = rep.factors_at(name)
    for kind in GAS_AND_YOUNG_KINDS:
        assert factors[kind] < 0.01, (name, kind, factors[kind])


def test_the_same_fills_show_those_kinds_where_the_gas_and_the_young_stars_are(audit_run):
    """Without this the test above could pass for the wrong reason (a rate that never rolled)."""
    census, _fill, _names = audit_run
    for name in ("thin-disk-arm", "satellite"):
        for kind in ("molecular-cloud", "supernova-remnant"):
            assert census[name].observed(kind) > 0, (name, kind)


def test_a_core_collapse_remnant_needs_a_young_population_and_an_old_one_makes_only_type_ia():
    for name in OLD_REGIONS:
        assert rep.factors_at(name)["supernova-remnant"] < 0.01
    # The Type Ia share of a place with no young stars is one (GEN.236), so whatever remnant it did make has no core.
    assert rep.factors_at("thin-disk-arm")["supernova-remnant"] > 0.5


# --- Every D or E region holds its members inside its radius, at any distance of the reach ------------------------------


REACH_POINTS_PC = [(8.0e3, 0.0, 0.0), (5.0e4, 3.0e4, 1.0e3), (-1.5e5, 1.0e5, 20.0), (0.0, -2.6e5, 100.0), (2.4e5, 1.0e5, 0.0)]


def test_every_d_and_e_region_holds_its_members_inside_its_radius_wherever_it_sits():
    seen = collections.Counter()
    for letter, regions in (("D", D_REGIONS), ("E", E_REGIONS)):
        for number, nebula in enumerate(regions):
            point = REACH_POINTS_PC[number % len(REACH_POINTS_PC)]
            members = nebula_contents.emission_members(SEED, nebula, "O5V", point)
            for member in members:
                assert math.dist(member.position_pc, point) <= ly_to_pc(nebula.radius_ly) * (1 + 1e-9), (letter, point)
                assert evolution.main_sequence_lifetime_gy(member.mass_sol) > nebula.age_myr / 1000.0
            seen[letter] += len(members)
    assert seen["D"] > 50 and seen["E"] > 50


def test_a_region_holds_the_same_kinds_of_members_wherever_it_is_in_the_galaxy():
    """Only the draws move with the seed star's point (its seed is the point), never the counts by type."""
    for nebula in E_REGIONS[:3] + D_REGIONS[:3]:
        counts = {point: hii.host_counts([m.mass_sol for m in nebula_contents.emission_members(SEED, nebula, "O5V", point)])
                  for point in REACH_POINTS_PC}
        assert len({tuple(sorted(found.items())) for found in counts.values()}) == 1, nebula.hosts


# --- The pre-placed share of every population ---------------------------------------------------------------------------


FLOOR_SOL = tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL
"""The default bright-star floor, 18,000 L_sun."""
MASS_CUT_SOL = 14.0
"""The default mass pass limit: every star born at or above it is placed."""


def _independent_massive_share(low_age_gy, high_age_gy):
    """The share of the living stars of a population, ages uniform in the window, that were born at or above the mass
    pass limit: the IMF over the lifetimes, from the textbook mass-luminosity-lifetime tables and not from the code
    that computes the share."""
    breaks, slopes = tuning.IMF_BREAKS_SOL, tuning.IMF_SLOPES
    factors = [1.0]
    for i in range(1, len(slopes)):
        factors.append(factors[-1] * breaks[i] ** (slopes[i] - slopes[i - 1]))
    low, high = math.log(breaks[0]), math.log(breaks[-1])
    steps = 6000
    total = massive = 0.0
    for k in range(steps):
        mass = math.exp(low + (k + 0.5) * (high - low) / steps)
        i = max(j for j in range(len(slopes)) if breaks[j] <= mass)
        weight = factors[i] * mass ** (-slopes[i]) * mass          # per unit of ln(mass)
        lifetime = evolution.remnant_formation_time_gy(mass)
        alive = min(1.0, max(0.0, (lifetime - low_age_gy) / (high_age_gy - low_age_gy)))
        total += weight * alive
        if mass >= MASS_CUT_SOL:
            massive += weight * alive
    return massive / total


def test_every_population_has_a_share_between_zero_and_one_and_the_old_ones_have_none():
    shares = {name: population.placed_star_fraction(FLOOR_SOL, MASS_CUT_SOL, name)
              for name in population.POPULATIONS}
    assert all(0.0 <= share <= 1.0 for share in shares.values())
    for name in ("old", "bulge", "halo", "satellite", "nuclear_cluster", "nuclear_disk"):
        assert shares[name] == 0.0, name
    # 10 to 13 Gyr populations hold no star at or above the floor and none born above 14 solar masses, except a globular cluster: its own pre-placed stars
    # start at the faint edge of its horizontal branch (GEN.161, `tuning.CLUSTER_BRIGHT_FLOOR_SOL`), so its giants and branch stars are placed at plan time.
    for name, (low, _high) in tuning.STELLAR_POPULATION_AGE_RANGES_GY.items():
        if low >= 3.0 and name != "cluster":
            assert shares[name] == 0.0, name
    assert shares["cluster"] > 0.0


def _independent_cluster_bright_share(fe_h=-1.5):
    """The share of a globular cluster's living stars at or above the horizontal branch's faint edge, from counting the stars of a coeval population and not
    from the code that computes the share. Stars leave the main sequence at the turn-off at the rate the IMF at the turn-off mass times the speed the turn-off
    mass falls with age gives (mid-age of the cluster window, the main-sequence lifetime of the birth [Fe/H]); the stars above the floor are those that are
    in a phase brighter than it, whose time is the horizontal branch's (`HORIZONTAL_BRANCH_DURATION_GY`), the red giants brighter than it (1 in 1.44 +- 0.06 of
    the branch's stars, [R3879]) and the asymptotic giants (`AGB_MS_FRACTION` of the main-sequence lifetime). The living stars are the IMF integrated over
    the stars still alive at an age uniform in the window."""
    low_age, high_age = tuning.STELLAR_POPULATION_AGE_RANGES_GY["cluster"]
    breaks, slopes = tuning.IMF_BREAKS_SOL, tuning.IMF_SLOPES
    factors = [1.0]
    for i in range(1, len(slopes)):
        factors.append(factors[-1] * breaks[i] ** (slopes[i] - slopes[i - 1]))

    def per_mass(mass):
        i = max([j for j in range(len(slopes)) if breaks[j] <= mass * (1.0 + 1.0e-12)] or [0])
        return factors[i] * mass ** (-slopes[i])
    steps, low, high = 6000, math.log(breaks[0]), math.log(breaks[-1])
    living = 0.0
    for k in range(steps):
        mass = math.exp(low + (k + 0.5) * (high - low) / steps)
        alive = min(1.0, max(0.0, (evolution.remnant_formation_time_gy(mass, fe_h) - low_age) / (high_age - low_age)))
        living += per_mass(mass) * mass * alive * (high - low) / steps
    age = 0.5 * (low_age + high_age)
    lower, upper = 0.3, 1.5                             # the turn-off mass: the star whose main sequence ends at this age (the lifetime falls with mass)
    for _ in range(80):
        middle = 0.5 * (lower + upper)
        lower, upper = (middle, upper) if evolution.main_sequence_lifetime_gy(middle, fe_h) > age else (lower, middle)
    turn_off, h = 0.5 * (lower + upper), 1.0e-3
    slope = (evolution.main_sequence_lifetime_gy(turn_off * (1.0 + h), fe_h) - evolution.main_sequence_lifetime_gy(turn_off * (1.0 - h), fe_h)) / (2.0 * h * turn_off)
    leaving_per_gy = per_mass(turn_off) * abs(1.0 / slope)
    branch = tuning.HORIZONTAL_BRANCH_DURATION_GY
    giants_above_branch = branch / 1.44
    asymptotic = tuning.AGB_MS_FRACTION * age
    return leaving_per_gy * (branch + giants_above_branch + asymptotic) / living


def test_a_clusters_pre_placed_share_is_its_branch_and_the_giants_over_the_branchs_floor():
    """GEN.161: the cluster's share at its floor of 35 L_sun matches the independent count of the stars brighter than it, within the spread of its inputs (the
    red giants to branch stars ratio 1.44 +- 0.06, and the turn-off of one mid-age taken for the window)."""
    share = population.placed_star_fraction(FLOOR_SOL, MASS_CUT_SOL, "cluster")
    assert share == pytest.approx(_independent_cluster_bright_share(), rel=0.10)
    assert share == pytest.approx(population.bright_star_fraction(tuning.CLUSTER_BRIGHT_FLOOR_SOL, "cluster"), rel=1e-12)


def test_the_young_populations_share_follows_the_mass_function_over_the_lifetimes():
    young = population.massive_star_fraction(MASS_CUT_SOL, "young")
    low, high = tuning.STELLAR_POPULATION_AGE_RANGES_GY["young"]
    assert young == pytest.approx(_independent_massive_share(low, high), rel=0.4)
    assert 1.0e-4 < young < 1.0e-3


def test_the_share_falls_with_the_age_of_the_window():
    young = population.placed_star_fraction(FLOOR_SOL, MASS_CUT_SOL, "young")
    intermediate = population.placed_star_fraction(FLOOR_SOL, MASS_CUT_SOL, "intermediate")
    assert young > intermediate > 0.0
    assert population.massive_star_fraction(MASS_CUT_SOL, "intermediate") == 0.0   # nothing born above 14 outlives 100 Myr
    for earlier, later in zip(tuning.OPEN_CLUSTER_POPULATIONS[:5], tuning.OPEN_CLUSTER_POPULATIONS[1:6]):
        assert population.massive_star_fraction(MASS_CUT_SOL, earlier) >= population.massive_star_fraction(
            MASS_CUT_SOL, later) - 1.0e-12


def test_a_cells_pre_placed_share_is_the_mixs_weighted_share_of_its_populations(audit_run):
    for name in rep.REGION_POINTS:
        context = brightStars.FillContext(rep.profile.profile_at(rep.REGION_POINTS[name][0], rep.WITH_COMPONENTS,
                                                                 rep.EDGE, rep.E1),
                                          min_luminosity_sol=FLOOR_SOL, star_mass_limit_sol=MASS_CUT_SOL)
        total = sum(context.densities.values())
        wanted = sum(density * population.placed_star_fraction(FLOOR_SOL, MASS_CUT_SOL, pop,
                                                               brightStars.age_laws_of(context.shape))
                     for pop, density in context.densities.items()) / total
        assert context.bright_share() == pytest.approx(wanted, rel=1e-9, abs=1e-15)
        cluster_share = population.bright_star_fraction(tuning.CLUSTER_BRIGHT_FLOOR_SOL, "cluster")
        if name == "globular-cluster":
            # A cluster's own stars are placed from the horizontal branch's floor (GEN.161): the share is the mix's, none of whose populations keeps more than the cluster.
            assert 0.5 * cluster_share < context.bright_share() <= cluster_share * (1.0 + 1.0e-9), name
        elif name in OLD_REGIONS + ("thick-disk",):
            assert context.bright_share() < 1.0e-5, name


# --- No cluster star is counted twice between the passes and the fill ---------------------------------------------------


CLUSTER_FLOOR_SOL = 300.0
"""A floor the cluster's giants reach (a few percent of an old population's stars are this bright)."""


def test_the_stars_the_pass_places_and_the_stars_the_fill_draws_add_up_to_the_clusters_systems():
    cluster = rep.CLUSTER
    fractions = brightStars.band_fractions(CLUSTER_FLOOR_SOL)
    bands = [scatter_limits.PassBand(CLUSTER_FLOOR_SOL, None, None)]
    share = brightStars.placed_star_fraction(CLUSTER_FLOOR_SOL, None, "cluster")
    systems = component_pass.systems_of(cluster)
    wanted = systems * share
    assert wanted > 100
    placed = len(list(component_pass.scatter_stars([cluster], rep.EDGE, 11, fractions, bands)))
    assert placed == pytest.approx(wanted, abs=4.0 * math.sqrt(wanted))
    # The fill of a cell inside the cluster draws what the pass left: its share removed from the cell's own count.
    context = brightStars.FillContext(rep.profile.profile_at(cluster.center_pc, rep.WITH_COMPONENTS, rep.EDGE, rep.E1),
                                      min_luminosity_sol=CLUSTER_FLOOR_SOL)
    assert context.bright_share() == pytest.approx(share, rel=1e-9)
    drawn_by_the_fill = systems * (1.0 - context.bright_share())
    assert placed + drawn_by_the_fill == pytest.approx(systems, abs=4.0 * math.sqrt(wanted))


def test_the_pass_places_the_same_stars_on_every_call_and_never_in_a_filled_sector():
    cluster = rep.CLUSTER
    fractions = brightStars.band_fractions(CLUSTER_FLOOR_SOL)
    bands = [scatter_limits.PassBand(CLUSTER_FLOOR_SOL, None, None)]
    first = list(component_pass.scatter_stars([cluster], rep.EDGE, 11, fractions, bands))
    assert first == list(component_pass.scatter_stars([cluster], rep.EDGE, 11, fractions, bands))
    skipped = {tuple(row[:3]) for row in first[:40]}
    kept = list(component_pass.scatter_stars([cluster], rep.EDGE, 11, fractions, bands, skip_addresses=skipped))
    assert kept and not any(tuple(row[:3]) in skipped for row in kept)
