# planetgen/generation/run_sector.py

"""
Run Sector
==========

The `sector` command: one or more independent sectors, each a shared
`SpaceSector` of systems built through `run_system`, with a sparse,
science-based population of exotic phenomena, star-hosted nebulae, fast
stars and (at the galactic center) the nucleus. `generate_sector` is
also what the galaxy command runs per sector.
"""

import copy
import hashlib
import math
import os
import sys
import time
from collections import Counter

from planetgen.db import store
from planetgen.galaxy import scales
from planetgen.galaxy import seed as galaxySeed
from planetgen.generation import directives
from planetgen.generation import cloud as cloudContext
from planetgen.generation import fill_stages
from planetgen.generation import multistar, nebula_contents, scatter_companions, stored_trees
from planetgen.generation import system_pass
from planetgen.generation import bright_stars as brightStars
from planetgen.galaxy import remnant_density, ejecta_field, nebula_field
from planetgen.galaxy.frame import galaxy_to_local_pc, local_to_galaxy_pc
from planetgen.galaxy.system_position import random_unit_vector
from planetgen.physics import constants as physical_constants
from planetgen.physics import stellar_evolution
from planetgen import tuning as program_constants
from planetgen.util import draw
from planetgen.util import log
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.galaxy.sector import SpaceSector, _sample_poisson_count
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.system import StarSystem
from planetgen.names.wordsalad import generate_sector_name
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.util.random import log_uniform
from planetgen.generation import phenomenon_scatter
from planetgen.generation.phenomena import registry
from planetgen.generation.star_labels import describe_types, star_label
from planetgen.generation import followed, run_common
from planetgen.generation import run_phenomenon
from planetgen.generation import run_system
from planetgen.queue import commands, redisqueue, runrecord, units, wrapper


SECTOR_PHENOMENON_KINDS = tuple(
    (row.kind, row.phenomenon_type, registry.resolve(row.builder)) for row in registry.rows(sector_roll=True)
)
"""tuple: `(program_constants.PHENOMENON_DENSITY_PC3 key, phenomenon type,
factory)` for each kind `generate_sector_phenomena` rolls, read from the
phenomena registry (`registry.Kind.sector_roll`; GEN.250). A factory takes
`(config, galactic_center_dist_ly)`; only black holes and neutron stars
use the distance (their Hill sphere and galactic orbit). A planetary
nebula is generated around its own new hot white dwarf system
(`_add_planetary_nebula`). Runaway and hypervelocity stars are flags on
ordinary systems (`flag_fast_stars`), and emission and reflection nebulae
grow around a sector's own hot stars (`add_star_hosted_nebulae`). In a
galaxy-placed sector the "molecular-cloud" roll gives way to the galaxy's
cloud field (`nebulaField`, GEN.47). Diffuse gas (classes A-B) fills about
half the disk, so it is the background, not generated."""

POPULATION_AGED_KINDS = ("brown-dwarf", "neutron-star", "black-hole")
"""tuple: The kinds whose age follows the population the sector's fill gives them (GEN.553): they are made with
`FillContext.apply_population` and take their age from that population's window and age law."""


def iter_sector_configs(args):
    """
    Yields one `SystemConfig` per system in the sector, one at a time, sharing
    the same value options across all of them (via `build_system_config`, so
    this can never silently drift from what those options mean for a single
    system; sector runs take no forcing options, GEN.51). A huge count
    (`--num-systems 1000000000`, or a huge `--density` draw) never builds them
    all up front -- `generate_sector` stops pulling once the sector is full.
    If `--min-habitable` was given, `HABITABLE_WORLD = True` is set on that
    many randomly-chosen configs; the indices are chosen (and any conflict
    reported) before the first config is yielded.

    Args:
        args (argparse.Namespace): Parsed arguments, with `args.num_systems`
            already resolved to a concrete count -- either the user's explicit
            `--num-systems`, or (see `generate_sector`) a per-sector value
            sampled from `--density`.

    Raises:
        SystemExit: Under `--strict`, if `--min-habitable` exceeds a
            `--density`-driven count (an explicit `--num-systems` is checked
            earlier by `validate_shared_generation_args`); without it, a
            warning and the sector gets `--min-habitable` systems (GEN.81).
    """
    count = args.num_systems
    if args.min_habitable > count:
        # GEN.81: the sector gets the habitable systems asked for, so it
        # holds that many systems at least.
        run_common._refuse_or_warn(args, f"--min-habitable ({args.min_habitable}) exceeds this sector's drawn system "
                              f"count ({count}) (with --density the count is drawn per sector); it gets "
                              f"{args.min_habitable} systems instead.")
        count = args.min_habitable
    if count <= 0:
        return

    first = run_system.build_system_config(args)
    forced = set()
    if args.min_habitable > 0:
        forced = set(draw.sample(range(count), k=args.min_habitable))

    for i in range(count):
        config = first if i == 0 else run_system.build_system_config(args)
        if i in forced:
            config.HABITABLE_WORLD = True
        yield config


def sector_star_count(sector):
    """How many stars `sector`'s systems hold (a binary counts two) --
    what every per-star phenomenon rate multiplies."""
    return sum(len(getattr(entry.star_system, "stars", None) or [None]) for entry in sector.entries)


def generate_sector_phenomena(sector, args, galactic_center_dist_ly=None, cloud_field=None, scattered=False,
                              kind_factors=None, type_ia_share=None, ambient_density_cm3=None, neighbour_hosts=(),
                              fill=None):
    """
    Populates an already-built `sector` with its exotic phenomena, each
    kind in `SECTOR_PHENOMENON_KINDS` sampled independently by a Poisson
    draw (`spaceSector._sample_poisson_count`, the same mechanism
    `--density`'s own system count uses) whose mean is the kind's rate
    per star (`program_constants.phenomenon_rate_per_star`, from Boss's
    research density at the solar neighborhood) times the sector's star
    count. A star count already tracks the local stellar density, so a
    denser sector gets proportionally more. A local-density sector gets
    about 58 rogue planets and 2 brown dwarfs; rarer kinds mostly none.

    Black holes and neutron stars are real, stellar-mass gravitating
    bodies, so they're added via `SpaceSector.add_phenomenon`'s
    Hill-sphere-aware placement (never within a neighboring star system's
    or another compact remnant's own Hill sphere); every other phenomenon
    type has no comparable gravitational footprint at this generator's
    scale and is simply placed at a random point in the sector's cube --
    see that method's own docstring.

    Args:
        sector (SpaceSector): An already-populated sector (every star
            system already added via `add_system`, so their Hill spheres
            are in place for a black hole's/neutron star's own placement
            check) to add phenomena to, in place.
        args (argparse.Namespace): Parsed arguments (none is read here).
        galactic_center_dist_ly (float, optional): As in `generate_sector`
            -- threaded into a standalone black hole's/neutron star's own
            Hill-sphere/galactic-orbit calculation, exactly like every star
            in the sector already gets.
        cloud_field (tuple, optional): For a galaxy-placed sector, the
            `nebula_field.clouds_reaching` arguments `(galaxy_seed, shape,
            center_pc, reach_pc)`: its molecular clouds then come from the
            galaxy's cloud field (GEN.47, into `sector.field_nebulae`)
            instead of this sector's own `"molecular-cloud"` roll.
        scattered (bool): The galaxy's phenomenon scatter ran (GEN.100):
            the kinds in `phenomenon_scatter.SCATTERED_KINDS` were built
            from their scatter rows (`add_scattered_phenomena`), and those
            below its mass cut are drawn by `add_below_cut_remnants`
            (GEN.168), so none is rolled here.
        type_ia_share (float, optional): The sector profile's `type_ia_share` (GEN.319): the share of this sector's
            supernova remnants that are Type Ia. `None`: the solar circle's.
        ambient_density_cm3 (float, optional): The interior density of the superbubble the sector's profile places it
            in (`GasProfile.superbubble`, GEN.298): the supernova remnants it rolls are sized by it. `None`: outside a
            bubble, each remnant's own class density.
        kind_factors (dict, optional): The sector profile's `kind_factors`
            (GEN.236): each kind's mean is its rate per star times the star
            count times its factor, and the star-hosted nebulae's chance is
            scaled by the `"nebula"` factor. A kind without one, and a
            sector with none given (a standalone sector), keeps factor 1.

        neighbour_hosts (iterable): The emission regions (D and E) lit by pre-placed hot stars in the neighbouring sectors,
            `(nebula, host_name, host_position_pc)` each (`FillContext.emission_hosts`, GEN.321): this sector builds the
            ionizing members of each whose points fall inside it.
        fill (brightStars.FillContext, optional): For a galaxy-placed sector, its fill: a brown dwarf takes a
            population from it, and so its age follows that population's star formation (GEN.553).

    Returns:
        list: The newly created `SectorPhenomenonEntry` instances. May hold
             fewer than the Poisson draw sampled for a massive type (a
             black hole/neutron star) if the sector was too full/small to
             fit it without violating another massive object's Hill sphere
             -- that one draw is silently skipped (see `ValueError` below)
             rather than crashing the whole sector's generation over a
             single rare phenomenon that didn't fit.
    """
    star_count = sector_star_count(sector)
    new_entries = []
    kind_factors = kind_factors or {}

    for kind, phenomenon_type, factory in SECTOR_PHENOMENON_KINDS:
        if kind == "molecular-cloud" and cloud_field is not None:
            if sector.cloud_field != cloud_field:
                cloudContext.begin(sector, cloud_field, nebula_field.clouds_reaching(*cloud_field))
            fill_stages.count(sector, fill_stages.NEBULAE, kind, len(sector.field_nebulae))
            continue
        if scattered and kind in phenomenon_scatter.SCATTERED_KINDS:
            continue
        # Each kind draws from a stream of its own (stage E), so adding a kind moves no other kind's draws.
        stage = fill_stages.NEBULAE if kind == "molecular-cloud" else fill_stages.PER_STAR_KINDS
        with fill_stages.stage(sector, stage, kind):
            rate = (remnant_density.rate_per_star(kind_factors, kind) if kind in remnant_density.KINDS
                    else program_constants.phenomenon_rate_per_star(kind) * kind_factors.get(kind, 1.0))
            count = _sample_poisson_count(rate * star_count)
            for _ in range(count):
                phenomenon_config = SystemConfig()
                phenomenon_config.TYPE_IA_SHARE = type_ia_share
                phenomenon_config.AMBIENT_DENSITY_CM3 = ambient_density_cm3
                if kind in POPULATION_AGED_KINDS and fill is not None:
                    fill.apply_population(phenomenon_config)   # GEN.553: its age follows its population
                phenomenon = factory(phenomenon_config, galactic_center_dist_ly)
                if kind == "planetary-nebula":
                    entry = _add_planetary_nebula(sector, args, phenomenon, galactic_center_dist_ly)
                    if entry is not None:
                        new_entries.append(entry)
                        fill_stages.count(sector, stage, kind)
                        fill_stages.count(sector, stage, "star-system")
                    continue
                try:
                    new_entries.append(sector.add_phenomenon(phenomenon, phenomenon_type))
                    fill_stages.count(sector, stage, kind)
                except ValueError:
                    # Only a massive type (black-hole/neutron-star) can raise
                    # here (see SpaceSector._random_position) -- no room left
                    # to place it without violating another massive object's
                    # Hill sphere. Skip just this one draw rather than this
                    # phenomenon type, or the whole sector.
                    continue

    own_hosts = []
    with fill_stages.stage(sector, fill_stages.NEBULAE, "hosted"):
        hosted = add_star_hosted_nebulae(sector, args, chance_factor=kind_factors.get("nebula", 1.0),
                                         hosts_out=own_hosts)
    fill_stages.count(sector, fill_stages.NEBULAE, "nebula", len(hosted))
    new_entries.extend(hosted)
    if cloud_field is not None:
        with fill_stages.stage(sector, fill_stages.NEBULAE, "protostars"):
            add_cloud_protostars(sector, args, cloud_field, galactic_center_dist_ly)
        with fill_stages.stage(sector, fill_stages.NEBULAE, "emission-members"):
            add_emission_members(sector, args, cloud_field, own_hosts, neighbour_hosts, galactic_center_dist_ly)
    return new_entries


def _spectral_code(star):
    """A star's spectral code (`"G2V"`), the first word of `Star.type`."""
    return (getattr(star, "type", "") or "").split(" ", 1)[0]


def _add_planetary_nebula(sector, args, nebula, galactic_center_dist_ly=None, position=None):
    """
    Adds `nebula` (a planetary nebula) around a new system whose star is
    its hot central white dwarf (`PLANETARY_NEBULA_CENTRAL_STAR_TYPES`),
    placed like any other system, or at `position` (a scattered
    nebula's stored point, kept where it was put). Returns the nebula's
    entry, or `None` when the sector has no room left for the system.
    """
    config = SystemConfig()
    config.STAR_TYPE = draw.choice(program_constants.PLANETARY_NEBULA_CENTRAL_STAR_TYPES)
    nebula_contents.planetary_central_config(config)  # GEN.320: a close companion in 20%
    system = StarSystem(system_config=config, galactic_center_dist_ly=galactic_center_dist_ly)
    try:
        if position is None:
            system_entry = sector.add_system(system, system_config=config)
        else:
            system_entry = sector.add_preplaced_system(system, position, system_config=config)
    except ValueError:
        return None
    log.debug(f"Sector {sector.name!r}: planetary nebula {nebula.name!r} around new central star "
              f"{config.STAR_TYPE} {system.name!r}")
    return sector.add_phenomenon(nebula, "nebula", position=system_entry.position)


def add_star_hosted_nebulae(sector, args, chance_factor=1.0, hosts_out=None):
    """
    Grows emission and reflection nebulae around `sector`'s own hot
    stars (GEN.10, GEN.150): each system whose placed star (`StarSystem.placed_star`, the star its place drew, which a
    heavier companion does not displace, GEN.710) is a main-sequence
    star matching a `NEBULA_HOST_RULES` row rolls that row's chance (the star
    is still in its natal cloud), and on a hit gets a nebula of one of the row's
    classes centered on it, the classes the star is young enough for
    (`NEBULA_CLASS_MAX_AGE_MYR`). The nebula's age is the star's, and the star is
    the host that sizes it (`Nebula(host=, age_myr=)`).

    Args:
        sector (SpaceSector): The populated sector.
        args (argparse.Namespace): Parsed arguments.
        chance_factor (float): The sector profile's `"nebula"` factor (GEN.236, the cold gas fraction): the rules'
            chances are scaled by it, so a star in a gas-poor place stays bare.
        hosts_out (list, optional): Receives `(nebula, host_name, local_position_ly)` for each nebula made, the seed
            star of an emission region being where its members are placed from (GEN.321).

    Returns:
        list: The new nebulae's `SectorPhenomenonEntry` instances.
    """
    entries = []
    for system_entry in list(sector.entries):
        host = system_entry.star_system.placed_star
        if host is None:
            continue
        nebula = nebula_contents.hosted_nebula_of(host, chance_factor)
        if nebula is None:
            continue
        log.debug(f"Sector {sector.name!r}: {_spectral_code(host)} star {system_entry.star_system.name!r} "
                  f"lights class {nebula.nebula_class} nebula {nebula.name!r}")
        entries.append(sector.add_phenomenon(nebula, "nebula", position=system_entry.position))
        if hosts_out is not None:
            hosts_out.append((nebula, nebula_contents.host_name_of(host), system_entry.position))
    return entries


def _build_member(sector, args, local, mass_sol, age_yr, seed, galactic_center_dist_ly, label, member_tree=None):
    """
    One ordinary young system of `mass_sol` and `age_yr` at the sector-local point `local`, built from `seed`, or None when
    the point is not in `sector` or the sector's cell is closed to it. `member_tree` is the plan the system is built from
    (`stored_trees.drawn`), where a group left a tree; left out, the system draws its own companions.
    """
    if not sector.contains(local):
        return None
    config = run_system.build_system_config(args)
    config.POPULATION = "young"
    age_gy = age_yr / 1.0e9
    with draw.bound(seed):
        state = stellar_evolution.evolve_star(mass_sol, age_gy, population=config.POPULATION,
                                              galactic_center_dist_ly=galactic_center_dist_ly)
        params = stellar_evolution.star_params(mass_sol, age_gy, state)
        system = StarSystem(system_config=config, galactic_center_dist_ly=galactic_center_dist_ly,
                            primary_star_params=params, member_tree=member_tree)
    try:
        entry = sector.add_preplaced_system(system, local, system_config=config)
    except ValueError:
        return None
    entry.system_seed = seed
    fill_stages.count(sector, fill_stages.NEBULAE, label)
    return entry


def _group_plan(tree):
    """The plan a member's system is built from where its group had decayed: the hierarchy the group left (`None` where
    the member is in no group, or in a live one, whose mates are systems of their own and which builds as any star)."""
    if tree is None or nebula_contents.is_live_group(tree):
        return None
    return stored_trees.drawn(tree)


def add_group_stars(sector, args, member, tree, center_pc, galactic_center_dist_ly, label):
    """
    Builds the group mates of a member in a live young group (GEN.328, `nebula_contents.member_group`): each star of the
    group's tree other than the member is an ordinary young system of its own birth mass and the member's age, at the
    member's point plus the difference of the two offsets in the group (the member keeps the point it was given). The
    sector builds the mates whose points fall in it, so a group across a sector boundary is made once, each star by the
    sector that holds it. Every entry carries `group_seed` (the member's seed, which names the group), `group_decay_age_yr` and
    `group_offset_au`, and so does the member's own, which `_add_members` marks.

    Returns:
        list: The mates' `SectorSystemEntry` instances.
    """
    entries = []
    home = tree.members[0].offset_au
    for mate in tree.members[1:]:
        if mate.kind not in stored_trees.PAIR_STAR_KINDS:
            continue    # a brown dwarf or remnant is not a star a system is built around
        point_pc = tuple(member.position_pc[i] + (mate.offset_au[i] - home[i]) / physical_constants.AU_PER_PARSEC
                         for i in range(3))
        local = tuple(pc_to_ly(value) for value in galaxy_to_local_pc(center_pc, point_pc))
        entry = _build_member(sector, args, local, mate.initial_mass_solar, member.age_yr,
                              f"{member.seed}:group:{mate.member_no}", galactic_center_dist_ly, label)
        if entry is not None:
            mark_group_star(entry, member, mate)
            entries.append(entry)
    return entries


def mark_group_star(entry, member, group_member):
    """Records on `entry` the group its star is in: the group's name (the seed of the member that rolled it), its decay age and
    the star's offset from the group's centre, AU."""
    entry.group_seed = member.seed
    entry.group_decay_age_yr = group_member.decay_age_yr
    entry.group_offset_au = group_member.offset_au


def _add_members(sector, args, members, center_pc, galactic_center_dist_ly, label):
    """
    Builds the `members` (`nebula_contents.CloudMember`s) whose points fall inside `sector`, each as an ordinary young
    system at its point from its own seed. A member that finds the sector's cell closed to it is skipped, not the
    nebula. A massive member may be born in a young group (`nebula_contents.member_group`, GEN.328): its mates are built
    where they fall (`add_group_stars`), by the same draw in every sector that reaches the member, and where the group had
    decayed the member's system is the hierarchy it left.

    Returns:
        list: The new systems' `SectorSystemEntry` instances.
    """
    entries = []
    for member in members:
        tree = nebula_contents.member_group(member)
        live = nebula_contents.is_live_group(tree)
        if live:
            entries.extend(add_group_stars(sector, args, member, tree, center_pc, galactic_center_dist_ly, "group-mate"))
        local = tuple(pc_to_ly(value) for value in galaxy_to_local_pc(center_pc, member.position_pc))
        entry = _build_member(sector, args, local, member.mass_sol, member.age_yr, member.seed, galactic_center_dist_ly,
                              label, member_tree=_group_plan(tree))
        if entry is None:
            continue
        if live:
            mark_group_star(entry, member, tree.members[0])
        entries.append(entry)
    return entries


def add_cloud_protostars(sector, args, cloud_field, galactic_center_dist_ly=None):
    """
    Builds the protostars of the field clouds that reach `sector` (GEN.320): each Bok globule (P) and star-forming core
    (Q) draws its members from its own seed (`nebula_contents.cloud_members`), and the sector builds the ones whose
    points fall inside it, each as an ordinary young system at its point (`_add_members`).

    Args:
        sector (SpaceSector): The populated sector, with `field_nebulae` set.
        args (argparse.Namespace): Parsed arguments.
        cloud_field (tuple): `(galaxy_seed, shape, center_pc, reach_pc)` as `generate_sector_phenomena` took it.
        galactic_center_dist_ly (float, optional): As in `generate_sector`.

    Returns:
        list: The new systems' `SectorSystemEntry` instances.
    """
    galaxy_seed, _shape, center_pc, _reach = cloud_field
    entries = []
    for nebula, cloud_center in sector.field_nebulae:
        entries.extend(_add_members(sector, args, nebula_contents.cloud_members(galaxy_seed, nebula, cloud_center),
                                    center_pc, galactic_center_dist_ly, "protostar"))
    return entries


def add_emission_members(sector, args, cloud_field, own_hosts, neighbour_hosts=(), galactic_center_dist_ly=None):
    """
    Builds the ionizing members of the emission regions (D and E) that reach `sector` (GEN.321): those lit by a star in
    this sector (`own_hosts`, `add_star_hosted_nebulae`'s `hosts_out`, the star's point sector-local) and those lit by a
    pre-placed star in a neighbouring one (`neighbour_hosts`, galaxy-frame points). Each region's members are
    `nebula_contents.emission_members` of its seed star's point, the same in every sector that reaches it, and the sector
    builds the ones whose points fall inside it. The seed star itself is the host and already real.

    Args:
        sector (SpaceSector): The populated sector.
        args (argparse.Namespace): Parsed arguments.
        cloud_field (tuple): `(galaxy_seed, shape, center_pc, reach_pc)` as `generate_sector_phenomena` took it.
        own_hosts (list): `(nebula, host_name, local_position_ly)` for this sector's own hosted nebulae.
        neighbour_hosts (iterable): `(nebula, host_name, host_position_pc)` for the regions of other sectors' stars.
        galactic_center_dist_ly (float, optional): As in `generate_sector`.

    Returns:
        list: The new systems' `SectorSystemEntry` instances.
    """
    galaxy_seed, _shape, center_pc, _reach = cloud_field
    regions = [(nebula, name, local_to_galaxy_pc(center_pc, tuple(ly_to_pc(value) for value in local)))
               for nebula, name, local in own_hosts]
    regions.extend(neighbour_hosts)
    entries = []
    for nebula, host_name, host_position_pc in regions:
        members = nebula_contents.emission_members(galaxy_seed, nebula, host_name, host_position_pc)
        entries.extend(_add_members(sector, args, members, center_pc, galactic_center_dist_ly, "emission-member"))
    return entries


def flag_fast_stars(sector, galactic_center_dist_ly=None, scattered=False):
    """
    Marks some of `sector`'s systems as runaway or hypervelocity stars
    (`StarSystem.runaway_class`, `runaway_speed_kms` and a random
    `runaway_direction`, schema v37 and v61): an
    ordinary generated system moving unusually fast, not a separate
    phenomenon. Each system rolls a hypervelocity chance first -- the
    "hypervelocity-star" rate per star scaled by
    `(scales.hypervelocity_reference_radius_pc() / r)^2`, since the central black hole ejects
    them -- then the runaway chance (about 1.5%).

    Args:
        sector (SpaceSector): The populated sector, changed in place.
        galactic_center_dist_ly (float, optional): The sector's distance
            from the galactic center; `None` uses `scales.default_galactocentric_ly()`.
        scattered (bool): The galaxy's phenomenon scatter ran (GEN.100):
            its hypervelocity stars are real rows already, so none is
            flagged here. Runaway stars still are.

    Returns:
        int: How many systems were flagged.
    """
    if galactic_center_dist_ly is None:
        galactic_center_dist_ly = scales.default_galactocentric_ly()
    radius_pc = max(ly_to_pc(galactic_center_dist_ly), 1.0)
    hvs_chance = 0.0 if scattered else min(
        1.0, program_constants.phenomenon_rate_per_star("hypervelocity-star")
        * (scales.hypervelocity_reference_radius_pc() / radius_pc) ** 2)
    runaway_chance = program_constants.phenomenon_rate_per_star("runaway-star")
    flagged = 0
    for entry in sector.entries:
        system = entry.star_system
        if getattr(system, "runaway_class", None) is not None or getattr(system, "runaway_decided", False):
            continue  # a scattered hypervelocity star is flagged already, a pre-placed O or B star was decided by its seed
        if draw.random() < hvs_chance:
            system.runaway_class = "hypervelocity"
            system.runaway_speed_kms = log_uniform(*program_constants.HYPERVELOCITY_STAR_SPEED_RANGE_KMS)
            system.runaway_direction = random_unit_vector()
        elif draw.random() < runaway_chance:
            system.runaway_class = "runaway"
            system.runaway_speed_kms = log_uniform(*program_constants.RUNAWAY_STAR_SPEED_RANGE_KMS)
            system.runaway_direction = random_unit_vector()
            system.runaway_origin = multistar.runaway_origin_at(draw.random())
        else:
            continue
        flagged += 1
    return flagged


NUCLEUS_ADDRESS = phenomenon_scatter.NUCLEUS_ADDRESS
"""tuple: The `(shell, band, wedge)` of the one sector per galaxy that rolls
for an active nucleus (see `add_galactic_nucleus`); in a galaxy whose
phenomenon scatter ran, the nucleus is a scatter row instead."""


def add_galactic_nucleus(sector, args, galactic_center_dist_ly, mass_solar=None):
    """
    Adds the galaxy's central supermassive black hole to `sector` at the
    galactic center: an active `Quasar` with
    `program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE`, otherwise a quiescent
    supermassive `BlackHole` (like Sagittarius A*), so every galaxy has one.

    There is only ever one nucleus, and only at the origin.
    `generate_and_save_sector_at` calls this for exactly one sector per
    galaxy, `NUCLEUS_ADDRESS` (shell 0, the one core cell, which is
    centered on the origin, so picking it keeps it to a single roll).

    The core sector's center is the galactic origin, so the nucleus sits at
    `(-galactic_center_dist_ly, 0, 0)` in the sector's own frame (a
    distance of 0 there); `store.insert_sector` converts that back to the
    galaxy origin.

    Args:
        sector (SpaceSector): The core sector, already populated.
        args (argparse.Namespace): Parsed arguments (none is read here).
        galactic_center_dist_ly (float): The sector center's distance from
            the galactic center, in light-years.
        mass_solar (float, optional): The central black hole's mass (`density.central_black_hole_mass_solar`,
            GEN.402); without a stored galaxy, a galaxy is drawn for it.

    Returns:
        SectorPhenomenonEntry: The quasar's or the black hole's entry.
    """
    nucleus_config = SystemConfig()
    position = (-galactic_center_dist_ly, 0.0, 0.0)
    with fill_stages.stage(sector, fill_stages.FLAGS, "nucleus"):
        if draw.random() < program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE:
            fill_stages.count(sector, fill_stages.FLAGS, "quasar")
            return sector.add_phenomenon(Quasar(nucleus_config, black_hole_mass_solar=mass_solar), "quasar",
                                         position=position)
        log.debug(f"Sector {sector.name!r}: galactic nucleus is quiescent (supermassive black hole, no quasar)")
        black_hole = BlackHole(nucleus_config, mass_class="supermassive", mass_solar=mass_solar)
        fill_stages.count(sector, fill_stages.FLAGS, "black-hole")
        return sector.add_phenomenon(black_hole, "black-hole", position=position)


def _seeded_build(seed, build):
    """`build()` on the draw stream bound to `seed`, the sector's stream put
    back afterward: an object built from its scatter row's seed is the same
    whenever and wherever its sector is filled."""
    with draw.bound(seed):
        return build()


def add_scattered_phenomena(sector, args, fill):
    """
    Builds a galaxy-placed sector's scattered phenomena (GEN.100): every
    unbuilt `phenomenon_scatter` row in the sector's cell becomes the
    object it stands for, from the row's own seed, at the row's stored
    point -- a black hole, neutron star, supernova remnant, planetary
    nebula (around its own new central star), hypervelocity star (a real
    system moving outward at its stored velocity), or the galaxy's nucleus.
    A sector places these before its own systems, so each of those keeps
    clear of a black hole's or neutron star's Hill sphere like any other.

    Returns:
        list: The new `SectorPhenomenonEntry` instances (a hypervelocity
            star adds a system entry instead).
    """
    new_entries = []
    for row in fill.phenomenon_rows:
        position = brightStars.local_position_ly(row, fill.center_pc)
        radius_pc = math.sqrt(sum((row[axis] / brightStars.MPC_PER_PC) ** 2
                                  for axis in ("position_x_mpc", "position_y_mpc", "position_z_mpc")))
        dist_ly = pc_to_ly(radius_pc)
        cut = fill.below_cut.min_mass_solar if fill.below_cut is not None else None
        type_ia_share = getattr(getattr(fill, "profile", None), "type_ia_share", None)
        entry = _seeded_build(row["seed"], lambda: _build_scattered(sector, args, row, position, dist_ly, cut,
                                                                    type_ia_share, _bubble_density(fill), fill))
        if entry is None:
            continue
        entry.scatter_id = row["id"]
        entry.system_seed = row["seed"]
        fill_stages.count(sector, fill_stages.SCATTERED, row["kind"])
        if hasattr(entry, "phenomenon_type"):
            new_entries.append(entry)
    return new_entries


def _bubble_density(fill):
    """The interior density, cm^-3, of the superbubble the fill's sector profile places the sector in (GEN.298), or
    `None` outside a bubble or without a profile (a standalone sector)."""
    bubble = getattr(getattr(getattr(fill, "profile", None), "gas", None), "superbubble", None)
    return None if bubble is None else bubble.ambient_density_cm3


def _build_scattered(sector, args, row, position, dist_ly, min_mass_solar=None, type_ia_share=None,
                     ambient_density_cm3=None, fill=None):
    """One scatter row built into `sector` at `position` (sector-local
    light-years): its entry, or `None` when a planetary nebula's central
    system found no room. A neutron star's or black hole's mass is drawn
    at or above the scatter's cut `min_mass_solar` (GEN.166; `None`: the
    whole range). A supernova remnant is Type Ia with the share `type_ia_share` (GEN.319; `None`: the solar
    circle's). A supernova remnant is sized by `ambient_density_cm3` where the sector lies in a superbubble (GEN.298;
    `None`: its own class density). A neutron star or black hole takes its age from the population the sector's
    `fill` gives it (GEN.553)."""
    kind = row["kind"]
    config = SystemConfig()
    config.TYPE_IA_SHARE = type_ia_share
    config.AMBIENT_DENSITY_CM3 = ambient_density_cm3
    if kind in POPULATION_AGED_KINDS and fill is not None:
        fill.apply_population(config)
    if kind == "hypervelocity-star":
        return _add_hypervelocity_star(sector, args, row, position, dist_ly)
    if kind == "planetary-nebula":
        nebula = Nebula(config, nebula_type="planetary")
        return _add_planetary_nebula(sector, args, nebula, dist_ly, position=position)
    if kind == "quasar":
        return sector.add_phenomenon(Quasar(config, black_hole_mass_solar=row["mass_solar"]), "quasar",
                                     position=position)
    if kind == "black-hole":
        black_hole = BlackHole(config, galactic_center_dist_ly=dist_ly, mass_class=row["subtype"],
                               mass_solar=row["mass_solar"],
                               mass_range=phenomenon_scatter.mass_range(kind, row["subtype"], min_mass_solar, True))
        return sector.add_phenomenon(black_hole, "black-hole", position=position)
    if kind == "neutron-star":
        neutron_star = NeutronStar(config, galactic_center_dist_ly=dist_ly, mass_solar=row["mass_solar"],
                                   mass_range=phenomenon_scatter.mass_range(kind, row["subtype"], min_mass_solar, True),
                                   magnetar=row["subtype"] == phenomenon_scatter.MAGNETAR_SUBTYPE)
        return sector.add_phenomenon(neutron_star, "neutron-star", position=position)
    if kind == "supernova-remnant":
        return sector.add_phenomenon(SupernovaRemnant(config), "supernova-remnant", position=position)
    raise ValueError(f"Unknown scattered phenomenon kind: {kind!r}")


def add_below_cut_remnants(sector, args, fill, galactic_center_dist_ly=None):
    """
    The neutron stars and black holes the phenomenon scatter left below its
    mass cut (GEN.168): `phenomenon_scatter.below_cut_draws` counts them on
    the sector's own stream, and each is built with its mass below the cut
    and placed like a rolled one (anywhere in the cube clear of the Hill
    spheres already there), on that same stream.

    Returns:
        list: The new `SectorPhenomenonEntry` instances (one that found no
            room is skipped, as `generate_sector_phenomena` does).
    """
    below_cut = fill.below_cut
    draws, rng = phenomenon_scatter.below_cut_draws(below_cut.address, fill.center_pc, fill.shape,
                                                    below_cut.expected_stars, below_cut.min_mass_solar,
                                                    below_cut.seed)
    entries = []
    made = 0
    with draw.bound(rng):
        for kind, subtype, mass_range, count in draws:
            for _ in range(count):
                config = fill.apply_population(SystemConfig())
                if kind == "black-hole":
                    remnant = BlackHole(config, galactic_center_dist_ly=galactic_center_dist_ly,
                                        mass_class=subtype, mass_range=mass_range)
                else:
                    remnant = NeutronStar(config, galactic_center_dist_ly=galactic_center_dist_ly,
                                          mass_range=mass_range)
                made += 1
                try:
                    entry = sector.add_phenomenon(remnant, kind)
                except ValueError:
                    continue  # no room clear of the Hill spheres already placed
                # The remnant's system roll (stage C) is its own: its place in the draw names it.
                entry.system_seed = f"{below_cut.seed}:{made}"
                entries.append(entry)
                fill_stages.count(sector, fill_stages.BELOW_CUT, kind)
    return entries


def add_ejecta(sector, args, fill, galactic_center_dist_ly=None):
    """
    The neutron stars and black holes the galaxy ejected that stream through this sector (GEN.238): per kind, a Poisson
    count with the mean of the sector's profile (`SectorProfile.ejecta`, `ejecta_field`), each remnant built like any
    other (a stellar-mass black hole), placed like a rolled one (anywhere in the cube clear of the Hill spheres already
    there), and given a velocity of its own (`ejecta_field.draw_velocity_kms`: the galaxy's kick law for a neutron star) in place of an orbit. They carry no
    planets: a kick that unbinds a star strips them. A sector with no profile (a standalone one) has none.

    Returns:
        list: The new `SectorPhenomenonEntry` instances (one that found no room is skipped).
    """
    profile = getattr(fill, "profile", None)
    expected = getattr(profile, "ejecta", None) or {}
    entries = []
    for kind in ejecta_field.KINDS:
        mean = expected.get(kind, 0.0)
        if mean <= 0.0:
            continue
        with fill_stages.stage(sector, fill_stages.EJECTA, kind):
            for _ in range(_sample_poisson_count(mean)):
                config = fill.apply_population(SystemConfig())
                if kind == "black-hole":
                    remnant = BlackHole(config, galactic_center_dist_ly=galactic_center_dist_ly, mass_class="stellar")
                else:
                    remnant = NeutronStar(config, galactic_center_dist_ly=galactic_center_dist_ly)
                remnant.eject(ejecta_field.draw_velocity_kms(kind, getattr(profile, "shape", None)))
                try:
                    entries.append(sector.add_phenomenon(remnant, kind))
                except ValueError:
                    continue  # no room clear of the Hill spheres already placed
                fill_stages.count(sector, fill_stages.EJECTA, kind)
    return entries


def _add_hypervelocity_star(sector, args, row, position, dist_ly):
    """A scattered hypervelocity star built as an ordinary young system at
    its stored point, flagged with its speed and direction (the stored
    galactic velocity). Returns its system entry."""
    config = run_system.build_system_config(args)
    config.POPULATION = "young"
    system = StarSystem(system_config=config, galactic_center_dist_ly=dist_ly)
    velocity = (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])
    speed = math.sqrt(sum(component * component for component in velocity))
    system.runaway_class = "hypervelocity"
    system.runaway_speed_kms = speed
    # Its point holds at the time the scatter drew it; the orbit update flies it on from there (GEN.137).
    system.epoch_unix = row.get("epoch_unix")
    system.runaway_direction = tuple(component / speed for component in velocity) if speed > 0.0 else (1.0, 0.0, 0.0)
    return sector.add_preplaced_system(system, position, system_config=config)


def _flag_scattered_runaway(system, row):
    """
    Marks the system built around a pre-placed O or B star as a runaway when its seed says it is one (GEN.329,
    `scatter_companions.runaway_draw`: 20% of O stars and 7.5% of B stars, two thirds by the supernova channel): the same
    answer the scatter used for the star's tree, with its speed and direction drawn from a stream of the row's seed. The
    system is marked decided either way, so the fill's field roll (`flag_fast_stars`) leaves it alone.
    """
    if row["yerkes_class"] == "VII":
        return
    system.runaway_decided = True
    runaway = scatter_companions.runaway_draw(row["seed"], row["initial_mass_sol"])
    if runaway is None:
        return
    with draw.bound(draw.Stream(f"{row['seed']}:runaway")):
        system.runaway_class = "runaway"
        system.runaway_speed_kms = log_uniform(*program_constants.RUNAWAY_STAR_SPEED_RANGE_KMS)
        system.runaway_direction = random_unit_vector()
    system.runaway_origin = runaway.origin


def _add_preplaced_systems(sector, args, fill, galactic_center_dist_ly):
    """Builds a full system around each of the sector's pre-placed bright
    stars (`fill.bright_rows`) and places it first, at its stored point.
    Each one's companion, planets and moons come from its own stored
    `seed`, without disturbing the sector's draw stream. A row with members is built from the tree it stored
    (`fill.stored_plan`, GEN.315), a row without is one star: neither draws a companion of its own."""
    for row in fill.bright_rows:
        cfg = run_system.build_system_config(args)
        cfg.POPULATION = row["population"]
        cfg.PLACE_TERMS = fill.place_terms
        with draw.bound(row["seed"]):
            system = StarSystem(system_config=cfg, galactic_center_dist_ly=galactic_center_dist_ly,
                                primary_star_params=brightStars.star_params(row),
                                member_tree=fill.stored_plan(row["id"]))
        _flag_scattered_runaway(system, row)
        entry = sector.add_preplaced_system(system, brightStars.local_position_ly(row, fill.center_pc),
                                            system_config=cfg)
        entry.bright_star_id = row["id"]
        fill_stages.count(sector, fill_stages.PRE_PLACED, "star-system")


def generate_sector(args, galactic_center_dist_ly=None, cell=None, fill=None, cloud_field=None):
    """
    Builds a fully populated `SpaceSector` from parsed args, without
    rendering, printing, or saving anything -- the shared core `run_sector`
    and the `galaxy` subcommand both build on.

    Args:
        args (argparse.Namespace): Parsed arguments from the `sector`
            subcommand (or an equivalently-shaped namespace the `galaxy`
            subcommand builds itself -- see
            `add_shared_generation_options`/`validate_shared_generation_args`
            for the option surface this function actually reads, via
            `iter_sector_configs`).
        galactic_center_dist_ly (float, optional): This sector's distance
            from the galactic center, in light-years, threaded down into
            every generated system's Hill-sphere calculation (see
            `Star.calculate_system_perimeter`'s docstring). `None` (the
            default) falls back to the radius that holds half of the galaxy's disc stars
            (`galaxy.scales.default_galactocentric_ly`) -- the `sector` subcommand (no galaxy context) always
            calls this with the default, so it keeps producing "unplaced"
            sectors.
        cell (galaxy.spherical.SectorCell, optional): The sector's real
            shell, band and wedge cell, in light-years, for a galaxy-placed sector --
            systems and phenomena are then placed inside it rather than in
            a cube (see `SpaceSector.cell`).
        fill (brightStars.FillContext, optional): For a galaxy-placed
            sector: its pre-placed bright stars are built and placed
            first, and every other system takes its age from the stellar
            population mix there and stays below the bright-star
            threshold; a `--density` count shrinks by the bright stars'
            share so the expected total is unchanged.
        cloud_field (tuple, optional): See `generate_sector_phenomena`.

    Returns:
        tuple: `(sector_name, SpaceSector)` -- `sector_name` is
              `args.sector_name` or a freshly generated one (see
              `generate_sector_name`); the `SpaceSector` has every system
              that fit added (`SpaceSector.add_system`'s Hill-sphere-based
              random placement -- see the "Capacity" note below for what
              happens when not all of them do), plus a realistically
              sparse population of exotic phenomena (see
              `generate_sector_phenomena`). When `args.density` drove the
              system count (see below) and both that Poisson draw and
              `generate_sector_phenomena`'s own independent draws came
              back completely empty, one system is force-added anyway --
              see the "guaranteed non-empty" note below.

    Capacity
    --------
    A sector's cube can only physically hold so many systems before every
    point left is within some existing system's Hill sphere -- real space
    is finite, and `SpaceSector.add_system`'s random placement
    (`_random_position`) raises `ValueError` once it can't find room for
    one more within `program_constants.SECTOR_MAX_PLACEMENT_ATTEMPTS`
    tries. `--num-systems`/`--density` (and, compounding, `--min-habitable`
    forcing extra large stars with their own larger Hill spheres) can ask
    for more systems than a `program_constants.DEFAULT_SECTOR_EDGE_LY`
    cube can hold at realistic stellar spacing. Rather than let that
    `ValueError` propagate (previously discarding every system already
    generated for this sector, including the expensive planet/moon
    generation behind each one) or keep burning full `StarSystem`
    generations on placements that are no longer going to fit, this loop
    stops as soon as one system can't be placed and returns the sector as
    it stands -- with fewer systems than requested, not zero. See
    `sector_generation_summary_lines`'s "actual vs. expected" density
    line, which already accounts for this (it was previously only ever
    exercised by `--density`'s own Poisson undersampling).
    """
    sector_name = args.sector_name or generate_sector_name()
    sector = SpaceSector(name=sector_name, cell=cell)
    fill_stages.begin(sector)       # every stage below draws from its own stream (GEN.251)
    if cloud_field is not None:     # GEN.357: the clouds are known before the first system is placed
        cloudContext.begin(sector, cloud_field, nebula_field.clouds_reaching(*cloud_field))
    scattered = fill is not None and fill.phenomena_scattered
    if fill is not None:
        with fill_stages.stage(sector, fill_stages.PRE_PLACED):
            _add_preplaced_systems(sector, args, fill, galactic_center_dist_ly)
    scattered_entries = []
    if scattered:
        with fill_stages.stage(sector, fill_stages.SCATTERED):
            scattered_entries = add_scattered_phenomena(sector, args, fill)
    with fill_stages.stage(sector, fill_stages.SYSTEM_PASS):
        system_pass.run(sector, args, scattered_entries, galactic_center_dist_ly)

    with fill_stages.stage(sector, fill_stages.SYSTEMS):
        # `--density` resolves to a concrete system count per sector (this
        # sector's own volume, sampled fresh each call) rather than once at
        # parse time -- this is what lets each sector under `--num-sectors`/
        # the `galaxy` subcommand vary independently instead of sharing one
        # fixed count. A copy avoids mutating the caller's shared `args`
        # namespace, since this function runs once per sector.
        density_driven = args.density is not None
        if density_driven:
            args = copy.copy(args)
            # An astronomical --density (1e308) can overflow the mean to inf;
            # any count that large just fills the sector (see "Capacity").
            mean = min(sector.expected_system_count() * args.density, sys.float_info.max)
            if fill is not None:
                mean *= 1.0 - fill.bright_share()
            args.num_systems = _sample_poisson_count(mean)

        total = args.num_systems
        for i, cfg in enumerate(iter_sector_configs(args)):
            if fill is not None:
                fill.apply(cfg)
            with log.timed_phase(f"generate system {i + 1}/{total}"):
                system = StarSystem(system_config=cfg, galactic_center_dist_ly=galactic_center_dist_ly)

            try:
                with log.timed_phase(f"place system {i + 1}/{total}"):
                    sector.add_system(system, system_config=cfg)
                fill_stages.count(sector, fill_stages.SYSTEMS, "star-system")
            except ValueError:
                # No room left for another system's Hill sphere in this
                # sector's cube -- see this function's own "Capacity"
                # docstring note. Every following config would almost
                # certainly fail the same way (this sector only gets fuller
                # from here), so stop generating and placing altogether
                # rather than pay for `total - i - 1` more full
                # StarSystem generations just to discard them too.
                log.normal(
                    f"Sector '{sector_name}': ran out of room after placing {len(sector.entries)} of "
                    f"{total} requested systems -- the {sector.edge_ly:.1f} ly cube has no space left "
                    f"that clears every already-placed system's Hill sphere. Returning the sector as-is "
                    f"rather than the full requested count."
                )
                break

    with log.timed_phase("generate_sector_phenomena"):
        generate_sector_phenomena(sector, args, galactic_center_dist_ly=galactic_center_dist_ly,
                                  cloud_field=cloud_field, scattered=scattered,
                                  kind_factors=getattr(getattr(fill, "profile", None), "kind_factors", None),
                                  type_ia_share=getattr(getattr(fill, "profile", None), "type_ia_share", None),
                                  ambient_density_cm3=_bubble_density(fill),
                                  neighbour_hosts=getattr(fill, "emission_hosts", ()), fill=fill)
        if scattered and fill.below_cut is not None:
            with fill_stages.stage(sector, fill_stages.BELOW_CUT):
                below_cut_entries = add_below_cut_remnants(sector, args, fill, galactic_center_dist_ly)
            with fill_stages.stage(sector, fill_stages.SYSTEM_PASS, "below-cut"):
                system_pass.run(sector, args, below_cut_entries, galactic_center_dist_ly)
        add_ejecta(sector, args, fill, galactic_center_dist_ly)
        with fill_stages.stage(sector, fill_stages.FLAGS, "fast-stars"):
            flagged = flag_fast_stars(sector, galactic_center_dist_ly=galactic_center_dist_ly, scattered=scattered)
        fill_stages.count(sector, fill_stages.FLAGS, "fast-star", flagged)

    capped = fill is not None and fill.capped
    if density_driven and not sector.entries and not sector.phenomena and not capped:
        # Guaranteed non-empty: a sector's own Poisson draws
        # (system count here, each phenomenon type in generate_sector_phenomena)
        # are independent, so all of them landing on zero simultaneously is
        # a real, expected outcome at low means (e.g. ~13% at mean 2) --
        # but a sector with literally nothing in it isn't useful to anyone
        # visiting it, so force exactly one system rather than leave it
        # empty. Only applies when the count came from --density (an
        # explicit --num-systems, including 0, is a deliberate request
        # this never overrides).
        with fill_stages.stage(sector, fill_stages.SYSTEMS, "fallback"):
            fallback_config = run_system.build_system_config(args)
            if fill is not None:
                fill.apply(fallback_config)
            fallback_system = StarSystem(system_config=fallback_config,
                                         galactic_center_dist_ly=galactic_center_dist_ly)
            sector.add_system(fallback_system, system_config=fallback_config)
        fill_stages.count(sector, fill_stages.SYSTEMS, "star-system")

    for line in fill_stages.describe(sector):
        log.debug(f"Sector {sector_name!r} fill, {line}")
    return sector_name, sector


def _summary_star_label(star):
    """`(label, plural)` for one system's star in the sector summary (UX.34)."""
    return star_label(getattr(star, "type", None), getattr(star, "yerkes_class", None))


def _log_summary(text):
    """Logs a multi-line summary one line per record, so each line gets
    the debug log's timestamp, process and source prefix (OPS.9)."""
    for line in text.splitlines():
        log.normal(line)


def sector_generation_summary_lines(sector, args_used):
    """
    Builds the per-sector status lines every sector-generating command
    (`sector`, `galaxy`'s three modes) prints right after a sector is
    saved: how many star systems of each spectral class and phenomena of
    each type were actually generated, plus how the sector's actual star
    count compares to what its own local density predicted -- so a run's
    output is enough on its own to sanity-check what came out of it,
    without a separate database query.

    Args:
        sector (SpaceSector): The freshly generated (and, by the time a
            caller has this, already saved -- so `sector.name` is final)
            sector.
        args_used (argparse.Namespace): The args actually passed to
            `generate_sector` to build this particular sector --
            `.density`/`.num_systems` reflect whichever one drove it
            (for `galaxy` mode, this is `_BatchDensity.resolve`'s own
            per-sector result, not necessarily the top-level parsed args).

    Returns:
        str: Two or three newline-joined, already-indented lines: systems
            by spectral class, phenomena by type (omitted when the sector
            has none), and actual vs. expected star density (`1.0` =
            real local stellar density, matching `--density`'s own
            convention).
    """
    system_types = Counter(_summary_star_label(entry.star_system.star) for entry in sector.entries)
    phenomenon_types = Counter(entry.phenomenon_type for entry in sector.phenomena)

    e_value = sector.expected_system_count()
    if args_used.density is not None:
        expected_density = args_used.density
    elif e_value > 0:
        expected_density = (args_used.num_systems or 0) / e_value
    else:
        expected_density = 0.0
    actual_density = (len(sector.entries) / e_value) if e_value > 0 else 0.0

    lines = []
    if system_types:
        types_str = describe_types(system_types)
    else:
        types_str = "none"
    lines.append(f"    Systems: {types_str}")

    if phenomenon_types:
        phenomena_str = ", ".join(
            f"{count} {run_phenomenon.TYPE_LABELS[ptype]}" for ptype, count in sorted(phenomenon_types.items())
        )
        lines.append(f"    Phenomena: {phenomena_str}")

    stars = sum(len(getattr(entry.star_system, "stars", None) or [None]) for entry in sector.entries)
    planets = sum(len(entry.star_system.all_planets) if hasattr(entry.star_system, "all_planets") else
                  len(getattr(entry.star_system, "planets", None) or ())
                  + len(getattr(entry.star_system, "secondary_planets", None) or ()) for entry in sector.entries)
    lines.append(f"    Totals: {len(sector.entries):,} star systems, {stars:,} stars, {planets:,} planets, "
                 f"{len(sector.phenomena):,} phenomena")
    lines.append(
        f"    Star density: actual {actual_density:.2f}x local, expected {expected_density:.2f}x local"
    )

    return "\n".join(lines)


class UnplacedSectorsFollower(followed.Follower):
    """
    The sectors of a `sector` run as its one stage (PERF.119, PERF.202): the `generate-unplaced-sector` command's units
    generate and save them, and this side logs each as it is saved and keeps the run's counts.
    """

    key = "unplaced"

    def __init__(self, args):
        super().__init__(args)
        self.label = commands.build_steps("sector", args)[0].label

    def params(self):
        return commands.build_steps("sector", self.args)[0].params

    def on_result(self, unit, result):
        args = self.args
        mysql_config = store.mysql_config_from_args(args)
        elsewhere = result["pid"] != os.getpid()     # a worker's own RUN_COUNTS die with it; the run's are here
        for saved in result["sectors"]:
            if elsewhere:
                run_common.RUN_COUNTS["sectors"] += 1
                run_common.RUN_COUNTS["systems"] += saved["systems"]
                run_common.RUN_COUNTS["phenomena"] += saved["phenomena"]
            run_common._record_sector(args, saved, saved["seconds"])
            phenomena_note = f", {saved['phenomena']} phenomena" if saved["phenomena"] else ""
            log.normal(
                f"Saved sector '{saved['name']}' to the database (sector_id={saved['sector_id']}, "
                f"{saved['systems']} systems{phenomena_note}, "
                f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port})."
            )
            _log_summary(saved["summary"])


def run_sector(args):
    """
    Generates `args.num_sectors` sectors (1 by default), each one a
    fully populated `SpaceSector` (see `generate_sector` -- one
    `SystemConfig` per system via `iter_sector_configs`, a full
    `StarSystem` from each, all placed in the sector via
    `SpaceSector.add_system`'s Hill-sphere-based placement) saved into
    the same database. No `galactic_center_dist_ly` is passed to
    `generate_sector` here -- the `sector` subcommand has no galaxy
    context (see the `galaxy` subcommand for that), so every sector it
    produces is "unplaced", regardless of `--num-sectors`.

    The run is one command through the wrapper (PERF.119, `generate-unplaced-sector`): the wrapper cuts the sectors
    into units, a worker generates and saves each, and this process logs each as it is saved. Each sector draws from
    a stream of its own, seeded from the run's and the sector's number, so the sectors are the same at any worker
    count.

    Each sector is saved to the database; only a short status line and
    summary per saved sector is printed, so a large `--num-sectors` run
    doesn't flood the console with rendered text nobody asked to see.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "sector"`).
    """
    mysql_config = store.mysql_config_from_args(args)
    try:
        run_common._check_estimate(args, [args] * args.num_sectors, f"{args.num_sectors} unplaced sector(s)")
    except run_common._EstimateOnly as exc:
        run_common._print_estimate(exc)
        return

    args.database = wrapper.database_args(mysql_config)
    args.run_seed = draw.getrandbits(128)
    follower = UnplacedSectorsFollower(args)
    try:
        followed.run_followed(args, [follower], follower.label, galaxy=False)
    except (redisqueue.Unavailable, runrecord.Busy, units.RunStopped) as exc:
        log.error(f"error: the sectors could not be generated: {exc}")
        raise SystemExit(1) from exc

    if args.num_sectors > 1:
        log.normal(f"Generated {args.num_sectors} sectors.")


UNPLACED_UNITS_MAX = 256
"""int: The most units a `sector` run is cut into: its sectors in static chunks of consecutive numbers, the same at any
worker count."""


def expand_unplaced_sectors(params):
    """
    The cut of a `sector` run (PERF.119; `generate-unplaced-sector`): `params["count"]` sectors numbered from 0, in about
    `UNPLACED_UNITS_MAX` static chunks of consecutive numbers. A unit's cost is its sectors.
    """
    count = params["count"]
    size = max(1, -(-count // UNPLACED_UNITS_MAX))
    units = []
    for start in range(0, count, size):
        numbers = list(range(start, min(start + size, count)))
        units.append(({"database": params["database"], "args": params["args"], "run_seed": params["run_seed"],
                       "numbers": numbers}, float(len(numbers))))
    return commands.Expansion(units, {"sectors": count})


def _task_seed(run_seed, number):
    """The seed of unplaced sector `number` of a run: the run's seed and the number, hashed, so it does not depend on the
    worker count or the order the units finish in."""
    digest = hashlib.sha256(f"{run_seed}:unplaced-{number}".encode("utf-8")).digest()
    return int.from_bytes(digest[:16], "big")


def unplaced_sectors_unit(unit_args, report):
    """
    The unit of a `sector` run: generates and saves each sector of the chunk (`_unplaced_sector_task`), each on its
    own stream (`_task_seed`).

    Returns:
        dict: `pid` (this process: a unit run by the process that started the run has counted its sectors there
            already) and `sectors`, one entry per sector in order: `number`, `sector_id`, `name`, `systems`, `stars`,
            `density`, `phenomena`, `summary` and `seconds`.
    """
    from planetgen.generation import run_galaxy
    args = run_galaxy._sector_args_of(unit_args)
    numbers = unit_args["numbers"]
    done = []
    for number in numbers:
        started = time.monotonic()
        with draw.bound(_task_seed(unit_args["run_seed"], number)):
            saved = _unplaced_sector_task(args)
        saved.update(number=number, seconds=time.monotonic() - started)
        done.append(saved)
        report(len(done), len(numbers))
    return {"sectors": done, "pid": os.getpid()}


def _forced_args(args, count):
    """`args` for the forced fill (GEN.180): `--min-habitable COUNT`, which forces a habitable world into that many of the
    sector's systems (and gives the sector that many systems at least, with a warning rather than a refusal even under
    `--strict`, as the draws have already failed)."""
    forced = copy.copy(args)
    forced.min_habitable = max(int(getattr(args, "min_habitable", 0) or 0), count)
    forced.strict = False
    return forced


def generate_directed_sector(args, galaxy_seed=None, address=None, **kwargs):
    """
    `generate_sector`, redrawn until the run's `--directive` minimums hold
    (GEN.96; see `directives.generate`). With none given it is one plain
    draw. A directive the closest draw still misses is logged as a warning.

    Returns:
        tuple: `(sector_name, SpaceSector)`, as `generate_sector`.
    """
    directive = directives.parse(getattr(args, "directive", None))
    if not directive:
        return generate_sector(args, **kwargs)
    sector, outcome = directives.generate(
        directive, lambda: generate_sector(args, **kwargs)[1], galaxy_seed=galaxy_seed, address=address,
        max_attempts=getattr(args, "directive_attempts", None) or directives.DEFAULT_ATTEMPTS,
        build_forced=lambda count: generate_sector(_forced_args(args, count), **kwargs)[1])
    where = f" {galaxySeed.address_text(address)}" if address is not None else ""
    log.normal(f"Sector{where}: {outcome.describe()}.")
    return sector.name, sector


def _unplaced_sector_task(args):
    """One `sector` subcommand sector, generated and saved -- a work queue
    task (see `_fill_sector_task`)."""
    _sector_name, sector = generate_directed_sector(args)
    sector_id = store.save_sector(sector, config=store.mysql_config_from_args(args),
                                  run_id=getattr(args, "generation_run_id", None))
    run_common._count_sector(sector)
    return {
        "sector_id": sector_id, "name": sector.name, "systems": len(sector.entries),
        "stars": sector_star_count(sector), "density": run_common._sector_density(args),
        "phenomena": len(sector.phenomena), "summary": sector_generation_summary_lines(sector, args),
    }
