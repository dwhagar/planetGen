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
import random
import re
import sys
from collections import Counter

from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.galaxy import nebula_field
from planetgen.physics import constants
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula, choose_weighted_class
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.galaxy.sector import SpaceSector, _sample_poisson_count
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.system import StarSystem
from planetgen.names.wordsalad import generate_sector_name
from planetgen.physics.units import ly_to_pc
from planetgen.util.random import log_uniform
from planetgen.generation import run_common
from planetgen.generation import run_phenomenon
from planetgen.generation import run_population
from planetgen.generation import run_system


SECTOR_PHENOMENON_KINDS = (
    ("black-hole", "black-hole",
     lambda config, dist_ly: BlackHole(config, galactic_center_dist_ly=dist_ly)),
    ("neutron-star", "neutron-star",
     lambda config, dist_ly: NeutronStar(config, galactic_center_dist_ly=dist_ly)),
    ("planetary-nebula", "nebula", lambda config, dist_ly: Nebula(config, nebula_type="planetary")),
    ("molecular-cloud", "nebula", lambda config, dist_ly: Nebula(config, nebula_type="dark")),
    ("supernova-remnant", "supernova-remnant", lambda config, dist_ly: SupernovaRemnant(config)),
    ("rogue-planet", "rogue-planet", lambda config, dist_ly: RoguePlanet(config)),
    ("brown-dwarf", "rogue-planet", lambda config, dist_ly: RoguePlanet(config, mass_bin="brown-dwarf")),
    ("comet", "comet", lambda config, dist_ly: InterstellarComet(config)),
    ("asteroid-field", "asteroid-field", lambda config, dist_ly: AsteroidField(config)),
)
"""tuple: `(program_constants.PHENOMENON_DENSITY_PC3 key, phenomenon type,
factory)` for each kind `generate_sector_phenomena` rolls. A factory takes
`(config, galactic_center_dist_ly)`; only black holes and neutron stars
use the distance (their Hill sphere and galactic orbit). A planetary
nebula is generated around its own new hot white dwarf system
(`_add_planetary_nebula`). Runaway and hypervelocity stars are flags on
ordinary systems (`flag_fast_stars`), and emission and reflection nebulae
grow around a sector's own hot stars (`add_star_hosted_nebulae`). In a
galaxy-placed sector the "molecular-cloud" roll gives way to the galaxy's
cloud field (`nebulaField`, GEN.47). Diffuse gas (classes A-B) fills about
half the disk, so it is the background, not generated."""


def build_sector_configs(args):
    """
    Builds one `SystemConfig` per system in the sector, sharing the same
    value options across all of them (via `build_system_config`, so this
    can never silently drift from what those options mean for a single
    system; sector runs take no forcing options, GEN.51), then -- if
    `--min-habitable` was given -- sets `HABITABLE_WORLD = True` on that
    many randomly-chosen configs.

    Args:
        args (argparse.Namespace): Parsed arguments, with
            `args.num_systems` already resolved to a concrete count --
            either the user's explicit `--num-systems`, or (see
            `generate_sector`) a per-sector value sampled from `--density`.

    Returns:
        list: A list of `num_systems` `SystemConfig` instances, one per
              system the sector will contain.

    Raises:
        SystemExit: Under `--strict`, if `--min-habitable` exceeds a
            `--density`-driven count (an explicit `--num-systems` is
            checked earlier by `validate_shared_generation_args`); without
            it, a warning and the sector gets `--min-habitable` systems
            (GEN.81).
    """
    return list(iter_sector_configs(args))


def iter_sector_configs(args):
    """
    `build_sector_configs`, lazily: yields the same `args.num_systems`
    configs one at a time, so a huge count (`--num-systems 1000000000`,
    or a huge `--density` draw) never builds them all up front --
    `generate_sector` stops pulling once the sector is full. Every config
    comes from the same `args`; the `--min-habitable` indices are chosen
    (and any conflict reported) before the first config is yielded.

    Raises:
        SystemExit: As `build_sector_configs`.
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
        forced = set(random.sample(range(count), k=args.min_habitable))

    for i in range(count):
        config = first if i == 0 else run_system.build_system_config(args)
        if i in forced:
            config.HABITABLE_WORLD = True
        yield config


def sector_star_count(sector):
    """How many stars `sector`'s systems hold (a binary counts two) --
    what every per-star phenomenon rate multiplies."""
    return sum(len(getattr(entry.star_system, "stars", None) or [None]) for entry in sector.entries)


def generate_sector_phenomena(sector, args, galactic_center_dist_ly=None, cloud_field=None):
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
        args (argparse.Namespace): Parsed arguments; only `args.markdown`
            is consulted (threaded into each phenomenon's own
            `SystemConfig`, matching every star system's own config).
        galactic_center_dist_ly (float, optional): As in `generate_sector`
            -- threaded into a standalone black hole's/neutron star's own
            Hill-sphere/galactic-orbit calculation, exactly like every star
            in the sector already gets.
        cloud_field (tuple, optional): For a galaxy-placed sector, the
            `nebula_field.clouds_reaching` arguments `(galaxy_seed, shape,
            center_pc, reach_pc)`: its molecular clouds then come from the
            galaxy's cloud field (GEN.47, into `sector.field_nebulae`)
            instead of this sector's own `"molecular-cloud"` roll.

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

    for kind, phenomenon_type, factory in SECTOR_PHENOMENON_KINDS:
        if kind == "molecular-cloud" and cloud_field is not None:
            sector.field_nebulae = nebula_field.clouds_reaching(*cloud_field)
            continue
        count = _sample_poisson_count(program_constants.phenomenon_rate_per_star(kind) * star_count)
        for _ in range(count):
            phenomenon_config = SystemConfig()
            phenomenon_config.MARKDOWN = args.markdown
            phenomenon = factory(phenomenon_config, galactic_center_dist_ly)
            if kind == "planetary-nebula":
                entry = _add_planetary_nebula(sector, args, phenomenon, galactic_center_dist_ly)
                if entry is not None:
                    new_entries.append(entry)
                continue
            try:
                new_entries.append(sector.add_phenomenon(phenomenon, phenomenon_type))
            except ValueError:
                # Only a massive type (black-hole/neutron-star) can raise
                # here (see SpaceSector._random_position) -- no room left
                # to place it without violating another massive object's
                # Hill sphere. Skip just this one draw rather than this
                # phenomenon type, or the whole sector.
                continue

    new_entries.extend(add_star_hosted_nebulae(sector, args))
    return new_entries


def _spectral_code(star):
    """A star's spectral code (`"G2V"`), the first word of `Star.type`."""
    return (getattr(star, "type", "") or "").split(" ", 1)[0]


def _add_planetary_nebula(sector, args, nebula, galactic_center_dist_ly=None):
    """
    Adds `nebula` (a planetary nebula) around a new system whose star is
    its hot central white dwarf (`PLANETARY_NEBULA_CENTRAL_STAR_TYPES`),
    placed like any other system. Returns the nebula's entry, or `None`
    when the sector has no room left for the system.
    """
    config = SystemConfig()
    config.MARKDOWN = args.markdown
    config.STAR_TYPE = random.choice(program_constants.PLANETARY_NEBULA_CENTRAL_STAR_TYPES)
    config.BINARY_SYSTEM = False
    system = StarSystem(system_config=config, galactic_center_dist_ly=galactic_center_dist_ly)
    try:
        system_entry = sector.add_system(system, system_config=config)
    except ValueError:
        return None
    log.debug(f"Sector {sector.name!r}: planetary nebula {nebula.name!r} around new central star "
              f"{config.STAR_TYPE} {system.name!r}")
    return sector.add_phenomenon(nebula, "nebula", position=system_entry.position)


def add_star_hosted_nebulae(sector, args):
    """
    Grows emission and reflection nebulae around `sector`'s own hot
    stars (GEN.10): each system whose primary is a main-sequence
    star matching a `NEBULA_HOST_RULES` row rolls that row's chance, and
    on a hit gets a nebula of one of the row's classes centered on it.

    Returns:
        list: The new nebulae's `SectorPhenomenonEntry` instances.
    """
    entries = []
    for system_entry in list(sector.entries):
        stars = getattr(system_entry.star_system, "stars", None) or []
        if not stars:
            continue
        code = _spectral_code(stars[0])
        match = re.fullmatch(r"([OBAFGKM])([0-9])V", code)
        if match is None:
            continue
        letter, subclass = match.group(1), int(match.group(2))
        for rule_letter, low, high, classes, chance in program_constants.NEBULA_HOST_RULES:
            if letter != rule_letter or not (low <= subclass <= high):
                continue
            if random.random() < chance:
                config = SystemConfig()
                config.MARKDOWN = args.markdown
                nebula_class = choose_weighted_class(classes, "Star-hosted nebula class")
                nebula = Nebula(config, nebula_class=nebula_class)
                log.debug(f"Sector {sector.name!r}: {code} star {system_entry.star_system.name!r} "
                          f"lights class {nebula_class} nebula {nebula.name!r}")
                entries.append(sector.add_phenomenon(nebula, "nebula", position=system_entry.position))
            break
    return entries


def flag_fast_stars(sector, galactic_center_dist_ly=None):
    """
    Marks some of `sector`'s systems as runaway or hypervelocity stars
    (`StarSystem.runaway_class` and `runaway_speed_kms`, schema v37): an
    ordinary generated system moving unusually fast, not a separate
    phenomenon. Each system rolls a hypervelocity chance first -- the
    "hypervelocity-star" rate per star scaled by
    `(HVS_REFERENCE_RADIUS_PC / r)^2`, since the central black hole ejects
    them -- then the runaway chance (about 1.5%).

    Args:
        sector (SpaceSector): The populated sector, changed in place.
        galactic_center_dist_ly (float, optional): The sector's distance
            from the galactic center; `None` uses
            `constants.GALACTIC_CENTER_DISTANCE_LY`.

    Returns:
        int: How many systems were flagged.
    """
    if galactic_center_dist_ly is None:
        galactic_center_dist_ly = constants.GALACTIC_CENTER_DISTANCE_LY
    radius_pc = max(ly_to_pc(galactic_center_dist_ly), 1.0)
    hvs_chance = min(1.0, program_constants.phenomenon_rate_per_star("hypervelocity-star")
                     * (program_constants.HVS_REFERENCE_RADIUS_PC / radius_pc) ** 2)
    runaway_chance = program_constants.phenomenon_rate_per_star("runaway-star")
    flagged = 0
    for entry in sector.entries:
        system = entry.star_system
        if random.random() < hvs_chance:
            system.runaway_class = "hypervelocity"
            system.runaway_speed_kms = log_uniform(*program_constants.HYPERVELOCITY_STAR_SPEED_RANGE_KMS)
        elif random.random() < runaway_chance:
            system.runaway_class = "runaway"
            system.runaway_speed_kms = log_uniform(*program_constants.RUNAWAY_STAR_SPEED_RANGE_KMS)
        else:
            continue
        flagged += 1
    return flagged


NUCLEUS_ADDRESS = (0, 0, 0)
"""tuple: The `(ring, layer, slot)` of the one sector per galaxy that rolls
for an active nucleus (see `add_galactic_nucleus`)."""


def add_galactic_nucleus(sector, args, galactic_center_dist_ly):
    """
    Adds the galaxy's central supermassive black hole to `sector` at the
    galactic center: an active `Quasar` with
    `program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE`, otherwise a quiescent
    supermassive `BlackHole` (like Sagittarius A*), so every galaxy has one.

    There is only ever one nucleus, and only at the origin.
    `generate_and_save_sector_at` calls this for exactly one sector per
    galaxy, `NUCLEUS_ADDRESS` (ring 0, layer 0, slot 0 -- every ring-0,
    layer-0 cell has the galactic axis as its inner edge and the plane
    through its middle, so each touches the origin, and picking one keeps
    it to a single roll).

    The sector's local +X axis points radially out from the galactic axis
    to the sector's own center, and a layer-0 center sits on the plane
    (`galaxyGeometry.sector_orientation`), so the center sits at
    `(-galactic_center_dist_ly, 0, 0)` in the sector's own frame;
    `store.insert_sector` converts that back to the galaxy origin.

    Args:
        sector (SpaceSector): The core sector, already populated.
        args (argparse.Namespace): Parsed arguments; only `args.markdown`
            is consulted.
        galactic_center_dist_ly (float): The sector center's distance from
            the galactic center, in light-years.

    Returns:
        SectorPhenomenonEntry: The quasar's or the black hole's entry.
    """
    nucleus_config = SystemConfig()
    nucleus_config.MARKDOWN = args.markdown
    position = (-galactic_center_dist_ly, 0.0, 0.0)
    if random.random() < program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE:
        return sector.add_phenomenon(Quasar(nucleus_config), "quasar", position=position)
    log.debug(f"Sector {sector.name!r}: galactic nucleus is quiescent (supermassive black hole, no quasar)")
    black_hole = BlackHole(nucleus_config, mass_class="supermassive")
    return sector.add_phenomenon(black_hole, "black-hole", position=position)


def _add_preplaced_systems(sector, args, fill, galactic_center_dist_ly):
    """Builds a full system around each of the sector's pre-placed bright
    stars (`fill.bright_rows`) and places it first, at its stored point.
    Each one's companion, planets and moons come from its own stored
    `seed`, without disturbing the run's random state."""
    for row in fill.bright_rows:
        cfg = run_system.build_system_config(args)
        cfg.POPULATION = row["population"]
        state = random.getstate()
        random.seed(row["seed"])
        try:
            system = StarSystem(system_config=cfg, galactic_center_dist_ly=galactic_center_dist_ly,
                                primary_star_params=brightStars.star_params(row))
        finally:
            random.setstate(state)
        entry = sector.add_preplaced_system(system, brightStars.local_position_ly(row, fill.center_pc),
                                            system_config=cfg)
        entry.bright_star_id = row["id"]


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
            `build_sector_configs`).
        galactic_center_dist_ly (float, optional): This sector's distance
            from the galactic center, in light-years, threaded down into
            every generated system's Hill-sphere calculation (see
            `Star.calculate_system_perimeter`'s docstring). `None` (the
            default) falls back to the fixed `GALACTIC_CENTER_DISTANCE_LY`
            constant -- the `sector` subcommand (no galaxy context) always
            calls this with the default, so it keeps producing "unplaced"
            sectors.
        cell (galaxyGeometry.SectorCell, optional): The sector's real
            cylindrical cell, in light-years, for a galaxy-placed sector --
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
    if fill is not None:
        _add_preplaced_systems(sector, args, fill, galactic_center_dist_ly)

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
                                  cloud_field=cloud_field)
        flag_fast_stars(sector, galactic_center_dist_ly=galactic_center_dist_ly)

    if density_driven and not sector.entries and not sector.phenomena:
        # Guaranteed non-empty: a sector's own Poisson draws
        # (system count here, each phenomenon type in generate_sector_phenomena)
        # are independent, so all of them landing on zero simultaneously is
        # a real, expected outcome at low means (e.g. ~13% at mean 2) --
        # but a sector with literally nothing in it isn't useful to anyone
        # visiting it, so force exactly one system rather than leave it
        # empty. Only applies when the count came from --density (an
        # explicit --num-systems, including 0, is a deliberate request
        # this never overrides).
        fallback_config = run_system.build_system_config(args)
        if fill is not None:
            fill.apply(fallback_config)
        fallback_system = StarSystem(system_config=fallback_config, galactic_center_dist_ly=galactic_center_dist_ly)
        sector.add_system(fallback_system, system_config=fallback_config)

    return sector_name, sector


_GIANT_YERKES = ("III", "II")


_SUPERGIANT_YERKES = ("IB", "IAB", "IA", "IA+", "0")


def _summary_star_label(star):
    """
    `(label, plural)` for one system's star in the sector summary (UX.34):
    white dwarfs, black holes and neutron stars by name, giants and
    supergiants under their letter, others by spectral letter alone.
    """
    yerkes = getattr(star, "yerkes_class", None)
    if yerkes in ("D", "VII"):
        return "white dwarf", "white dwarfs"
    if yerkes == "BH":
        return "black hole", "black holes"
    if yerkes == "NS":
        return "neutron star", "neutron stars"
    letter = (star.type or "?")[0]
    if yerkes in _GIANT_YERKES:
        return f"{letter}-type giant", f"{letter}-type giants"
    if yerkes in _SUPERGIANT_YERKES:
        return f"{letter}-type supergiant", f"{letter}-type supergiants"
    return f"{letter}-type", f"{letter}-type"


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
        types_str = ", ".join(f"{count} {label if count == 1 else plural}"
                              for (label, plural), count in sorted(system_types.items()))
    else:
        types_str = "none"
    lines.append(f"    Systems: {types_str}")

    if phenomenon_types:
        phenomena_str = ", ".join(
            f"{count} {run_phenomenon.TYPE_LABELS[ptype]}" for ptype, count in sorted(phenomenon_types.items())
        )
        lines.append(f"    Phenomena: {phenomena_str}")

    lines.append(
        f"    Star density: actual {actual_density:.2f}x local, expected {expected_density:.2f}x local"
    )

    return "\n".join(lines)


def run_sector(args):
    """
    Generates `args.num_sectors` sectors (1 by default), each one a
    fully populated `SpaceSector` (see `generate_sector` -- one
    `SystemConfig` per system via `build_sector_configs`, a full
    `StarSystem` from each, all placed in the sector via
    `SpaceSector.add_system`'s Hill-sphere-based placement) saved into
    the same database. No `galactic_center_dist_ly` is passed to
    `generate_sector` here -- the `sector` subcommand has no galaxy
    context (see the `galaxy` subcommand for that), so every sector it
    produces is "unplaced", regardless of `--num-sectors`.

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

    def saved(result, seconds, _weight):
        if queue.parallel:
            run_common.RUN_COUNTS["sectors"] += 1
            run_common.RUN_COUNTS["systems"] += result["systems"]
            run_common.RUN_COUNTS["phenomena"] += result["phenomena"]
        run_common._record_sector(args, result, seconds)
        phenomena_note = f", {result['phenomena']} phenomena" if result["phenomena"] else ""
        log.normal(
            f"Saved sector '{result['name']}' to the database (sector_id={result['sector_id']}, "
            f"{result['systems']} systems{phenomena_note}, "
            f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port})."
        )
        _log_summary(result["summary"])

    try:
        with run_common._work_queue(args, f"Sectors ({args.num_sectors} unplaced)") as queue:
            queue.expect(args.num_sectors)
            for index in range(args.num_sectors):
                queue.submit("sector", f"unplaced-{index}", _unplaced_sector_task, args, on_done=saved)
    finally:
        run_common._finish_stats(args)

    if args.num_sectors > 1:
        log.normal(f"Generated {args.num_sectors} sectors.")
    run_population.run_population_after(args)


def _unplaced_sector_task(args):
    """One `sector` subcommand sector, generated and saved -- a work queue
    task (see `_fill_sector_task`)."""
    _sector_name, sector = generate_sector(args)
    sector_id = store.save_sector(sector, config=store.mysql_config_from_args(args))
    run_common._count_sector(sector)
    return {
        "sector_id": sector_id, "name": sector.name, "systems": len(sector.entries),
        "stars": sector_star_count(sector), "density": run_common._sector_density(args),
        "phenomena": len(sector.phenomena), "summary": sector_generation_summary_lines(sector, args),
    }
