# planetgen/generation/run_phenomenon.py

"""
Run Phenomenon
==============

The `phenomenon` command: one exotic stellar phenomenon (a black hole,
neutron star, nebula, supernova remnant, rogue planet, interstellar comet
or asteroid field), on demand. The sector command reuses
`generate_phenomenon` to seed every sector with its own population of
them.
"""

import random

from planetgen.db import store
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.system import StarSystem
from planetgen.generation import run_common


TYPE_LABELS = {
    "black-hole": "black hole",
    "neutron-star": "neutron star",
    "nebula": "nebula",
    "supernova-remnant": "supernova remnant",
    "rogue-planet": "rogue planet",
    "comet": "interstellar comet",
    "asteroid-field": "asteroid field",
    "quasar": "quasar",
}
"""dict: `--type` value -> human-readable label, used in the `phenomenon`
subcommand's own status output (not part of the generated phenomenon's
own text)."""


def generate_phenomenon(phenomenon_type, system_config, anchor_system, name=None):
    """
    Builds one generated phenomenon object for `phenomenon_type`.

    Args:
        phenomenon_type (str): One of `program_constants.PHENOMENON_TYPE_CHOICES`.
        system_config (SystemConfig): The shared config to generate with.
        anchor_system (bool): Only consulted for `"black-hole"`/
            `"neutron-star"` -- if True, wraps the generated compact
            remnant in a full `StarSystem` (see
            `StarSystem.__init__`'s `compact_remnant` parameter) instead
            of returning it standalone.
        name (str, optional): An explicit name for the phenomenon (or, for
            an anchored system, for the compact remnant itself).

    Returns:
        The generated phenomenon: a `BlackHole`, `NeutronStar`, `StarSystem`
        (only when `anchor_system` is True), `Nebula`, `SupernovaRemnant`,
        `RoguePlanet`, `InterstellarComet`, `AsteroidField`, or `Quasar`.

    Raises:
        ValueError: If `phenomenon_type` isn't one of the recognized choices.
    """
    if phenomenon_type == "black-hole":
        remnant = BlackHole(system_config, name=name)
        return StarSystem(system_config=system_config, compact_remnant=remnant) if anchor_system else remnant
    if phenomenon_type == "neutron-star":
        remnant = NeutronStar(system_config, name=name)
        return StarSystem(system_config=system_config, compact_remnant=remnant) if anchor_system else remnant
    if phenomenon_type == "nebula":
        return Nebula(system_config, name=name)
    if phenomenon_type == "supernova-remnant":
        return SupernovaRemnant(system_config, name=name)
    if phenomenon_type == "rogue-planet":
        return RoguePlanet(system_config, name=name)
    if phenomenon_type == "comet":
        return InterstellarComet(system_config, name=name)
    if phenomenon_type == "asteroid-field":
        return AsteroidField(system_config, name=name)
    if phenomenon_type == "quasar":
        return Quasar(system_config, name=name)

    raise ValueError(f"Unknown phenomenon type: {phenomenon_type!r}")


def run_phenomenon(args):
    """
    Generates one exotic phenomenon and saves it to the database.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "phenomenon"`).
    """
    phenomenon_type = args.type or random.choice(program_constants.RANDOM_PHENOMENON_TYPE_CHOICES)

    if args.sector_id is not None:
        # Checked before generating anything, as a clean one-line error.
        conn = store.get_connection(store.mysql_config_from_args(args))
        try:
            store.get_sector_galaxy_position(conn, args.sector_id)
        except ValueError as exc:
            log.error(f"Error: --sector-id {args.sector_id}: {exc}")
            raise SystemExit(1) from exc
        finally:
            conn.close()

    system_config = SystemConfig()
    system_config.MARKDOWN = args.markdown
    if args.num_orbits is not None:
        system_config.NUM_ORBITS = args.num_orbits

    phenomenon = generate_phenomenon(phenomenon_type, system_config, args.anchor_system, name=args.name)

    mysql_config = store.mysql_config_from_args(args)
    phenomenon_id = store.save_phenomenon(phenomenon, system_config, phenomenon_type, config=mysql_config,
                                         sector_id=args.sector_id)
    run_common.RUN_COUNTS["phenomena"] += 1
    log.normal(f"Saved {TYPE_LABELS[phenomenon_type]} to the database (id={phenomenon_id}, "
               f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")
