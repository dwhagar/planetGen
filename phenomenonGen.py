"""
phenomenonGen.py
================

Generates a single exotic stellar phenomenon -- a black hole, neutron
star, nebula, supernova remnant, rogue planet, interstellar comet, or
standalone asteroid field -- kept deliberately separate from
`systemGen.py`'s normal system generation.

Per this feature's design, these phenomena are NOT part of normal system
generation odds: `StarSystem._generate_planets`'s per-slot rolls never
produce one, and `systemGen.py`/`sectorGen.py` never reference this
module. They are reachable only through this script's own, rarer,
on-demand `--type` choice (uniformly random among all seven when omitted).

A black hole or neutron star may optionally anchor a full `StarSystem`
(`--anchor-system`) -- e.g. a pulsar with a fallback-disk planet (real
examples exist, PSR B1257+12) -- reusing all of `StarSystem`'s existing
orbit-placement/rendering logic via `compactRemnant.py`'s `Star` subclass
design (see that module's docstring). The other five phenomena are always
standalone; `--anchor-system` doesn't apply to them.

A nebula or standalone asteroid field may optionally be placed in the
galaxy (`--sector-id`), near an already galaxy-placed sector -- see
`stellarObjects._db.compute_phenomenon_placement`/schema.sql's "v18"
header note for why this is a galaxy-frame sphere rather than a
sector-relative offset. Every other phenomenon type stays unplaced;
`--sector-id` doesn't apply to them.
"""

import argparse
import logging
import os
import random
import secrets
import sys

# stellarObjects lives at src/stellarObjects (src layout) -- add src/ to the
# import path so this keeps working without requiring `pip install .`
# first, matching how html/'s CGI scripts fall back to a no-install layout.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from stellarObjects.asteroidFieldData import AsteroidField
from stellarObjects.compactRemnant import BlackHole, NeutronStar
from stellarObjects.config import SystemConfig
from stellarObjects import _db, program_constants
from stellarObjects._version import VersionAction, version_banner
from stellarObjects.nebulaData import Nebula
from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet
from stellarObjects.supernovaRemnantData import SupernovaRemnant
from stellarObjects.systemData import StarSystem

# Suppress transformers warnings
logging.getLogger("transformers").setLevel(logging.ERROR)

TYPE_LABELS = {
    "black-hole": "black hole",
    "neutron-star": "neutron star",
    "nebula": "nebula",
    "supernova-remnant": "supernova remnant",
    "rogue-planet": "rogue planet",
    "comet": "interstellar comet",
    "asteroid-field": "asteroid field",
}
"""dict: `--type` value -> human-readable label, used in this script's own
status output (not part of the generated phenomenon's own text)."""


def process_args():
    """
    Parses command-line arguments for generating one exotic phenomenon.

    Unlike `systemGen.py`'s tri-state `+name`/`-name` generation options
    (which bias/force outcomes within a normally-generated system), this
    script has no such options: each `--type` always generates, since
    reaching this generation mode at all is already the "forced" choice
    (see module docstring).

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    additional_info = [
        "Additional Information:",
        "Generates a single exotic stellar phenomenon, kept separate from systemGen.py's normal system",
        "generation. Omitting --type picks uniformly at random among all seven phenomena. --anchor-system",
        "(black-hole/neutron-star only) builds a full star system around the compact remnant instead of",
        "describing it standalone -- disk-physics-driven planet generation naturally tends toward zero",
        "planets around a dark remnant, so pass --num-orbits (or edit the generated SystemConfig) to force",
        "orbiting bodies if you want them. --sector-id (nebula/asteroid-field only) places the generated",
        "phenomenon in the galaxy near an already galaxy-placed sector, instead of leaving it unplaced.",
    ]
    additional_info = " ".join(additional_info)

    parser = argparse.ArgumentParser(
        description="Exotic Stellar Phenomenon Generation Options",
        epilog=additional_info)

    parser.add_argument('--version', action=VersionAction, banner=version_banner('phenomenonGen.py'))

    parser.add_argument('--type', type=str, choices=list(program_constants.PHENOMENON_TYPE_CHOICES),
                         help="The kind of phenomenon to generate. Omit to pick uniformly at random.")

    parser.add_argument('--anchor-system', action='store_true',
                         help="Only valid with --type black-hole or --type neutron-star: build a full star "
                              "system (with orbiting planets, if any) anchored by the compact remnant, "
                              "instead of describing it standalone.")

    parser.add_argument('--num-orbits', type=int,
                         help="With --anchor-system, force an exact number of orbital slots around the "
                              "compact remnant (see this script's epilog).")

    parser.add_argument('--sector-id', type=int,
                         help="Only valid with --type nebula or --type asteroid-field: place the generated "
                              "phenomenon in the galaxy near this already-generated, already galaxy-placed "
                              "sector (see stellarObjects._db.compute_phenomenon_placement). Omit to generate "
                              "it unplaced, as before.")

    # Output to a file
    parser.add_argument('--output', '-o', type=str, help="Output to a file.")

    # Database persistence
    _db.add_mysql_connection_args(parser)

    # Output in Markdown format
    parser.add_argument('--markdown', '-m', action='store_true', help="Output in Markdown format.")

    # Name
    parser.add_argument('--name', type=str, help="Force the name of the generated phenomenon.")

    args = parser.parse_args()

    if args.num_orbits is not None and not args.anchor_system:
        parser.error("--num-orbits requires --anchor-system.")
    if args.num_orbits is not None and args.num_orbits < 0:
        parser.error("--num-orbits must be zero or a positive integer.")
    if args.anchor_system and args.type not in (None, "black-hole", "neutron-star"):
        parser.error("--anchor-system is only valid with --type black-hole or --type neutron-star.")
    if args.sector_id is not None and args.type not in (None, "nebula", "asteroid-field"):
        parser.error("--sector-id is only valid with --type nebula or --type asteroid-field.")

    return args


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
        `RoguePlanet`, `InterstellarComet`, or `AsteroidField`.

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

    raise ValueError(f"Unknown phenomenon type: {phenomenon_type!r}")


def main():
    """
    The main entry point for exotic phenomenon generation.

    Parses arguments, picks a phenomenon type (`--type`, or uniformly at
    random if omitted), generates it via `generate_phenomenon`, writes its
    string representation to `--output` or stdout, then persists it to the
    database via `stellarObjects._db.save_phenomenon` (mirrors
    `systemGen.py main`'s own unconditional save).
    """
    # Seed the random number generator with a cryptographically secure seed
    random.seed(secrets.randbits(128))

    args = process_args()
    phenomenon_type = args.type or random.choice(program_constants.PHENOMENON_TYPE_CHOICES)

    system_config = SystemConfig()
    system_config.MARKDOWN = args.markdown
    if args.num_orbits is not None:
        system_config.NUM_ORBITS = args.num_orbits

    phenomenon = generate_phenomenon(phenomenon_type, system_config, args.anchor_system, name=args.name)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(str(phenomenon))
    else:
        print(phenomenon)

    mysql_config = _db.mysql_config_from_args(args)
    phenomenon_id = _db.save_phenomenon(phenomenon, system_config, phenomenon_type, config=mysql_config,
                                         sector_id=args.sector_id)
    print(f"Saved {TYPE_LABELS[phenomenon_type]} to the database (id={phenomenon_id}, "
          f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")


if __name__ == "__main__":
    main()
