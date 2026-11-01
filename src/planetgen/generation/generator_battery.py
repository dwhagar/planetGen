# planetgen/generation/generator_battery.py

"""
The generator battery (OPS.57)
==============================

`generator_battery_digest()` runs a fixed battery of tiny canned cases, with no database and no network, and returns
one digest of the lot together with a digest for each case, so a tool can say that the generator still makes what it
made (OPS.15, OPS.14's mismatch check, TEST.77's tier A) and, when it does not, which case first differs
(docs/design/generation/generation-determinism.md section 5.4).

The cases, each from a fixed seed:

- `sector-a` and `sector-b`: a sector of three systems through `generate_sector`, under the unit seed of two sector
  addresses of one galaxy seed;
- `phenomenon`: one rogue planet from `generate_phenomenon`;
- `nebula-shape`: a nebula shape from `nebula_shape.draw_shape`, its scale and 50 `contains` answers;
- `bright-star-scatter`: the bright stars of a tiny outline of a toy galaxy (`bright_stars.scatter`);
- `population-stream`: ten stars of the bright-star population model (`star_population.sample_bright_stars`);
- `unit-seed`: `galaxy.seed.unit_seed` for a known address, which must never change.

A case is turned into a tree of plain values (`tree`: objects by their fields, sets sorted, the wall-clock fields of
`db.fingerprint_tree.CLOCK_KEYS` left out) and hashed with the encoder of GEN.304 (`util.fingerprint_encoding`) into
a `full`, a `structural` (floats written as `~`) and a `floats` digest. The overall digest is over the cases' full
digests in the order of `CASES`. A structural difference is a change of the code; a difference in the floats alone is
a change of the platform or its libraries (docs/design/generation/reproducible-galaxies.md).

Nothing in a case reads the clock, the environment or a stored file (only `kept_battery_digest` may look for a copy a test run kept, PERF.209), and every case runs on a stream of its own, so the
result does not depend on what ran before or on the hash seed.
"""

import concurrent.futures
import dataclasses
import math
import multiprocessing

from planetgen import tuning
from planetgen.db.fingerprint_tree import CLOCK_KEYS
from planetgen.util import draw
from planetgen.util import fingerprint_encoding as encoding

GALAXY_SEED = bytes(range(1, 17))
"""bytes: The galaxy seed of the sector cases and of the known answer."""

SECTOR_ADDRESSES = ((7, 1, 2), (12, -2, 0))
"""tuple: The `(shell, band, wedge)` of the two sector cases."""

SECTOR_SYSTEMS = 3
"""int: Systems in each sector case."""

NEBULA_SEED = 0x0EB0
"""int: The stream seed of the nebula case."""

NEBULA_POINTS = 50
"""int: Points the nebula case asks `contains` about."""

PHENOMENON_ADDRESS = (1, 0, 0)
"""tuple: The address whose unit seed the phenomenon case runs under."""

SCATTER_THRESHOLD_SOL = 2500.0
"""float: The scatter case's luminosity threshold, solar luminosities."""

SCATTER_EXPECTED_SYSTEMS = 2000.0
"""float: Systems per sector at density 1 in the scatter case's toy galaxy."""

SCATTER_EXTENTS = ((0, 0), (1, 1), (2, 2))
"""tuple: `(shell, highest band)` of the scatter case's outline: a ball three shells deep."""

SCATTER_SEED = 11
"""int: The scatter case's seed."""

POPULATION_SEED = 0x9091
"""int: The stream seed of the population case."""

POPULATION_STARS = 10
"""int: Stars the population case draws."""

KNOWN_UNIT_SEED = (bytes(range(16)), "sector", (3, 1, 0))
"""tuple: `(galaxy seed, kind, address)` of the known-answer case."""

CASES = ("sector-a", "sector-b", "phenomenon", "brown-dwarf", "nebula-shape", "bright-star-scatter", "population-stream", "unit-seed")
"""tuple: The cases, in the order the overall digest takes them."""


@dataclasses.dataclass(frozen=True)
class CaseDigest:
    """The three digests of one case (`util.fingerprint_encoding.Digests`) as hex text."""

    full: str
    structural: str
    floats: str


@dataclasses.dataclass(frozen=True)
class BatteryDigest:
    """
    The result of `generator_battery_digest`.

    Attributes:
        overall (str): The hex digest over every case's full digest, in the order of `CASES`.
        cases (dict): Case name to `CaseDigest`, in the order of `CASES`.
    """

    overall: str
    cases: dict

    def first_difference(self, other):
        """The first case, in the order of `CASES`, whose full digest differs from `other`'s, or `None`."""
        for name in CASES:
            if self.cases[name].full != other.cases[name].full:
                return name
        return None


def tree(value, seen=None):
    """
    `value` as a tree of the encoder's types: an object by its class name and its fields (callables left out), a set
    as a sorted list, a dict with its keys as text, a number from numpy as a Python one, and an object met again as
    its class name only. The clocks of `CLOCK_KEYS` are left out of every dict.
    """
    if seen is None:
        seen = set()
    if value is None or isinstance(value, (bool, str, bytes)):
        return value
    if hasattr(value, "tolist") and callable(value.tolist):
        return tree(value.tolist(), seen)
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, (list, tuple)):
        return [tree(item, seen) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((tree(item, seen) for item in value), key=lambda item: encoding.encode(item))
    if isinstance(value, dict):
        return {str(key): tree(item, seen) for key, item in value.items() if key not in CLOCK_KEYS}
    if id(value) in seen:
        return f"<{type(value).__name__}>"
    fields = getattr(value, "__dict__", None)
    if fields is None:
        slots = getattr(type(value), "__slots__", ())
        fields = {name: getattr(value, name) for name in slots if hasattr(value, name)}
    seen.add(id(value))
    return {"class": type(value).__name__,
            **{key: tree(item, seen) for key, item in fields.items() if not callable(item) and key not in CLOCK_KEYS}}


def _sector_case(address):
    """A sector of `SECTOR_SYSTEMS` systems under the unit seed of `address`, as its own dictionary."""
    from planetgen.cli import generate as generate_cli
    from planetgen.galaxy import seed as galaxySeed
    from planetgen.generation import run_sector

    parser, command_parsers = generate_cli.build_parser()
    args = parser.parse_args(["sector", "--num-systems", str(SECTOR_SYSTEMS), "--name", "Battery"])
    generate_cli.validate_shared_generation_args(args, command_parsers["sector"])
    generate_cli.validate_sector_args(args, command_parsers["sector"])
    with galaxySeed.seeded(GALAXY_SEED, "sector", address):
        _name, sector = run_sector.generate_sector(args)
    return tree(sector.to_dict())


def _phenomenon_case():
    """One rogue planet under a unit seed of its own."""
    from planetgen.galaxy import seed as galaxySeed
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.run_phenomenon import generate_phenomenon

    with galaxySeed.seeded(GALAXY_SEED, "phenomenon", PHENOMENON_ADDRESS):
        return tree(generate_phenomenon("rogue-planet", SystemConfig(), False, name="Battery Rogue"))


def _brown_dwarf_case():
    """One brown dwarf under a unit seed of its own: its record and the spectral type it is drawn (GEN.760)."""
    from planetgen.galaxy import seed as galaxySeed
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.rogue import RoguePlanet

    with galaxySeed.seeded(GALAXY_SEED, "brown-dwarf", PHENOMENON_ADDRESS):
        dwarf = RoguePlanet(SystemConfig(), name="Battery Dwarf", mass_bin="brown-dwarf")
        return {"record": tree(dwarf), "spectral_type": dwarf.spectral_type}


def _nebula_case():
    """A nebula shape, its dictionary and which of `NEBULA_POINTS` points it contains."""
    from planetgen.galaxy import nebula_shape

    stream = draw.Stream(NEBULA_SEED)
    shape = nebula_shape.draw_shape(stream)
    points = [[stream.uniform(-1.0, 1.0) for _ in range(3)] for _ in range(NEBULA_POINTS)]
    return {"shape": tree(shape.to_dict()), "scale": shape.scale, "inside": [bool(shape.contains(p)) for p in points]}


def _toy_shape():
    """A small galaxy read in relative units: density 1 at its reference radius between two arms."""
    from planetgen.galaxy import density, structure

    shape = density.build_galaxy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0, bulge_amplitude=2.0,
        arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4, structure=structure.REFERENCE_STRUCTURE)
    radius = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc
    angle = density._interarm_angle(radius, shape)
    position = (radius * math.cos(angle), radius * math.sin(angle), 0.0)
    return shape._replace(k_norm=shape.k_norm / density.relative_density(position, shape))


def _scatter_case():
    """The bright stars of the toy galaxy's tiny outline."""
    from planetgen.generation import bright_stars
    from planetgen.physics.units import ly_to_pc

    edge_pc = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
    rows = list(bright_stars.scatter(_toy_shape(), list(SCATTER_EXTENTS), edge_pc, SCATTER_EXPECTED_SYSTEMS,
                                     SCATTER_THRESHOLD_SOL, SCATTER_SEED))
    return [tree(row) for row in rows]


def _population_case():
    """Stars of the bright-star population model from a stream of its own."""
    from planetgen.generation import star_population

    stream = draw.Stream(POPULATION_SEED)
    return tree(star_population.sample_bright_stars(POPULATION_STARS, SCATTER_THRESHOLD_SOL, "old", stream))


def _unit_seed_case():
    """The unit seed of a known address, as an integer."""
    from planetgen.galaxy import seed as galaxySeed

    return galaxySeed.unit_seed(*KNOWN_UNIT_SEED)


_BUILDERS = {
    "sector-a": lambda: _sector_case(SECTOR_ADDRESSES[0]),
    "sector-b": lambda: _sector_case(SECTOR_ADDRESSES[1]),
    "phenomenon": _phenomenon_case,
    "brown-dwarf": _brown_dwarf_case,
    "nebula-shape": _nebula_case,
    "bright-star-scatter": _scatter_case,
    "population-stream": _population_case,
    "unit-seed": _unit_seed_case,
}


def case_value(name):
    """The value case `name` makes, before it is hashed (for a test or a report that shows what changed)."""
    return _BUILDERS[name]()


def generator_battery_digest(progress=None):
    """
    Runs the battery and returns its `BatteryDigest`. No database, no network, no clock; most of the time is the bright-star scatter
    case, which builds the sampling tables of the stellar populations.

    Args:
        progress (callable, optional): Called as `progress(done, total)` before the first case (`done` 0) and after each
            case, so the update can show a bar.
    """
    cases = {}
    if progress:
        progress(0, len(CASES))
    for done, name in enumerate(CASES, 1):
        found = encoding.digests(case_value(name))
        cases[name] = CaseDigest(found.full, found.structural, found.floats)
        if progress:
            progress(done, len(CASES))
    overall = encoding.digest({"battery": [[name, cases[name].full] for name in CASES]})
    return BatteryDigest(overall, cases)


_running = None
"""BatteryDigest or None: What `running_digest` has worked out in this process (or `start_running_digest` was handed)."""


def running_digest(progress=None):
    """
    The running code's `BatteryDigest`, worked out once a process: what a plan records with the galaxy and an update
    records in the version history (OPS.58), so the five seconds are paid once however many call it.

    Args:
        progress (callable, optional): `progress(done, total)` after each case, when the digest is worked out here.
    """
    global _running
    if _running is None:
        _running = kept_battery_digest(progress)
    return _running


def kept_battery_digest(progress=None):
    """`generator_battery_digest()`, or the copy a test run kept for the same code and data (`battery_digest_cache`, PERF.209)."""
    import functools

    from planetgen.generation import battery_digest_cache

    return battery_digest_cache.digest(functools.partial(generator_battery_digest, progress), BatteryDigest, CaseDigest)


def start_running_digest():
    """
    Starts working out `running_digest()` in a process of its own (PERF.176), for a caller that will need it later and
    has other work to do meanwhile: the battery uses no database and nothing of this process, so the digest is the
    one this process would make, five seconds that would otherwise be waited for on a free core.

    Returns:
        callable: `finish()`, which waits for the process and keeps its digest as this process's `running_digest()`. If
            the process could not be started or died, `running_digest()` works the digest out here when it is asked for.
    """
    if _running is not None:
        return lambda: None
    try:
        pool = concurrent.futures.ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn"))
        future = pool.submit(kept_battery_digest)
    except Exception:  # noqa: BLE001 -- no process to be had: the digest is worked out where it is asked for
        return lambda: None

    def finish():
        global _running
        try:
            _running = future.result()
        except Exception:  # noqa: BLE001 -- the digest is then worked out where it is asked for
            pass
        finally:
            pool.shutdown()

    return finish


def first_difference(stored, found):
    """
    The first case, in the order of `CASES`, whose full digest differs between two `{case: digest}` dicts (as `cases_of`
    gives and a record keeps), or `None` when none does. A case in only one of them is not compared, so a record from a
    battery with other cases says nothing of them.
    """
    for name in CASES:
        if name in stored and name in found and stored[name] != found[name]:
            return name
    return None


def cases_of(digest):
    """`{case: full digest}` of a `BatteryDigest`, the form a record keeps."""
    return {name: digest.cases[name].full for name in CASES}
