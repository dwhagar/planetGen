import argparse
import json
import logging
import random
import secrets

from stellarObjects.config import SystemConfig
from stellarObjects import program_constants
from stellarObjects.systemData import StarSystem

# Suppress transformers warnings
logging.getLogger("transformers").setLevel(logging.ERROR)

# Tri-state options settable with `+name` (force True) / `-name` (force False).
# Each maps directly onto a `SystemConfig` attribute of the same name.
TRISTATE_OPTIONS = [
    ("habitable_world", "HABITABLE_WORLD", "a habitable world"),
    ("asteroid_belt", "ASTEROID_BELT", "an asteroid belt"),
    ("large_star", "LARGE_STAR", "a large, massive star"),
    ("moons", "MOONS", "moons on every planet"),
    ("max_planets", "MAX_PLANETS", "the maximum number of orbital objects"),
    ("intelligent_life", "INTELLIGENT_LIFE", "intelligent life on a planet"),
    ("binary_system", "BINARY_SYSTEM", "a binary (P-type) star system"),
    ("planets", "PLANETS", "at least one planet or asteroid belt"),
]


class TristateAction(argparse.Action):
    """
    Sets `namespace.dest` to True when invoked as `+name`, or False when
    invoked as `-name`. Leaving the option off the command line leaves the
    `default` (None) in place, meaning "let the generator decide".
    """

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, option_string.startswith('+'))


def process_args():
    """
    Parses command-line arguments for customizing the star system generation.

    Most of the boolean generation options use a `+name`/`-name` tri-state
    syntax rather than separate on/off flags: `+name` forces that feature to
    be present, `-name` forces it to be absent, and omitting it leaves the
    outcome up to the generator's normal random logic. These are (see
    `TRISTATE_OPTIONS`):
        - `+habitable_world` / `-habitable_world`
        - `+asteroid_belt` / `-asteroid_belt`
        - `+large_star` / `-large_star`
        - `+moons` / `-moons`
        - `+max_planets` / `-max_planets`
        - `+intelligent_life` / `-intelligent_life`
        - `+binary_system` / `-binary_system`
        - `+planets` / `-planets`

    The remaining, non-tri-state options:
    - `--system-file` / `-f`: Loads system generation options from a JSON
      file. Any tri-state or value option also given on the command line
      overrides the corresponding value from the file.
    - `--output` / `-o`: Specifies a file path to write the output to.
    - `--markdown` / `-m`: Formats the output in Markdown instead of the default
      wikitext.
    - `--star-type`: Allows specifying a precise star type to generate, such as
      'G2V'.
    - `--num-orbits`: Sets an exact number of orbital slots (planets and
      asteroid belts combined) to generate.
    - `--name`: Specifies a name for the star system, overriding the default
      random generation.
    - `--age`: Specifies the age of the star system, either "young" or "old".
    - `--flavor-chance-system`: Overrides the default system-level flavor text chance.
    - `--flavor-chance-planet`: Overrides the default planet-level flavor text chance.
    - `--max-planet-flavor`: Sets the maximum flavor text total for planets to 99.

    Per-orbit slot specifications (exact planet class, asteroid belt vs.
    planet, exact moon counts) are only settable via `--system-file`'s
    `"slots"` key -- see the module-level `SystemConfig.SLOTS` docstring.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments
                            as attributes. Each attribute corresponds to a specific
                            flag or option and its value.
    """
    additional_info = [
        "Additional Information:",
        "This tool is designed as a personal tool for the Molten Aether FFRP game, the output is designed to be",
        "simply cut and paste from the program output into the wiki, see https://wiki.moltenaether.com for wiki and",
        "game information.  Most generation options use a +name/-name syntax: +habitable_world forces a habitable",
        "world to be generated, -habitable_world forbids one, and omitting the option leaves it to chance. Use",
        "--system-file to load a full system specification (including exact per-orbit slot contents) from a JSON",
        "file. Forcing a habitable world and an asteroid belt will automatically force a large star to ensure there",
        "is room for both objects. The most common stars are small dwarf stars which make a smaller star system.",
        "Forcing a large star as well as the maximum number of planets will cause the generated system to be very",
        "large with an extremely high number of planets. Do not assume that just because it is generated here, it",
        "is accurate or possible, such large systems may require editing as some worlds may end up saying they are",
        "several thousand AU's from the central star."
    ]

    additional_info = " ".join(additional_info)

    parser = argparse.ArgumentParser(
        description="System Generation Options",
        epilog=additional_info,
        prefix_chars='-+')

    for name, _attr, description in TRISTATE_OPTIONS:
        parser.add_argument(f'-{name}', f'+{name}', dest=name, action=TristateAction,
                            nargs=0, default=None,
                            help=f"+{name} forces the system to have {description}; "
                                 f"-{name} forces the system to not have {description}.")

    # Load system options from a JSON file
    parser.add_argument('--system-file', '-f', type=str,
                        help="Load system generation options from a JSON file. Command-line options "
                             "override the values it sets.")

    # Output to a file
    parser.add_argument('--output', '-o', type=str, help="Output to a file.")

    # Output in Markdown format
    parser.add_argument('--markdown', '-m', action='store_true', help="Output in Markdown format.")

    # Star Type
    parser.add_argument('--star-type', type=str,
                        help="Force the generation of a specific star type (e.g., G2V).")

    # Number of Orbits
    parser.add_argument('--num-orbits', type=int,
                        help="Force an exact number of orbital slots (planets and asteroid belts "
                             "combined) to be generated.")

    # System Name
    parser.add_argument('--name', type=str,
                        help="Force the name of the star system.")

    # System Age
    parser.add_argument('--age', type=str, choices=['young', 'old'],
                        help="Specify the age of the star system (young or old).")

    # Override Flavor Chance System
    parser.add_argument('--flavor-chance-system', type=float,
                        help="Override the default FLAVOR_CHANCE_SYSTEM constant.")

    # Override Flavor Chance Planet
    parser.add_argument('--flavor-chance-planet', type=float,
                        help="Override the default FLAVOR_CHANCE_PLANET constant.")

    # Max Planet Flavor
    parser.add_argument('--max-planet-flavor', action='store_true',
                        help="Sets the maximum flavor text total for planets to 99.")

    args = parser.parse_args()

    # Argument validation logic
    if args.planets is False and (args.moons or args.max_planets or args.habitable_world):
        parser.error("-planets cannot be combined with +moons, +max_planets, or +habitable_world.")

    if args.star_type and args.large_star:
        parser.error("--star-type cannot be combined with +large_star.")

    if args.intelligent_life is not None and args.habitable_world is False:
        parser.error("+intelligent_life/-intelligent_life cannot be combined with -habitable_world.")

    if args.num_orbits is not None and args.num_orbits < 0:
        parser.error("--num-orbits must be zero or a positive integer.")

    if args.num_orbits is not None and args.planets is False:
        parser.error("--num-orbits cannot be combined with -planets.")

    if args.flavor_chance_system is not None and not (0.0 <= args.flavor_chance_system <= 1.0):
        parser.error("--flavor-chance-system must be a float between 0.0 and 1.0.")

    if args.flavor_chance_planet is not None and not (0.0 <= args.flavor_chance_planet <= 1.0):
        parser.error("--flavor-chance-planet must be a float between 0.0 and 1.0.")

    return args


def load_system_file(path):
    """
    Loads a system generation specification from a JSON file.

    The JSON file may contain any of the following keys, each corresponding
    to a `SystemConfig` attribute of the same name (see that class's
    docstrings for details): `star_type`, `name`, `age`, `num_orbits`,
    `slots`, `habitable_world`, `asteroid_belt`, `large_star`, `moons`,
    `max_planets`, `intelligent_life`, `binary_system`, `planets`,
    `markdown`, `flavor_chance_system`, `flavor_chance_planet`,
    `max_planet_flavor`, `output`.

    `slots`, if present, is a list whose entries are either `null` or an
    object with a required `"type"` ("planet" or "asteroid_belt") and, for
    planets, optional `"planet_class"` (e.g. "M") and `"moons"` (an exact
    moon count) keys, e.g.:

        {
          "star_type": "G2V",
          "num_orbits": 5,
          "habitable_world": true,
          "slots": [
            {"type": "planet", "planet_class": "M", "moons": 1},
            {"type": "asteroid_belt"},
            null,
            {"type": "planet", "planet_class": "J", "moons": 4},
            null
          ]
        }

    Args:
        path (str): Path to the JSON system specification file.

    Returns:
        dict: The parsed JSON content.
    """
    with open(path, 'r') as f:
        return json.load(f)


def apply_system_file(system_config, data):
    """
    Applies a system specification (as loaded by `load_system_file`) onto a
    `SystemConfig` instance.

    Args:
        system_config (SystemConfig): The config object to update in place.
        data (dict): The parsed JSON system specification.
    """
    simple_keys = [
        "star_type", "name", "age", "num_orbits", "slots",
        "habitable_world", "asteroid_belt", "large_star", "moons",
        "max_planets", "intelligent_life", "binary_system", "planets",
        "markdown",
    ]
    key_to_attr = {key: key.upper() for key in simple_keys}

    for key, attr in key_to_attr.items():
        if key in data:
            setattr(system_config, attr, data[key])

    if "flavor_chance_system" in data and data["flavor_chance_system"] is not None:
        program_constants.FLAVOR_CHANCE_SYSTEM = data["flavor_chance_system"]

    if "flavor_chance_planet" in data and data["flavor_chance_planet"] is not None:
        program_constants.FLAVOR_CHANCE_PLANET = data["flavor_chance_planet"]

    if data.get("max_planet_flavor"):
        program_constants.MAX_FLAVOR_TOTAL = 99
        program_constants.FLAVOR_CHANCE_PLANET = 1


def build_system_config(args):
    """
    Builds a fully-configured `SystemConfig` from parsed command-line
    arguments, applying the same precedence `main()` always has: a
    `--system-file` JSON spec (if given) sets the baseline, then any
    explicitly-given tri-state or value option on the command line overrides
    it, then the cross-option normalization rules (intelligent life implies
    a habitable world; a habitable world plus an asteroid belt implies a
    large star) are applied last.

    Shared with `sectorGen.py`, which builds one `SystemConfig` per system
    in a sector this same way, so the two scripts can never silently drift
    apart on what a given set of options actually means.

    Args:
        args (argparse.Namespace): Parsed arguments from `process_args()`
                                   (or a namespace shaped the same way).

    Returns:
        tuple: (SystemConfig, output_path) -- output_path is `args.output`,
              or the `--system-file`'s own `"output"` key if `args.output`
              wasn't given.

    Raises:
        SystemExit: If a habitable world and an asteroid belt are both
                   forced while a large star is explicitly forbidden (see
                   TRISTATE_OPTIONS's `large_star`) -- that combination has
                   no room to be satisfied.
    """
    system_config = SystemConfig()
    output_path = args.output

    if args.system_file:
        file_data = load_system_file(args.system_file)
        apply_system_file(system_config, file_data)
        if output_path is None and file_data.get("output"):
            output_path = file_data["output"]

    # Command-line tri-state options override anything set by --system-file.
    for name, attr, _description in TRISTATE_OPTIONS:
        value = getattr(args, name)
        if value is not None:
            setattr(system_config, attr, value)

    # Command-line value options override anything set by --system-file.
    if args.markdown:
        system_config.MARKDOWN = True
    if args.star_type is not None:
        system_config.STAR_TYPE = args.star_type
    if args.num_orbits is not None:
        system_config.NUM_ORBITS = args.num_orbits
    if args.name is not None:
        system_config.NAME = args.name
    if args.age is not None:
        system_config.AGE = args.age

    if args.flavor_chance_system is not None:
        program_constants.FLAVOR_CHANCE_SYSTEM = args.flavor_chance_system

    if args.flavor_chance_planet is not None:
        program_constants.FLAVOR_CHANCE_PLANET = args.flavor_chance_planet

    if args.max_planet_flavor:
        program_constants.MAX_FLAVOR_TOTAL = 99
        program_constants.FLAVOR_CHANCE_PLANET = 1

    # Intelligent life requires (and forbids the absence of) a habitable world.
    if system_config.INTELLIGENT_LIFE is not None:
        system_config.HABITABLE_WORLD = True

    if system_config.HABITABLE_WORLD is True and system_config.ASTEROID_BELT is True:
        if system_config.LARGE_STAR is False:
            raise SystemExit("Error: forcing both a habitable world and an asteroid belt requires a large "
                              "star; -large_star cannot be combined with +habitable_world and +asteroid_belt.")
        system_config.LARGE_STAR = True

    return system_config, output_path


def main():
    """
    The main entry point for the star system generation script.

    This function orchestrates the entire process of generating a star system.
    It parses command-line arguments, builds a `SystemConfig` from them via
    `build_system_config` (see that function for the JSON-file/CLI precedence
    and cross-option normalization rules it applies), and instantiates the
    `StarSystem` class, which uses that config to guide the generation of the
    star and any accompanying planets.

    Finally, the script checks if an output file was specified. If so, it writes
    the generated system's string representation to that file. Otherwise, it
    prints the output directly to the console.
    """
    # Seed the random number generator with a cryptographically secure seed
    random.seed(secrets.randbits(128))

    args = process_args()
    system_config, output_path = build_system_config(args)
    system = StarSystem(system_config=system_config)

    if output_path:
        with open(output_path, 'w') as f:
            f.write(str(system))
    else:
        print(system)

if __name__ == "__main__":
    main()
