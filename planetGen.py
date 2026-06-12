import argparse
import logging
from stellarObjects import StarSystem, config

# Suppress transformers warnings
logging.getLogger("transformers").setLevel(logging.ERROR)

def process_args():
    """
    Parses command-line arguments for customizing the star system generation.

    This function sets up an argument parser to handle various flags and options
    that allow the user to control the star system creation process. It includes
    options to force the inclusion of specific features like a habitable world or
    an asteroid belt, control the size and composition of the system, and specify
    output formats.

    The parser includes the following arguments:
    - `--force-habitable-world` / `-fhw`: Ensures at least one habitable planet
      is generated.
    - `--force-asteroid-belt` / `-fab`: Ensures at least one asteroid belt is
      generated.
    - `--force-large-star` / `-fls`: Forces the generation of a larger, more
      massive star.
    - `--force-moons` / `-fm`: Increases the likelihood of planets having moons.
    - `--force-max-planets` / `-fmp`: Generates the maximum number of planets
      the system can support.
    - `--absurd`: Creates an extremely large and dense system for creative
      purposes.
    - `--output` / `-o`: Specifies a file path to write the output to.
    - `--markdown` / `-m`: Formats the output in Markdown instead of the default
      wikitext.
    - `--no-planets` / `-np`: Skips the planet generation process entirely,
      creating only a star.
    - `--star-type`: Allows specifying a precise star type to generate, such as
      'G2V'.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments
                            as attributes. Each attribute corresponds to a specific
                            flag or option and its value.
    """
    additional_info = [
        "Additional Information:",
        "This tool is designed as a personal tool for the Molten Aether FFRP game, the output is designed to be",
        "simply cut and paste from the program output into the wiki, see https://wiki.moltenaether.com for wiki and",
        "game information.  Using commands to force a habitable world and an asteroid belt will automatically force",
        "a large star to ensure there is room for both objects. The most common stars are small dwarf stars which",
        "make a smaller star system.  Forcing a large star as well as maximized planets will cause the generated",
        "system to be very large with an extremely high and sometimes absurd number of planets.  Do not assume that",
        "just because it is generated here, it is accurate or possible, such large systems may require editing as",
        "some worlds may end up saying they are several thousand AU's from the central star."
    ]

    additional_info = " ".join(additional_info)

    parser = argparse.ArgumentParser(
        description="System Generation Options",
        epilog=additional_info)

    # Force Habitable World
    parser.add_argument('--force-habitable-world', '-fhw', action='store_true',
                        help="Force the generation of a habitable world in the system.")

    # Force Asteroid Belt
    parser.add_argument('--force-asteroid-belt', '-fab', action='store_true',
                        help="Force the generation of an asteroid belt in the system.")

    # Force Large Star
    parser.add_argument('--force-large-star', '-fls', action='store_true',
                        help="Force the generation of a large star.")

    # Force Lots of Moons
    parser.add_argument('--force-moons', '-fm', action='store_true',
                        help="Force the generation of lots of moons in a system.")

    # Force Maximum Planets
    parser.add_argument('--force-max-planets', '-fmp', action='store_true',
                        help="Force the system to maximized the number of planets.")

    # Make the biggest and most absurd system possible.
    parser.add_argument('--absurd', action='store_true',
    help="Force the system to generate the largest star possible w/ max planets and moons.")

    # Output to a file
    parser.add_argument('--output', '-o', type=str, help="Output to a file.")

    # Output in Markdown format
    parser.add_argument('--markdown', '-m', action='store_true', help="Output in Markdown format.")

    # No Planets
    parser.add_argument('--no-planets', '-np', action='store_true',
                        help="Do not generate any planets, only the star.")

    # Star Type
    parser.add_argument('--star-type', type=str,
                        help="Force the generation of a specific star type (e.g., G2V).")


    args = parser.parse_args()

    if args.no_planets and (args.force_moons or args.force_max_planets or args.absurd or args.force_habitable_world or args.force_asteroid_belt):
        parser.error("--no-planets cannot be combined with --force-moons, --force-max-planets, --absurd, --force-habitable-world, or --force-asteroid-belt.")

    if args.star_type and args.force_large_star:
        parser.error("--star-type cannot be combined with --force-large-star.")

    return args

def main():
    """
    The main entry point for the star system generation script.

    This function orchestrates the entire process of generating a star system. It
    begins by parsing command-line arguments to retrieve user-specified options.
    It then processes these options, adjusting generation parameters accordingly.
    For instance, it can force the creation of a large star if both a habitable
    world and an asteroid belt are requested, or it can create an "absurd" system
    with maximum planets and moons.

    The function then instantiates the `StarSystem` class, which uses the
    centralized `config` object to guide the generation of the star and any
    accompanying planets.

    Finally, the script checks if an output file was specified. If so, it writes
    the generated system's string representation to that file. Otherwise, it

    prints the output directly to the console.
    """
    args = process_args()

    config.MARKDOWN = args.markdown
    config.FORCE_HABITABLE_WORLD = args.force_habitable_world
    config.FORCE_ASTEROID_BELT = args.force_asteroid_belt
    config.FORCE_LARGE_STAR = args.force_large_star
    config.FORCE_MOONS = args.force_moons
    config.FORCE_MAX_PLANETS = args.force_max_planets
    config.ABSURD = args.absurd
    config.NO_PLANETS = args.no_planets
    config.STAR_TYPE = args.star_type

    if config.FORCE_HABITABLE_WORLD and config.FORCE_ASTEROID_BELT and not config.ABSURD:
        config.FORCE_LARGE_STAR = True
    elif config.ABSURD:
        config.FORCE_MAX_PLANETS = True
        config.FORCE_MOONS = True

    system = StarSystem()

    if args.output:
        with open(args.output, 'w') as f:
            f.write(str(system))
    else:
        print(system)

if __name__ == "__main__":
    main()