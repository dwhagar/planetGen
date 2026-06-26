import argparse
import logging
from stellarObjects import constants
from stellarObjects.systemData import StarSystem
from stellarObjects.config import SystemConfig

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
    - `--name`: Specifies a name for the star system, overriding the default
      random generation.
    - `--age`: Specifies the age of the star system, either "young" or "old".
    - `--force-intelligent-life` / `-fil`: Ensures at least one planet with
      intelligent life is generated. Implies --force-habitable-world.
    - `--no-intelligent-life` / `-nil`: Ensures no planet with intelligent life
      is generated. Implies --force-habitable-world.
    - `--no-habitable-world` / `-nhw`: Ensures no habitable world is generated
      in the system.
    - `--flavor-chance-system`: Overrides the default system-level flavor text chance.
    - `--flavor-chance-planet`: Overrides the default planet-level flavor text chance.
    - `--max-planet-flavor`: Sets the maximum flavor text total for planets to 99.

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

    # System Name
    parser.add_argument('--name', type=str,
                        help="Force the name of the star system.")

    # System Age
    parser.add_argument('--age', type=str, choices=['young', 'old'],
                        help="Specify the age of the star system (young or old).")

    # Force Intelligent Life
    parser.add_argument('--force-intelligent-life', '-fil', action='store_true',
                        help="Force the generation of intelligent life on a planet. Implies --force-habitable-world.")

    # No Intelligent Life
    parser.add_argument('--no-intelligent-life', '-nil', action='store_true',
                        help="Ensure no intelligent life is generated on any planet. Implies --force-habitable-world.")

    # No Habitable World
    parser.add_argument('--no-habitable-world', '-nhw', action='store_true',
                        help="Ensures no habitable world is generated in the system.")

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
    if args.no_planets and (args.force_moons or args.force_max_planets or args.absurd or args.force_habitable_world):
        parser.error("--no-planets cannot be combined with --force-moons, --force-max-planets, --absurd, or --force-habitable-world.")

    if args.star_type and args.force_large_star:
        parser.error("--star-type cannot be combined with --force-large-star.")

    if args.force_intelligent_life and args.no_intelligent_life:
        parser.error("--force-intelligent-life and --no-intelligent-life cannot be used together.")

    if args.force_habitable_world and args.no_habitable_world:
        parser.error("--force-habitable-world and --no-habitable-world cannot be used together.")

    if args.flavor_chance_system is not None and not (0.0 <= args.flavor_chance_system <= 1.0):
        parser.error("--flavor-chance-system must be a float between 0.0 and 1.0.")

    if args.flavor_chance_planet is not None and not (0.0 <= args.flavor_chance_planet <= 1.0):
        parser.error("--flavor-chance-planet must be a float between 0.0 and 1.0.")

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

    # Override global constants at module level before initializing system configurations.
    if args.flavor_chance_system is not None:
        constants.FLAVOR_CHANCE_SYSTEM = args.flavor_chance_system

    if args.flavor_chance_planet is not None:
        constants.FLAVOR_CHANCE_PLANET = args.flavor_chance_planet

    if args.max_planet_flavor:
        constants.MAX_FLAVOR_TOTAL = 99
        constants.FLAVOR_CHANCE_PLANET = 1

    system_config = SystemConfig() # Create an instance of SystemConfig

    # Assign parsed arguments to the SystemConfig object
    system_config.MARKDOWN = args.markdown
    system_config.FORCE_HABITABLE_WORLD = args.force_habitable_world
    system_config.FORCE_ASTEROID_BELT = args.force_asteroid_belt
    system_config.FORCE_LARGE_STAR = args.force_large_star
    system_config.FORCE_MOONS = args.force_moons
    system_config.FORCE_MAX_PLANETS = args.force_max_planets
    system_config.ABSURD = args.absurd
    system_config.NO_PLANETS = args.no_planets
    system_config.STAR_TYPE = args.star_type
    system_config.NAME = args.name
    system_config.AGE = args.age
    system_config.FORCE_INT = args.force_intelligent_life
    system_config.NO_INT = args.no_intelligent_life
    system_config.NO_HABITABLE_WORLD = args.no_habitable_world

    # Apply conditional logic based on SystemConfig flags
    if system_config.FORCE_INT or system_config.NO_INT:
        system_config.FORCE_HABITABLE_WORLD = True

    if system_config.FORCE_HABITABLE_WORLD and system_config.FORCE_ASTEROID_BELT and not system_config.ABSURD:
        system_config.FORCE_LARGE_STAR = True
    elif system_config.ABSURD:
        system_config.FORCE_MAX_PLANETS = True
        system_config.FORCE_MOONS = True

    system = StarSystem(system_config=system_config) # Pass the SystemConfig instance

    if args.output:
        with open(args.output, 'w') as f:
            f.write(str(system))
    else:
        print(system)

if __name__ == "__main__":
    main()