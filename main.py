import argparse
from stellarObjects import StarSystem

def process_args():
    """Parses command-line arguments for boolean flags related to system generation.

    Returns:
        argparse.Namespace: An object containing parsed argument values as boolean attributes.
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

    args = parser.parse_args()
    return args

def main():
    """
    Main function to generate and display a star system.
    """

    args = process_args()

    if args.force_habitable_world and args.force_asteroid_belt and not args.absurd:
        args.force_large_star = True
    elif args.absurd:
        args.force_max_planets = True
        args.force_moons = True

    system = StarSystem(force_hab=args.force_habitable_world,
                        force_belt=args.force_asteroid_belt,
                        force_large=args.force_large_star,
                        force_moons=args.force_moons,
                        force_planets=args.force_max_planets,
                        absurd=args.absurd)
    print(system)

if __name__ == "__main__":
    main()