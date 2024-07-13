import argparse
from stellarObjects import StarSystem

def parse_args():
    """Parses command-line arguments for boolean flags related to system generation.

    Returns:
        argparse.Namespace: An object containing parsed argument values as boolean attributes.
    """
    parser = argparse.ArgumentParser(description="System Generation Options")

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
                        help="Force the system  to maximized the number of planets.")

    args = parser.parse_args()
    return args

def main():
    """
    Main function to generate and display a star system.
    """

    args = parse_args()

    if args.force_habitable_world and args.force_asteroid_belt:
        args.force_large_star = True

    system = StarSystem(force_hab=args.force_habitable_world,
                        force_belt=args.force_asteroid_belt,
                        force_large=args.force_large_star,
                        force_moons=args.force_moons,
                        force_planets=args.force_max_planets)
    print(system)

if __name__ == "__main__":
    main()