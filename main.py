import argparse
from stellarObjects import StarSystem

def parse_boolean_args():
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

    args = parser.parse_args()
    return args

def main():
    """
    Main function to generate and display a star system.
    """

    args = parse_boolean_args()

    if args.force_habitable_world and args.force_asteroid_belt:
        args.force_large_star = True

    system = StarSystem(force_hab=args.force_habitable_world, force_belt=args.force_asteroid_belt, force_large=args.force_large_star)
    print(system)

if __name__ == "__main__":
    main()