# tests/reproducible_sector_probe.py

"""
Run as a script by `test_reproducible_draws.py` (GEN.56): generates one
small sector on a fixed seed and prints the SHA-256 of a canonical dump
of everything in it, so the test can compare runs under different
`PYTHONHASHSEED` values and locales.
"""

import hashlib
import locale
import math
import sys


def canonical(value, seen=None):
    """`value` as text that depends only on its contents: attributes and
    dict items in sorted order, sets sorted, floats by `repr`, an object
    met again (a back reference) by its class name only."""
    if seen is None:
        seen = set()
    if value is None or isinstance(value, (bool, int, str)):
        return repr(value)
    if isinstance(value, float):
        return "nan" if math.isnan(value) else repr(value)
    if isinstance(value, (list, tuple)):
        return "[" + ",".join(canonical(item, seen) for item in value) + "]"
    if isinstance(value, (set, frozenset)):
        return "{" + ",".join(sorted(canonical(item, seen) for item in value)) + "}"
    if isinstance(value, dict):
        items = sorted((canonical(key, seen), canonical(item, seen)) for key, item in value.items())
        return "{" + ",".join(f"{key}:{item}" for key, item in items) + "}"
    if id(value) in seen:
        return f"<{type(value).__name__}>"
    fields = getattr(value, "__dict__", None)
    if fields is None:
        slots = getattr(type(value), "__slots__", ())
        fields = {name: getattr(value, name) for name in slots if hasattr(value, name)}
    if not fields:
        return f"<{type(value).__name__}>"
    seen.add(id(value))
    return type(value).__name__ + canonical({key: item for key, item in fields.items() if not callable(item)}, seen)


def main():
    locale.setlocale(locale.LC_ALL, "")
    from planetgen.cli import generate as generate_cli
    from planetgen.generation import run_sector
    from planetgen.util import draw

    parser, command_parsers = generate_cli.build_parser()
    args = parser.parse_args(["sector", "--num-systems", sys.argv[1] if len(sys.argv) > 1 else "4", "--name", "Probe"])
    command_parser = command_parsers["sector"]
    generate_cli.validate_shared_generation_args(args, command_parser)
    generate_cli.validate_sector_args(args, command_parser)
    with draw.bound(20261009):
        _name, sector = run_sector.generate_sector(args)
    print(hashlib.sha256(canonical(sector).encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
