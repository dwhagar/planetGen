#!/usr/bin/env python3
# planetgen.cli.maintenance

"""
Turns maintenance mode on or off (UX.112).

    python3 -m planetgen.cli.maintenance on [--eta "about two hours"]
    python3 -m planetgen.cli.maintenance off
    python3 -m planetgen.cli.maintenance status

While it is on, every page and API call answers 503 to anyone but a
signed-in admin. The switch is the `maintenance.json` file beside
`settings.json` (`PLANETGEN_MAINTENANCE_FILE` overrides it), so it works
without a database. The admin page has the same switch.
"""

import argparse
import sys

from planetgen.admin import maintenance


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("on", "off", "status"))
    parser.add_argument("--eta", default="", help="When it should finish, in words (shown to visitors)")
    args = parser.parse_args(argv)
    if args.command == "on":
        try:
            maintenance.turn_on(args.eta, by="command line")
        except ValueError as exc:
            print(exc, file=sys.stderr)
            return 1
    elif args.command == "off":
        maintenance.turn_off(by="command line")
    state = maintenance.read()
    if state["enabled"]:
        print("Maintenance mode is ON" + (f" (expected finish: {state['eta']})" if state["eta"] else "")
              + (f", since {state['since']}" if state["since"] else "") + ".")
    else:
        print("Maintenance mode is off.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
