#!/usr/bin/env python3
# planetgen.cli.maintenance

"""
Turns maintenance mode on or off (UX.112).

    python3 -m planetgen.cli.maintenance on [--eta "about two hours"]
    python3 -m planetgen.cli.maintenance off
    python3 -m planetgen.cli.maintenance status

While it is on, every page and API call answers 503 to anyone but a
signed-in admin. The switch is the `maintenance.json` file in the state
folder, `/var/lib/planetGen/state` (`PLANETGEN_STATE_DIR` or
`PLANETGEN_MAINTENANCE_FILE` override it), so it works without a database. The admin page has the same switch.
"""

import argparse
import sys

from planetgen.admin import maintenance
from planetgen.util import log


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("on", "off", "status"))
    parser.add_argument("--eta", default="", help="When it should finish, in words (shown to visitors)")
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser)
    try:
        if args.command == "on":
            maintenance.turn_on(args.eta, by="command line")
        elif args.command == "off":
            maintenance.turn_off(by="command line")
    except (ValueError, maintenance.MaintenanceWriteError) as exc:
        print(exc, file=sys.stderr)
        return 1
    state = maintenance.read()
    if state["enabled"]:
        print("Maintenance mode is ON" + (f" (expected finish: {state['eta']})" if state["eta"] else "")
              + (f", since {state['since']}" if state["since"] else "") + ".")
    else:
        print("Maintenance mode is off.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
