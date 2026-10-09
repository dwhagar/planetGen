#!/usr/bin/env python3
# planetgen.cli.warm_map

"""
Builds the Galaxy Map's opening view ahead of time (MAP.134), into the
tile cache the web site reads, so the first visitor after a release (or
after the cache was cleared) isn't the one who waits for it. update.sh
runs it in the background after every update; run it by hand after
clearing the tile cache. Run it as the web server's user so the files it
writes are readable by Apache.

    python3 -m planetgen.cli.warm_map

Prints how many tiles it built and how long that took. The work is the
same the first `/galaxy` visit does, so it stays correct however many
stars the galaxy holds; it never blocks an update because update.sh does
not wait for it.
"""

import sys

from planetgen.web.app import create_app
from planetgen.web.helpers import db_name
from planetgen.web.warmup import warm_opening_view


def main():
    app = create_app()
    with app.test_request_context("/galaxy"):
        try:
            result = warm_opening_view(db_name())
        except Exception as exc:  # noqa: BLE001 -- a failed warm-up only means the first visit builds the view
            print(f"warning: the Galaxy Map's opening view was not built ({exc}).", file=sys.stderr)
            sys.exit(1)
    print(f"Galaxy Map opening view: {result['tiles']} tiles ({result['cached']} already cached) "
          f"in {result['seconds']:.1f} s.")


if __name__ == "__main__":
    main()
