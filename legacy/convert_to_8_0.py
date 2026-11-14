#!/usr/bin/env python3
# legacy/convert_to_8_0.py

"""
Converts the data of the last release before 8.0 to the 8.0 form (OPS.177, Boss's
hard cut of 2026-11-12: no generation before 8.0 is compatible with 8.0 and above).

Run it by hand, on purpose, and only when you have data made before 8.0 that you
want to keep; nothing else in planetGen calls it. It does one thing:

* The generated data lives in the MySQL content database. This runs
  `bump_database_to_8_0.py` on it (DB.126), which rewrites its structure and the
  values the old steps moved, up to the oldest structure 8.0 and above accept. Take the
  same arguments, so `--status` shows the steps and `--yes` skips the question.
* `config.json` needs no conversion: a settings file without a `config_version` reads as
  version 1 (`planetgen.util.settings`), and every key of 8.0 has a default.

Systems drawn before 8.0 are not redrawn: the generator was reworked, so the same seed
now gives a different system (8.0 release notes). The converted rows keep what was stored.

Usage:
    python3 legacy/convert_to_8_0.py --status [--mysql-database NAME]
    python3 legacy/convert_to_8_0.py --yes    [--mysql-host H --mysql-port P --mysql-user U
                                               --mysql-password W --mysql-database NAME]
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bump_database_to_8_0  # noqa: E402


def main(argv=None):
    return bump_database_to_8_0.main(argv)


if __name__ == "__main__":
    sys.exit(main())
