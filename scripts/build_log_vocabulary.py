#!/usr/bin/env python3
"""
Writes the vocabulary block of `scripts/logrecord.sh` from the logging module (OPS.104): the severities, the ones that go to standard error,
the word of a stage and the divider (`planetgen.util.logrecord.shell_vocabulary`).

    python scripts/build_log_vocabulary.py            # write the block
    python scripts/build_log_vocabulary.py --check    # fail if the committed block differs

The shell cannot import Python before the virtual environment exists (`install.sh` runs first), so `scripts/logrecord.sh` stays the one shell
implementation of the display line; this keeps its words the module's. `src/tests/test_log_shell.py` fails when the block is out of date.
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SHELL = os.path.join(HERE, "logrecord.sh")
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from planetgen.util import logrecord  # noqa: E402


def rebuilt(text):
    """`text` (the shell file) with its generated block replaced by the module's."""
    start = text.index(logrecord.SHELL_BEGIN)
    end = text.index(logrecord.SHELL_END) + len(logrecord.SHELL_END)
    return text[:start] + logrecord.shell_vocabulary() + text[end:]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="build in memory and fail if the committed file differs")
    args = parser.parse_args()
    with open(SHELL, encoding="utf-8") as handle:
        text = handle.read()
    fresh = rebuilt(text)
    if args.check:
        if fresh != text:
            print("scripts/logrecord.sh: the vocabulary block is out of date; run python scripts/build_log_vocabulary.py", file=sys.stderr)
            return 1
        return 0
    with open(SHELL, "w", encoding="utf-8") as handle:
        handle.write(fresh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
