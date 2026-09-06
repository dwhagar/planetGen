# stellarObjects/_version.py

"""
Single Source of Truth for Version Information
================================================

Deliberately kept in its own tiny, dependency-free module rather than
directly in `stellarObjects/__init__.py`. `setup.py` needs to read the
version number at build time, but `stellarObjects/__init__.py` imports
`Planet`, which imports `names.py`, which imports `nltk` at module load time
-- and `setup.py` runs in an isolated build environment that doesn't have
`nltk` (or any other runtime dependency) installed yet (see `setup.py`'s own
comment on this). `setup.py` reads `__version__` out of this file as plain
text (regex, not import) to avoid ever triggering that chain; anything
importing this module normally (e.g. `stellarObjects/__init__.py`, or the
CLI scripts) does so safely, same as any other submodule.
"""

import argparse

__version__ = "5.1.0"

REPO_URL = "https://github.com/dwhagar/planetGen"

LICENSE_SUMMARY = "License: CC0 1.0 Universal (Public Domain Dedication) -- see LICENSE.md"


def version_banner(prog_name):
    """
    Builds the multi-line banner printed by each CLI script's `--version`
    flag: the program name and version, the project's repository URL, and a
    one-line license summary.

    Args:
        prog_name (str): The name of the CLI script to display (e.g.
                         "systemGen.py").

    Returns:
        str: The formatted, multi-line version banner.
    """
    return f"{prog_name} {__version__}\n{REPO_URL}\n{LICENSE_SUMMARY}"


class VersionAction(argparse.Action):
    """
    Prints a `--version` banner exactly as given and exits.

    Deliberately not argparse's own built-in `action='version'`: that action
    pipes the string through argparse's `HelpFormatter`, which treats it as
    flowing prose and rewraps it -- collapsing this project's intentional
    line breaks (program/version, repo URL, license summary each on their
    own line) into a re-wrapped paragraph. This action writes the banner
    straight to stdout instead, preserving it verbatim.
    """

    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS,
                 banner=None, help="show version, repository URL, and license, then exit"):
        super().__init__(option_strings=option_strings, dest=dest, default=default,
                         nargs=0, help=help)
        self.banner = banner

    def __call__(self, parser, namespace, values, option_string=None):
        print(self.banner)
        parser.exit()
