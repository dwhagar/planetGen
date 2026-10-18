#!/usr/bin/env python3
"""
Regenerates, or checks, the global documents (DOC.69): the pages built from the whole repository, which any lane's
merge can make stale and which only the TODO thread rebuilds.

Boss (2026-10-18, 11:24 PM EDT): "TODO should be the only one who regenerates global documents." A lane's pull
request never contains a regenerated global file, and a stale one never fails a lane's gate; the TODO thread runs this
after each batch of merges and the release stamp checks it.

The global documents are:

* `docs/design/architecture-map.html` (`scripts/build_architecture_map.py`): lists every module of the package;
* `docs/plan/*.html` (`scripts/build_todo_docs.py`): the plan pages, from `docs/TODO.md` and the notes beside it;
* `docs/index.html` (`scripts/build_docs_index.py`): every HTML page under `docs/` with its description, so it is
  built last.

The tests that compare each with its source (`pytest.mark.globaldocs`, `src/tests/conftest.py`) are left out of every
run unless `PLANETGEN_GLOBAL_DOCS=1`; `scripts/test.py ... --global-docs` sets it, and the release stamp runs
`scripts/test.py full --check --global-docs`.

Usage
=====
    python scripts/build_global_docs.py            # rewrite every global document (TODO's step after a merge batch)
    python scripts/build_global_docs.py --check    # fail, naming each one, if any is out of date
"""

import argparse
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")

GENERATORS = (
    ("docs/design/architecture-map.html", "build_architecture_map"),
    ("docs/plan/*.html", "build_todo_docs"),
    ("docs/index.html", "build_docs_index"),
)
"""tuple: `(what it writes, script under scripts/)` in the order they run: the index last, as it lists the others."""


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SCRIPTS, f"{name}.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(check=False, generators=GENERATORS, load=_load):
    """
    Runs every generator (`main()`, or `main(["--check"])` with `check`), each to the end, and returns the list of
    what is out of date (with `check`) or failed to build (without): empty when every one is current.

    Args:
        check (bool): Only compare, write nothing.
        generators (tuple): `(what it writes, script name)` pairs.
        load (callable): Script name to its module.
    """
    failed = []
    for written, name in generators:
        code = load(name).main(["--check"] if check else [])
        if code:
            failed.append(written)
    return failed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="fail if any global document is out of date")
    args = parser.parse_args(argv)
    failed = run(check=args.check)
    if failed:
        what = "out of date" if args.check else "not built"
        print(f"global documents {what}: {', '.join(failed)}; run python scripts/build_global_docs.py")
        return 1
    print("global documents are up to date" if args.check else "global documents rebuilt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
