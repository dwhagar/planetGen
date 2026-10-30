#!/usr/bin/env python3
"""
The README's list of external files (DOC.72): every dataset the setup step downloads, with its source, its files and their URLs, how to
fetch it, how to install it by hand, and what stops working while it is missing. Written from the catalogue (`planetgen.datasets.catalogue`),
so a dataset added later is in the README the next time this runs, and `src/tests/test_readme_external_files.py` fails until it has been.

    python3 scripts/readme_external_files.py            # print the section
    python3 scripts/readme_external_files.py --write    # put it between the markers in README.md
    python3 scripts/readme_external_files.py --check    # exit 1 when README.md differs
"""

import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

START = "<!-- external-files:start (written by scripts/readme_external_files.py) -->"
END = "<!-- external-files:end -->"


def _megabytes(size):
    return f"{size / 1e6:.1f} MB" if size >= 1e5 else f"{size / 1e3:.0f} kB"


def _dataset_block(dataset):
    total = sum(piece.size for piece in dataset.pieces)
    lines = [f"### `{dataset.key}`: {dataset.title}", "",
             f"- Source: {dataset.source}.",
             f"- Download: {len(dataset.pieces)} file{'s' if len(dataset.pieces) != 1 else ''}, {_megabytes(total)}. "
             f"{'Kept after setup: ' + ', '.join(f'`{name}`' for name in dataset.derived) + ' (the download is deleted once it is reduced).' if dataset.derived else 'All kept as downloaded.'}",
             f"- Fetch: `python3 -m planetgen.cli.datasets fetch {dataset.key}`.",
             f"- Replaces: {dataset.replaces}.",
             f"- With `--disregard-missing-data`: {dataset.fallback[0].upper() + dataset.fallback[1:]}." if dataset.fallback
             else "- With `--disregard-missing-data`: nothing falls back, the run is still refused.",
             "- Files (URL, then the file name it is kept under in the dataset's folder):"]
    lines += [f"  - <{piece.url}> as `{piece.name}`" for piece in dataset.pieces]
    return lines


def section():
    from planetgen.datasets import catalogue

    lines = [START, ""]
    for dataset in catalogue.ALL:
        lines += _dataset_block(dataset) + [""]
    lines += [END]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    text = section()
    readme = ROOT / "README.md"
    current = readme.read_text(encoding="utf-8")
    if START not in current or END not in current:
        print("README.md has no external-files markers", file=sys.stderr)
        return 1
    head, rest = current.split(START, 1)
    _old, tail = rest.split(END, 1)
    updated = head + text + tail
    if args.check:
        if updated != current:
            print("README.md's external files differ from the catalogue: run scripts/readme_external_files.py --write", file=sys.stderr)
            return 1
        return 0
    if args.write:
        readme.write_text(updated, encoding="utf-8")
        return 0
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
