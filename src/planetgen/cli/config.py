#!/usr/bin/env python3
# planetgen.cli.config

"""
Checks and documents the settings (ADM.42).

    python3 -m planetgen.cli.config check   # validate config.json and settings.json
    python3 -m planetgen.cli.config docs    # rewrite the files generated from the model

`check` reads the same files the app does, with unknown keys in
`config.json` as errors, and exits 1 naming every field that doesn't fit.
`docs` rewrites the option table in `docs/config.md` (between its
GENERATED markers), `config.json.example` and `config.schema.json` from
`planetgen.util.settings.Settings`; a test fails when they differ from
what it would write.
"""

import argparse
import os
import sys

from planetgen.util import logpaths, settings

DOCS_PATH = os.path.join(logpaths.PROJECT_ROOT, "docs", "config.md")
EXAMPLE_PATH = os.path.join(logpaths.PROJECT_ROOT, "config.json.example")
SCHEMA_PATH = os.path.join(logpaths.PROJECT_ROOT, "config.schema.json")
BEGIN = "<!-- BEGIN GENERATED OPTIONS (python -m planetgen.cli.config docs) -->"
END = "<!-- END GENERATED OPTIONS -->"


def render_docs(current):
    """`docs/config.md` with the generated table between its markers
    replaced; raises `ValueError` when the markers are missing."""
    start, end = current.find(BEGIN), current.find(END)
    if start < 0 or end < start:
        raise ValueError(f"docs/config.md needs the {BEGIN!r} and {END!r} lines")
    return current[:start + len(BEGIN)] + "\n\n" + settings.docs_table() + "\n" + current[end:]


def generated_files():
    """`{path: content}` of every generated file as it should be."""
    with open(DOCS_PATH, "r", encoding="utf-8") as f:
        docs = render_docs(f.read())
    return {DOCS_PATH: docs, EXAMPLE_PATH: settings.example_json(), SCHEMA_PATH: settings.schema_json()}


def run_check():
    try:
        config = settings._read_json_object(logpaths.CONFIG_PATH, "config.json")
        overlay_path = settings.default_web_settings_path()
        overlay = settings._read_json_object(overlay_path, "settings.json") if os.path.isfile(overlay_path) else None
        settings.build(config, overlay, strict=True)
    except settings.SettingsError as exc:
        print(f"Settings problem:\n{exc}", file=sys.stderr)
        return 1
    print("config.json" + (" and settings.json are" if overlay is not None else " is") + " valid.")
    return 0


def run_docs():
    for path, content in generated_files().items():
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"wrote {os.path.relpath(path, logpaths.PROJECT_ROOT)}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("check", "docs"))
    args = parser.parse_args(argv)
    return run_check() if args.command == "check" else run_docs()


if __name__ == "__main__":
    sys.exit(main())
