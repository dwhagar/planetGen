#!/usr/bin/env python3
# planetgen.cli.config

"""
Checks, documents and edits the settings.

    python3 -m planetgen.cli.config check                      # validate config.json
    python3 -m planetgen.cli.config docs                       # rewrite the files generated from the model
    python3 -m planetgen.cli.config get [OPTION]               # the value that applies, and where it comes from
    python3 -m planetgen.cli.config set OPTION=VALUE [...]     # save options in config.json
    python3 -m planetgen.cli.config unset OPTION [...]         # drop them, so the default applies
    python3 -m planetgen.cli.config restore BACKUP             # put the configuration of a backup archive back
                                                               # (`planetgen.cli.backup list` shows the archives)

`config.json` is the one file that holds all program configuration: your own copy of
`config.json.example`, which nothing ever overwrites. `check` reads it the way the app does, with unknown keys as
errors, and exits 1 naming every field that doesn't fit.
`docs` rewrites the option table in `docs/reference/config.md` (between its
GENERATED markers), `config.json.example` and `config.schema.json` from
`planetgen.util.settings.Settings`; a test fails when they differ from
what it would write.

`set`, `unset` and `restore` use the same writer as the configuration page
(`planetgen.util.settings_writer`): every value is checked against the model, an
the old file is backed up first and each change is written to the audit log. Whoever
can run the command owns the server, so every option may be set here, the database,
the secrets and the Owner-only ones included; the Admin page cannot. If `config.json`
does not exist the command says so and creates nothing. A value is JSON when
it parses as JSON for its type and text otherwise (`set site_name="My galaxy"`,
`set tile_cache.max_mb=300`, `set debug=true`); `OPTION=-` reads the value from
standard input, which keeps a secret out of the shell history.
"""

import argparse
import getpass
import json
import os
import sys

from planetgen.admin import activity_log, auth
from planetgen.util import logpaths, settings, settings_writer

DOCS_PATH = os.path.join(logpaths.PROJECT_ROOT, "docs", "reference", "config.md")
EXAMPLE_PATH = os.path.join(logpaths.PROJECT_ROOT, "config.json.example")
SCHEMA_PATH = os.path.join(logpaths.PROJECT_ROOT, "config.schema.json")
BEGIN = "<!-- BEGIN GENERATED OPTIONS (python -m planetgen.cli.config docs) -->"
END = "<!-- END GENERATED OPTIONS -->"


def render_docs(current):
    """`docs/reference/config.md` with the generated table between its markers
    replaced; raises `ValueError` when the markers are missing."""
    start, end = current.find(BEGIN), current.find(END)
    if start < 0 or end < start:
        raise ValueError(f"docs/reference/config.md needs the {BEGIN!r} and {END!r} lines")
    return current[:start + len(BEGIN)] + "\n\n" + settings.docs_table() + "\n" + current[end:]


def generated_files():
    """`{path: content}` of every generated file as it should be."""
    with open(DOCS_PATH, "r", encoding="utf-8") as f:
        docs = render_docs(f.read())
    return {DOCS_PATH: docs, EXAMPLE_PATH: settings.example_json(), SCHEMA_PATH: settings.schema_json()}


def run_check():
    try:
        config = settings._read_json_object(logpaths.CONFIG_PATH, "config.json")
        settings.build(config, strict=True)
    except settings.SettingsError as exc:
        print(f"Settings problem:\n{exc}", file=sys.stderr)
        return 1
    print("config.json is valid." if os.path.isfile(logpaths.CONFIG_PATH)
          else "There is no config.json: the defaults apply (copy config.json.example to config.json to change them).")
    for stray in unused_files():
        print(f"note: {stray} is not read any more; config.json is the one file for all configuration "
              "(move anything you still need into it, then delete the file).")
    return 0


def unused_files():
    """The `settings.json` files of the dropped overlay design that still exist: ignored, but named by `check`."""
    folders = {logpaths.PROJECT_ROOT, os.path.dirname(os.path.abspath(logpaths.CONFIG_PATH)), settings.state_dir(),
               "/var/lib/planetGen"}
    return sorted(path for path in (os.path.join(folder, "settings.json") for folder in folders)
                  if os.path.isfile(path))


def run_docs():
    for path, content in generated_files().items():
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"wrote {os.path.relpath(path, logpaths.PROJECT_ROOT)}")
    return 0


def _actor():
    try:
        return f"cli:{getpass.getuser()}"
    except Exception:  # noqa: BLE001
        return "cli"


def _parse_value(option, text):
    """`text` as the value of `option`: text for a text option, else JSON when it parses (a number, `true`, a
    list) and the text itself when not (a comma-separated list, a value the model will then refuse)."""
    fields = {path: meta for path, _info, meta in settings.iter_fields()}
    meta = fields.get(option)
    if meta is not None and meta["annotation"] is not str:
        try:
            return json.loads(text)
        except ValueError:
            return text
    return text


def _read_assignments(items):
    changes = {}
    for item in items:
        option, sep, text = item.partition("=")
        if not sep or not option:
            raise SystemExit(f"error: {item!r} is not OPTION=VALUE")
        if text == "-":
            text = sys.stdin.read().rstrip("\n")
        changes[option] = _parse_value(option, text)
    return changes


def _open_audit(args):
    """`(audit function, closer)`: writes `admin_audit_log` rows when the control database can be reached, and
    says so on standard error when it can't (the activity log gets the line either way)."""
    try:
        from planetgen.db.store import control_mysql_config, get_control_connection, mysql_config_from_args
        conn = get_control_connection(control_mysql_config(mysql_config_from_args(args)))
    except Exception as exc:  # noqa: BLE001 -- pymysql errors, a missing config, anything: the save still goes ahead
        print(f"warning: the control database could not be reached ({exc}); the change is in the activity "
              "log but not in the admin audit table.", file=sys.stderr)
        return None, lambda: None
    who = _actor()

    def audit(action, target, detail):
        auth.record_audit(conn, None, who, action, target=target, detail=detail)

    return audit, conn.close


def _show_result(result, verb):
    for option, message in result.errors.items():
        print(f"error: {option + ': ' if option else ''}{message}", file=sys.stderr)
    if result.errors:
        return 1
    if not result.changes:
        print("Nothing to change.")
        return 0
    for change in result.changes:
        print(f"{change.path}: {change.shown('old')!r} -> {change.shown('new')!r}")
    for option, variable in result.shadowed.items():
        print(f"note: {variable} is set in this environment and overrides {option}.")
    if result.backup:
        print(f"{verb}. Backup: {result.backup}")
    else:
        print(f"{verb}.")
    if result.restart:
        print("These need a restart to take effect: " + ", ".join(result.restart))
    for line in result.audit_errors:
        print(f"warning: audit row not written: {line}", file=sys.stderr)
    return 0


def _source(option, meta, config, environ):
    variable = meta.get("x-env")
    if variable and (environ.get(variable) or (meta.get("x-env-empty") and variable in environ)):
        return f"environment ({variable})"
    return "config.json" if settings_writer._get(config, option)[0] else "default"


def run_get(option):
    try:
        config = settings._read_json_object(logpaths.CONFIG_PATH, "config.json")
        current = settings.get_settings()
    except (settings.SettingsError, settings_writer.SettingsWriteError) as exc:
        print(f"Settings problem:\n{exc}", file=sys.stderr)
        return 1
    rows = [(path, meta) for path, _info, meta in settings.iter_fields() if option in (None, path)]
    if not rows:
        print(f"error: {option} is not a setting", file=sys.stderr)
        return 1
    for path, meta in rows:
        value = settings._lookup(current, path)
        shown = ("(set)" if value else "(empty)") if meta["x-secret"] else json.dumps(value)
        print(f"{path} = {shown}    [{_source(path, meta, config, os.environ)}]")
    return 0


def run_save(args, changes, verb):
    audit, close = _open_audit(args)
    try:
        result = settings_writer.save(changes, _actor(), role="host", audit=audit)
    except settings_writer.SettingsWriteError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    finally:
        close()
    return _show_result(result, verb)


def run_restore(args):
    audit, close = _open_audit(args)
    try:
        result = settings_writer.restore(args.backup, _actor(), role="host", audit=audit)
    except settings_writer.SettingsWriteError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    finally:
        close()
    return _show_result(result, "Restored")


def main(argv=None):
    from planetgen.db.store import add_mysql_connection_args
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("check", help="validate config.json")
    commands.add_parser("docs", help="rewrite the generated files")
    get = commands.add_parser("get", help="show the value that applies and where it comes from")
    get.add_argument("option", nargs="?")
    save = commands.add_parser("set", help="save options in config.json")
    save.add_argument("assignments", nargs="+", metavar="OPTION=VALUE")
    unset = commands.add_parser("unset", help="remove options from config.json (the default applies)")
    unset.add_argument("options", nargs="+", metavar="OPTION")
    restore = commands.add_parser("restore", help="put the configuration of a backup archive back")
    restore.add_argument("backup")
    for sub in (save, unset, restore):
        add_mysql_connection_args(sub)
    args = parser.parse_args(argv)
    if args.command == "check":
        return run_check()
    if args.command == "docs":
        return run_docs()
    if args.command == "get":
        return run_get(args.option)
    if args.command == "set":
        return run_save(args, _read_assignments(args.assignments), "Saved")
    if args.command == "unset":
        return run_save(args, {option: settings_writer.RESET for option in args.options}, "Saved")
    return run_restore(args)


if __name__ == "__main__":
    sys.exit(main())
