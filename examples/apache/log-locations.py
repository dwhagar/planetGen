"""
examples/apache/log-locations.py

Sets up planetGen's two log files on Linux and macOS (OPS.5), for
setup-debug-log.sh, which install.sh and update.sh run as root:

- the debug log: `PLANETGEN_LOG_FILE`, else `"log_file"` in config.json,
  else /var/log/planetgen.log (`appconfig.log_file_path`). Prepared
  whether or not `debug` is on, so turning it on later needs nothing
  else: its folder is made when it's missing, and the file itself is
  created, owned by the web server's user and group, mode 0660.
- the always-on activity log: `planetgen.log` in `PLANETGEN_LOG_DIR`,
  else `"log_dir"`, else the platform's folder (`appconfig.
  activity_log_path`). Its folder is owned by root (whoever runs this) and the web server's
  group, mode 2770 (the setgid bit keeps new files, the rotated ones
  among them, in that group), and the file is 0660 like the debug log.

Then it checks that the web server's user can really write each file
(every folder above it must let that user through) and that its group,
which CLI users who run the generator join, can too.

Anything it can't do (no rights, a path on a read-only or missing drive,
a folder in the way that is a file, a web server user that doesn't exist
yet) is a warning, never a failure (Boss, 2026-10-01): it prints what is
wrong and the exact commands that fix it, or how to point `log_file` and
`log_dir` somewhere writable, and exits 0. The program itself still
falls back as before when it can't open a log.

Like deploy-paths.py it runs as root, so it is run with `python3 -I` and
loads appconfig.py straight from its file (standard library only).

Usage:
    python3 -I examples/apache/log-locations.py <repo-dir> <user> <group>
"""

import grp
import importlib.util
import os
import pwd
import stat
import sys

IS_MACOS = sys.platform == "darwin"
UPDATE_HINT = "sudo ./update.sh"


def load_appconfig(repo_dir):
    path = os.path.join(repo_dir, "src", "planetgen", "util", "appconfig.py")
    spec = importlib.util.spec_from_file_location("planetgen_appconfig", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_of(env_var, config, key, defaults):
    """Where a log setting came from, for the report (`load_config` fills
    in `defaults`, so a value equal to its default reads as the default)."""
    if os.environ.get(env_var):
        return f"from {env_var}"
    if config.get(key) and config.get(key) != defaults.get(key):
        return f'from "{key}" in config.json'
    return "the default"


def lookup_identity(user, group):
    """`(uid, gid, the user's group ids)`, with `None` for a uid or gid
    that doesn't exist on this machine."""
    try:
        entry = pwd.getpwnam(user)
    except KeyError:
        entry = None
    try:
        gid = grp.getgrnam(group).gr_gid
    except KeyError:
        gid = None
    if entry is None:
        return None, gid, set()
    gids = {entry.pw_gid}
    gids.update(g.gr_gid for g in grp.getgrall() if user in g.gr_mem)
    return entry.pw_uid, gid, gids


def _allowed(info, uid, gids, bits):
    """Whether `uid` (with `gids`) gets `bits` (a mix of 4/2/1) on a file
    with `info`, by its owner/group/other mode bits. root always does."""
    if uid == 0:
        return True
    mode = stat.S_IMODE(info.st_mode)
    if uid is not None and info.st_uid == uid:
        granted = (mode >> 6) & 7
    elif info.st_gid in gids:
        granted = (mode >> 3) & 7
    else:
        granted = mode & 7
    return granted & bits == bits


def can_write(path, uid, gids):
    """
    Whether a user can append to the file `path`: every folder above it
    must let them through (x), and the file itself must be writable.
    `uid=None` asks for the group alone (a CLI user in that group).
    Returns `(ok, the path that stops them)`.
    """
    path = os.path.abspath(path)
    folder = os.path.dirname(path)
    parts = []
    current = folder
    while True:
        parts.append(current)
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    for part in reversed(parts):
        try:
            info = os.stat(part)
        except OSError:
            return False, part
        if not _allowed(info, uid, gids, 1):
            return False, part
    try:
        info = os.stat(path)
    except OSError:
        return False, path
    ok = _allowed(info, uid, gids, 2)
    return ok, (None if ok else path)


class LogSetup:
    """Collects what one run did, and what it couldn't, for the report."""

    def __init__(self, user, group):
        self.user = user
        self.group = group
        self.uid, self.gid, self.gids = lookup_identity(user, group)
        self.warnings = 0

    def warn(self, what, problems, commands, setting):
        self.warnings += 1
        err = sys.stderr
        print(f"warning: {what} isn't fully set up:", file=err)
        for problem in problems:
            print(f"  - {problem}", file=err)
        print("  The install/update carries on; the site still runs, but this log may not be written.", file=err)
        print("  Fix it with:", file=err)
        for command in commands:
            print(f"    {command}", file=err)
        print(f"  or point {setting} somewhere this server can write, then run {UPDATE_HINT} again.", file=err)

    _ACTIONS = {"makedirs": "create the folder", "_touch": "create", "chown": "set the owner of",
                "chmod": "set the mode of"}

    def _try(self, problems, action, *args):
        try:
            action(*args)
            return True
        except OSError as exc:
            verb = self._ACTIONS.get(action.__name__, action.__name__)
            problems.append(f"couldn't {verb} {args[0]}: {exc.strerror or exc}")
            return False

    def _owner(self, path, uid, gid, problems):
        if uid is None or gid is None:
            return False
        return self._try(problems, os.chown, path, uid, gid)

    def _make_folder(self, folder, problems):
        """Creates `folder` (and its parents) unless something on the way
        is a file rather than a folder, which is reported as such."""
        current = folder
        while not os.path.exists(current):
            parent = os.path.dirname(current)
            if parent == current:
                break
            current = parent
        if os.path.exists(current) and not os.path.isdir(current):
            problems.append(f"{current} is a file, not a folder")
            return
        self._try(problems, os.makedirs, folder, 0o755)

    def _identity_problems(self):
        problems = []
        if self.uid is None:
            hint = "" if IS_MACOS else " (is the web server installed? sudo apt install apache2)"
            problems.append(f"the web server's user '{self.user}' doesn't exist on this machine{hint}")
        if self.gid is None:
            problems.append(f"the web server's group '{self.group}' doesn't exist on this machine")
        return problems

    def _check_access(self, path, problems):
        if self.uid is None or self.gid is None:
            return
        ok, blocker = can_write(path, self.uid, self.gids)
        if not ok:
            problems.append(f"'{self.user}' can't write {path} ({blocker or path} is in the way)")
        ok, blocker = can_write(path, None, {self.gid})
        if not ok:
            problems.append(f"the '{self.group}' group (CLI users) can't write {path} "
                            f"({blocker or path} is in the way)")

    def _folder_commands(self, folder, path):
        """`sudo chmod o+x` lines for any folder above `path` the user
        can't pass through (a 0700 /srv, say)."""
        commands = []
        if self.uid is None:
            return commands
        current = folder
        chain = []
        while True:
            chain.append(current)
            parent = os.path.dirname(current)
            if parent == current:
                break
            current = parent
        for part in reversed(chain):
            try:
                info = os.stat(part)
            except OSError:
                break
            if stat.S_ISDIR(info.st_mode) and not _allowed(info, self.uid, self.gids, 1):
                commands.append(f"sudo chmod o+x {_quote(part)}")
        return commands

    def debug_log(self, path, source, debug_on):
        path = os.path.abspath(path)
        folder = os.path.dirname(path)
        problems = self._identity_problems()
        if not os.path.isdir(folder):
            self._make_folder(folder, problems)
        if os.path.isdir(folder) and not os.path.exists(path):
            self._try(problems, _touch, path)
        if os.path.isfile(path):
            if self._owner(path, self.uid, self.gid, problems):
                self._try(problems, os.chmod, path, 0o660)
        elif os.path.exists(path):
            problems.append(f"{path} exists but isn't a file")
        if not problems:
            self._check_access(path, problems)
        state = "on" if debug_on else "off; nothing is written until it's on"
        if problems:
            commands = [f"sudo mkdir -p {_quote(folder)}", f"sudo touch {_quote(path)}",
                        f"sudo chown {self.user}:{self.group} {_quote(path)}", f"sudo chmod 0660 {_quote(path)}"]
            self.warn(f"the debug log {path} ({source}; debug is {state})", problems,
                      commands + self._folder_commands(folder, path),
                      '"log_file" in config.json (or PLANETGEN_LOG_FILE)')
        else:
            print(f"Debug log: {path} ({source}; debug is {state})")
            print(f"  owned by {self.user}:{self.group}, mode 0660; add CLI users who run the generator "
                  f"to the {self.group} group so they can append")

    def activity_log(self, path, source):
        path = os.path.abspath(path)
        folder = os.path.dirname(path)
        problems = self._identity_problems()
        if not os.path.isdir(folder):
            self._make_folder(folder, problems)
        if os.path.isdir(folder):
            # Owned by whoever runs this (root, from install.sh/update.sh).
            if self._owner(folder, os.geteuid(), self.gid, problems):
                self._try(problems, os.chmod, folder, 0o2770)
            if not os.path.exists(path):
                self._try(problems, _touch, path)
        if os.path.isfile(path):
            if self._owner(path, self.uid, self.gid, problems):
                self._try(problems, os.chmod, path, 0o660)
        elif os.path.exists(path):
            problems.append(f"{path} exists but isn't a file")
        if not problems:
            self._check_access(path, problems)
        if problems:
            commands = [f"sudo mkdir -p {_quote(folder)}", f"sudo chown root:{self.group} {_quote(folder)}",
                        f"sudo chmod 2770 {_quote(folder)}", f"sudo touch {_quote(path)}",
                        f"sudo chown {self.user}:{self.group} {_quote(path)}", f"sudo chmod 0660 {_quote(path)}"]
            self.warn(f"the activity log {path} ({source})", problems,
                      commands + self._folder_commands(folder, path),
                      '"log_dir" in config.json (or PLANETGEN_LOG_DIR)')
        else:
            print(f"Activity log: {path} ({source}; folder root:{self.group} 2770, "
                  f"file {self.user}:{self.group} 0660)")


def _touch(path):
    with open(path, "a", encoding="utf-8"):
        pass
    os.chmod(path, 0o660)


def _quote(path):
    """A path as a shell word: quoted only when it needs it."""
    if path and all(c.isalnum() or c in "/._-+:@%" for c in path):
        return path
    return "'" + path.replace("'", "'\\''") + "'"


def main(argv):
    if len(argv) != 4:
        print("usage: python3 -I log-locations.py <repo-dir> <user> <group>", file=sys.stderr)
        return 2
    repo_dir, user, group = argv[1:]
    try:
        appconfig = load_appconfig(repo_dir)
        config = appconfig.load_config()
        debug_on = appconfig.debug_enabled(config)
        log_file = appconfig.log_file_path(config)
        activity = appconfig.activity_log_path(config)
    except Exception as exc:  # a broken config.json: say so, never stop the install
        print(f"warning: couldn't work out where the logs go ({type(exc).__name__}: {exc}).", file=sys.stderr)
        print(f"  Check config.json's \"log_file\" and \"log_dir\" (each a path, or leave them out for the "
              f"defaults), then run {UPDATE_HINT} again.", file=sys.stderr)
        return 0
    setup = LogSetup(user, group)
    setup.debug_log(log_file, source_of("PLANETGEN_LOG_FILE", config, "log_file", appconfig.DEFAULT_CONFIG), debug_on)
    setup.activity_log(activity, source_of("PLANETGEN_LOG_DIR", config, "log_dir", appconfig.DEFAULT_CONFIG))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
