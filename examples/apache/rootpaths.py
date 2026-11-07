"""
examples/apache/rootpaths.py

Which folders and files a root script (install.sh, update.sh and the scripts in this folder) may create, chown, move
files into, or write into a root-owned config for (SEC.44).

`config.json` is the web user's file: the Settings page rewrites it. So a folder it names (`tile_cache.dir`, `jobs.dir`,
`log_file`, `log_dir`) is a value an attacker who took over the web user can choose, and the scripts run as root,
make that folder, give it to the web user and write its name into /etc/logrotate.d. `check` is the one rule a path
read from config.json must pass first:

- only the plain characters of a path (letters, digits and `_ . / + @ : -`): no space, quote, newline, brace or glob,
  because the name goes into a logrotate or newsyslog file;
- absolute, already normal (no `..`, no `//`, no `.`), under one of the allowed parent folders (`ALLOWED_PARENTS`, the
  checkout, and anything the operator lists in `PLANETGEN_ROOT_PATH_PARENTS`, which only root's own environment can set),
  and strictly below it, not in one of the system folders under them (`PROTECTED`);
- no symbolic link on the way, and none at the end;
- every existing folder between the allowed parent and the path is root's (under /home and /Users: anyone's but the web
  user's, since a user may install in a home folder) and not writable by the web user (so the web user cannot swap a
  link in after the check); no name on the way starts with a dot (`~/.ssh`);
- an existing path is only taken over when it is already the web user's, or root's with the web group (what an earlier run
  left), or empty. /var/lib/mysql is none of those, so it is never re-owned.

A value that came from root's own environment (`PLANETGEN_TILE_CACHE_DIR`, `PLANETGEN_JOBS_DIR`, `PLANETGEN_LOG_FILE`,
`PLANETGEN_LOG_DIR`) is the operator's and is not checked: the web user cannot set root's environment (`SetEnv` in the
Apache vhost never reaches these scripts).

Standard library only: root runs this under `python3 -I`, and nothing is imported from the repo.
"""

import grp
import os
import pwd
import re
import stat

HOME_PARENTS = ("/home", "/Users")
"""tuple: Where users' home folders are (Linux, macOS). A planetGen installed in a home folder on a shared machine keeps its
folders there; the folders on the way are the user's, not root's, so they only have to be closed to the web user."""

ALLOWED_PARENTS = ("/var/lib", "/var/cache", "/var/log", "/srv", "/opt", "/usr/local/var",
                   "/Library/Logs", "/Library/Caches", "/Library/Application Support") + HOME_PARENTS
"""tuple: The folders a path from config.json may sit under (the checkout is allowed too)."""

PARENTS_ENV = "PLANETGEN_ROOT_PATH_PARENTS"
"""str: Extra allowed parents, separated by `os.pathsep`, from root's own environment."""

PROTECTED = frozenset({
    "mysql", "mariadb", "postgresql", "redis", "dpkg", "apt", "systemd", "docker", "containers", "snapd", "ssh",
    "sudo", "apache2", "nginx", "httpd", "journal", "private", "chrony", "logrotate", "cron", "ufw", "fail2ban",
    "auth.log", "syslog", "wtmp", "btmp", "lastlog", "kern.log", "messages", "secure", "dmesg",
})
"""frozenset: A name that must not be the first one under an allowed parent (a system service's folder or log)."""

PLAIN = re.compile(r"[A-Za-z0-9_./+@:\-]+")

ENV_OF = {"cache": "PLANETGEN_TILE_CACHE_DIR", "jobs": "PLANETGEN_JOBS_DIR", "log_file": "PLANETGEN_LOG_FILE",
          "log_dir": "PLANETGEN_LOG_DIR"}
"""dict: Which environment variable, when it holds the path, makes it the operator's choice."""

TRUSTED_UID = 0
"""int: Who owns the folders on the way: root. (Tests that do not run as root change it.)"""


class UnsafePath(Exception):
    """A path root will not use; the message says why, for the admin."""


def web_identity(user, group):
    """`(uid, gid)` of the web server's user and group; `None` for one that does not exist."""
    try:
        uid = pwd.getpwnam(user).pw_uid if user else None
    except KeyError:
        uid = None
    try:
        gid = grp.getgrnam(group).gr_gid if group else None
    except KeyError:
        gid = None
    return uid, gid


def _parents(extra):
    found = list(ALLOWED_PARENTS) + [p for p in extra if p]
    found += [p for p in os.environ.get(PARENTS_ENV, "").split(os.pathsep) if p]
    return [p.rstrip("/") or "/" for p in found]


def _web_gids(uid, gid):
    """Every group id the web server's user acts with: its main group, the one named, and any group that lists it."""
    found = {gid} if gid is not None else set()
    if uid is None:
        return found
    try:
        entry = pwd.getpwuid(uid)
    except KeyError:
        return found
    found.add(entry.pw_gid)
    found.update(g.gr_gid for g in grp.getgrall() if entry.pw_name in g.gr_mem)
    return found


def _web_can_write(info, uid, gid):
    """Whether the web server's user can write a folder: anyone can, or its group can."""
    mode = stat.S_IMODE(info.st_mode)
    return bool(mode & 0o002 or (mode & 0o020 and info.st_gid in _web_gids(uid, gid)))


def _is_empty(path):
    try:
        return not os.listdir(path)
    except OSError:
        return False


def check(path, label, uid=None, gid=None, extra_parents=()):
    """
    Returns `path` when root may use it, else raises `UnsafePath`.

    Args:
        path (str): A folder (`label` "cache", "jobs" or "log_dir") or a file ("log_file") read from config.json.
        label (str): What it is, for the message.
        uid, gid (int, optional): The web server's user and group (`web_identity`). Without them an existing folder
            is not judged by who owns it (a path only appended to, like update.log).
        extra_parents (iterable): Further allowed parents, such as the checkout.
    """
    trusted = os.environ.get(ENV_OF.get(label, ""))
    if trusted and trusted == path:
        return path
    if not isinstance(path, str) or not path or not PLAIN.fullmatch(path):
        raise UnsafePath(f"{label}: {path!r} has characters other than letters, digits and _ . / + @ : -")
    if not path.startswith("/") or (len(path) > 1 and os.path.normpath(path) != path):
        raise UnsafePath(f"{label}: {path} must be an absolute path with no '..', '.' or '//'")
    parent = next((p for p in sorted(_parents(extra_parents), key=len, reverse=True)
                   if path.startswith(p.rstrip("/") + "/")), None)
    if parent is None:
        raise UnsafePath(f"{label}: {path} is not under a folder the install scripts may use ({', '.join(ALLOWED_PARENTS)}"
                         f", the checkout, or those in {PARENTS_ENV}); put it there, or set the {ENV_OF.get(label, '')} "
                         f"environment variable in root's shell")
    relative = path[len(parent.rstrip("/")) + 1:]
    if any(part.startswith(".") for part in relative.split("/")):
        raise UnsafePath(f"{label}: {path} has a hidden folder or file name (one that starts with a dot)")
    home = parent in HOME_PARENTS
    if relative.split("/")[0] in PROTECTED:
        raise UnsafePath(f"{label}: {path} is a system service's folder or log")
    real_parent = os.path.realpath(parent)
    if os.path.realpath(path) != os.path.join(real_parent, relative):
        raise UnsafePath(f"{label}: {path} goes through a symbolic link")
    walk = real_parent
    for part in relative.split("/")[:-1]:
        walk = os.path.join(walk, part)
        try:
            info = os.lstat(walk)
        except FileNotFoundError:
            break
        theirs = info.st_uid == TRUSTED_UID or (home and (uid is None or info.st_uid != uid))
        if not stat.S_ISDIR(info.st_mode) or not theirs or _web_can_write(info, uid, gid):
            raise UnsafePath(f"{label}: {walk}, on the way to {path}, is not a folder that only "
                             f"{'its owner and root' if home else 'root'} can change")
    try:
        info = os.lstat(path)
    except FileNotFoundError:
        return path
    if stat.S_ISLNK(info.st_mode):
        raise UnsafePath(f"{label}: {path} is a symbolic link")
    is_file = label == "log_file"
    if not (stat.S_ISREG(info.st_mode) if is_file else stat.S_ISDIR(info.st_mode)):
        raise UnsafePath(f"{label}: {path} exists and is not a {'file' if is_file else 'folder'}")
    if uid is None and gid is None:
        return path
    mine = uid is not None and info.st_uid == uid
    left_by_us = info.st_uid == TRUSTED_UID and gid is not None and info.st_gid == gid
    empty = info.st_size == 0 if is_file else _is_empty(path)
    if not (mine or left_by_us or empty):
        raise UnsafePath(f"{label}: {path} already exists, is not the web server's, and is not empty; it was not "
                         f"taken over. If it really is planetGen's, give it to the web user by hand first")
    return path
