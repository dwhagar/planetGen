# planetgen/backup/core.py

"""
The backup archive (ADM.87)
===========================

Every backup is one `backups/<UTC time>.tar.bz2` taken at one moment, holding **parts**: the settings
configuration now (`settings/`, a copy of config.json), the control database dump (`control/`, ADM.88) and the DB repair parity data
(`parity/`, DB.40) in 9.0. A part is a small provider (`Part`: a name, `dump(directory)` and
`restore(directory, context)`) registered in one table, so a later part adds a provider and not a format.

    20261016T193012Z.tar.bz2
        manifest.json     {"format": 1, "taken_at", "taken_by", "reason", "version",
                           "parts": [{"name", "directory", "size", "files": [{"path", "size", "sha256"}]}]}
        settings/config.json

`take` writes one archive (a temporary file, `fsync`, rename; mode 0600), `list_archives` reads the manifests,
`verify` checks every hash, `restore` verifies first and restores the whole archive or the named parts, and
`prune` keeps the newest `settings.backups_kept` archives, never deleting the newest one that holds every
registered part. A restore takes a safety archive of the parts it is about to replace first. A part in an
archive that this build has no provider for is skipped and reported, and the archive keeps it: nothing here
rewrites an archive.

Pieces are read out of the archive by name; the archive is never extracted blindly: only the files the
manifest lists are written, to a temporary folder, after their hashes match, and names with a `..` or a
leading `/` are refused.

This module imports no Flask.
"""

import hashlib
import io
import json
import os
import re
import shutil
import tarfile
import tempfile
import threading
from dataclasses import dataclass, field
from datetime import datetime, timezone

try:  # POSIX; elsewhere the in-process lock below is all there is
    import fcntl
except ImportError:  # pragma: no cover
    fcntl = None

from planetgen._version import __version__
from planetgen.admin import activity_log

FORMAT = 1
MANIFEST = "manifest.json"
SUFFIX = ".tar.bz2"
ARCHIVE_MODE = 0o600
DEFAULT_KEPT = 20
"""int: How many archives are kept when `settings.backups_kept` can't be read."""

MAX_FILE_BYTES = 8 * 1024 ** 3
"""int: The largest single file read out of an archive (a guard against a crafted archive)."""

FIX_HINT = ("The backups folder must belong to the web server's user: run ./update.sh on the server "
            "(or `sudo examples/apache/create-cache-dir.sh`) and try again.")

_NAME = re.compile(r"^(\d{8}T\d{6}Z)(-\d+)?" + re.escape(SUFFIX) + r"$")
_PART = re.compile(r"^[a-z][a-z0-9_-]{0,31}$")
_FILE = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]{0,100}(/[A-Za-z0-9_][A-Za-z0-9_.-]{0,100}){0,6}$")
_lock = threading.Lock()


class BackupError(OSError):
    """A backup could not be taken, read or restored; the message says what and where."""


# ---------------------------------------------------------------------
# Parts
# ---------------------------------------------------------------------

@dataclass
class RestoreContext:
    """Who is restoring and how: handed to every part's `restore`."""
    who: str = None
    role: str = "admin"
    audit: object = None
    environ: object = None


class Part:
    """A provider for one kind of backed-up data. Subclasses set `name` (lower case, it is the archive's folder)
    and `label` (for a page), and implement `dump` and `restore`."""

    name = ""
    label = ""

    def dump(self, directory):
        """Writes this part's files into `directory` (an empty folder)."""
        raise NotImplementedError

    def restore(self, directory, context):
        """Puts the files in `directory` (a hash-checked copy of what `dump` wrote) back. May return a result
        object for the caller to show."""
        raise NotImplementedError


_PARTS = {}


def register(part):
    """Adds a provider (replacing one of the same name); returns it."""
    if not _PART.match(part.name):
        raise ValueError(f"{part.name!r} is not a usable part name")
    _PARTS[part.name] = part
    return part


def registered():
    """The registered parts, `{name: Part}`, in registration order."""
    return dict(_PARTS)


# ---------------------------------------------------------------------
# Folder, files and locks
# ---------------------------------------------------------------------

def backups_dir():
    """The folder the archives live in: `backups` in the web-owned state folder."""
    from planetgen.util import settings
    return os.path.join(settings.state_dir(), "backups")


def fsync_folder(folder):
    try:
        fd = os.open(folder, os.O_RDONLY)
    except OSError:  # not possible on some platforms; the rename itself is still atomic
        return
    try:
        os.fsync(fd)
    except OSError:
        pass
    finally:
        os.close(fd)


def match_owner(path, reference):
    """When running as root, gives `path` the owner and group of `reference` (the state folder belongs to the web
    server's user; a file root made there would otherwise be unreadable to the site)."""
    if os.name != "posix" or os.geteuid() != 0:
        return
    try:
        stat = os.stat(reference)
        os.chown(path, stat.st_uid, stat.st_gid)
    except OSError:
        pass


def atomic_write(path, data, mode=None):
    """Writes `data` (bytes) to `path` so a reader sees the old file or the whole new one: a temporary file in
    the same folder, `fsync`, the old file's mode (else `mode`, else 0640), `os.replace`, `fsync` of the folder.
    Raises `BackupError` naming the path."""
    folder = os.path.dirname(os.path.abspath(path))
    try:
        os.makedirs(folder, exist_ok=True)
        handle, temp = tempfile.mkstemp(prefix=".write-", dir=folder)
    except OSError as exc:
        raise BackupError(f"{path} can't be written ({exc.strerror or exc}). {FIX_HINT}") from exc
    try:
        with os.fdopen(handle, "wb") as out:
            out.write(data)
            out.flush()
            os.fsync(out.fileno())
        try:
            keep = os.stat(path).st_mode & 0o777
        except FileNotFoundError:
            keep = mode if mode is not None else 0o640
        os.chmod(temp, keep)
        match_owner(temp, folder)
        os.replace(temp, path)
    except BaseException as exc:
        if os.path.exists(temp):
            os.unlink(temp)
        if isinstance(exc, OSError):
            raise BackupError(f"{path} can't be written ({exc.strerror or exc}). {FIX_HINT}") from exc
        raise
    fsync_folder(folder)


class _Lock:
    """Exclusive over the backups folder: a lock file between processes, a thread lock inside one."""

    def __init__(self, folder):
        self.folder = folder
        self.handle = None

    def __enter__(self):
        _lock.acquire()
        try:
            if not os.path.isdir(self.folder):
                os.makedirs(self.folder, mode=0o700, exist_ok=True)
                match_owner(self.folder, os.path.dirname(self.folder))
            if fcntl is not None:
                lock_path = os.path.join(self.folder, ".lock")
                fresh = not os.path.exists(lock_path)
                self.handle = open(lock_path, "a")
                if fresh:
                    match_owner(lock_path, self.folder)
                fcntl.flock(self.handle, fcntl.LOCK_EX)
        except OSError as exc:
            self._close()
            _lock.release()
            raise BackupError(f"The backups folder {self.folder} can't be used ({exc.strerror or exc}). "
                              f"{FIX_HINT}") from exc
        except BaseException:
            self._close()
            _lock.release()
            raise
        return self

    def _close(self):
        if self.handle is not None:
            try:
                self.handle.close()
            finally:
                self.handle = None

    def __exit__(self, *exc_info):
        try:
            if self.handle is not None and fcntl is not None:
                fcntl.flock(self.handle, fcntl.LOCK_UN)
            self._close()
        finally:
            _lock.release()


def write_problem(folder=None):
    """Why this process can't write the backups folder, or `None` (a check made without writing)."""
    target = folder or backups_dir()
    while not os.path.isdir(target):
        parent = os.path.dirname(target)
        if not parent or parent == target:
            break
        target = parent
    if os.access(target, os.W_OK | os.X_OK):
        return None
    return f"{target} is not writable by this user. {FIX_HINT}"


# ---------------------------------------------------------------------
# Archives
# ---------------------------------------------------------------------

@dataclass
class PartInfo:
    name: str
    size: int
    files: list = field(default_factory=list)
    known: bool = True
    """False when this build has no provider for the part."""


@dataclass
class Archive:
    name: str
    path: str
    taken_at: str
    taken_by: str
    reason: str
    version: str
    parts: list
    size: int
    complete: bool = False
    """True when it holds every registered part."""
    readable: bool = True

    @property
    def part_names(self):
        return [p.name for p in self.parts]


@dataclass
class Verification:
    ok: bool
    problems: list = field(default_factory=list)
    manifest: dict = field(default_factory=dict)


@dataclass
class RestoreReport:
    archive: str
    restored: list = field(default_factory=list)
    skipped: list = field(default_factory=list)
    """Parts in the archive this build has no provider for: reported, never dropped."""
    results: dict = field(default_factory=dict)
    safety: str = None
    """The archive taken of the parts before they were replaced."""


def _utc(when=None):
    return (when or datetime.now(timezone.utc)).astimezone(timezone.utc)


def _hash_file(path):
    digest, size = hashlib.sha256(), 0
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
            size += len(block)
    return digest.hexdigest(), size


def _order(name):
    stamp, n = _NAME.match(name).groups()
    return stamp, int(n[1:]) if n else 1


def take(parts=None, who=None, reason="backup", directory=None, when=None, keep=None):
    """
    Writes one archive of `parts` (names; default every registered part) and prunes.

    Args:
        who (str): Who is making it (a user name, `cli:<user>`).
        reason (str): Why (`config.update`, `manual`, `before restore of ...`).
        keep (int): Archives to keep; default the option `settings.backups_kept`.

    Returns:
        Archive: The new archive.

    Raises:
        BackupError: A part is not registered, a part's `dump` failed, or the folder can't be written; nothing is
            left behind.
    """
    folder = directory or backups_dir()
    chosen = list(registered()) if parts is None else list(parts)
    unknown = [name for name in chosen if name not in _PARTS]
    if unknown or not chosen:
        raise BackupError("No backup part named " + ", ".join(repr(n) for n in unknown) if unknown
                          else "There is nothing to back up: no part is registered.")
    when = _utc(when)
    with _Lock(folder):
        work = tempfile.mkdtemp(prefix=".work-", dir=folder)
        temp = None
        try:
            manifest = {"format": FORMAT, "taken_at": when.strftime("%Y-%m-%dT%H:%M:%SZ"), "taken_by": who,
                        "reason": reason, "version": __version__, "parts": []}
            for name in chosen:
                part_dir = os.path.join(work, name)
                os.mkdir(part_dir)
                try:
                    _PARTS[name].dump(part_dir)
                except BackupError:
                    raise
                except Exception as exc:  # noqa: BLE001 -- name the part, leave nothing behind
                    raise BackupError(f"The {name} part could not be backed up: {exc}") from exc
                files = []
                for root, _dirs, names in os.walk(part_dir):
                    for filename in sorted(names):
                        full = os.path.join(root, filename)
                        digest, size = _hash_file(full)
                        files.append({"path": os.path.relpath(full, part_dir).replace(os.sep, "/"),
                                      "size": size, "sha256": digest})
                files.sort(key=lambda f: f["path"])
                manifest["parts"].append({"name": name, "directory": name, "size": sum(f["size"] for f in files),
                                          "files": files})
            handle, temp = tempfile.mkstemp(prefix=".archive-", dir=folder)
            os.close(handle)
            with tarfile.open(temp, "w:bz2") as tar:
                data = (json.dumps(manifest, indent=2) + "\n").encode("utf-8")
                info = tarfile.TarInfo(MANIFEST)
                info.size, info.mtime, info.mode = len(data), int(when.timestamp()), 0o600
                tar.addfile(info, io.BytesIO(data))
                for part in manifest["parts"]:
                    for item in part["files"]:
                        arc = f"{part['directory']}/{item['path']}"
                        info = tar.gettarinfo(os.path.join(work, part["directory"], *item["path"].split("/")), arc)
                        info.uid = info.gid = 0
                        info.uname = info.gname = ""
                        info.mode = 0o600
                        with open(os.path.join(work, part["directory"], *item["path"].split("/")), "rb") as src:
                            tar.addfile(info, src)
            with open(temp, "rb") as written:
                os.fsync(written.fileno())
            base = when.strftime("%Y%m%dT%H%M%SZ")
            n, name = 1, f"{base}{SUFFIX}"
            while os.path.exists(os.path.join(folder, name)):
                n += 1
                name = f"{base}-{n}{SUFFIX}"
            os.chmod(temp, ARCHIVE_MODE)
            match_owner(temp, folder)
            os.replace(temp, os.path.join(folder, name))
            temp = None
            fsync_folder(folder)
        except OSError as exc:
            if isinstance(exc, BackupError):
                raise
            raise BackupError(f"A backup could not be written in {folder} ({exc.strerror or exc}). "
                              f"{FIX_HINT}") from exc
        finally:
            shutil.rmtree(work, ignore_errors=True)
            if temp and os.path.exists(temp):
                os.unlink(temp)
        archive = _archive_row(folder, name)
        _prune(folder, keep)
    activity_log.event("DB", "backup.take", user=who, archive=name, reason=reason, parts=",".join(chosen))
    return archive


def _read_manifest(path):
    with tarfile.open(path, "r:bz2") as tar:
        member = next(iter(tar), None)
        if member is None or member.name != MANIFEST or not member.isfile() or member.size > 64 * 1024 * 1024:
            raise ValueError("the manifest is missing or not first")
        manifest = json.loads(tar.extractfile(member).read().decode("utf-8"))
    if not isinstance(manifest, dict) or manifest.get("format") != FORMAT or not isinstance(manifest.get("parts"), list):
        raise ValueError("the manifest is not format 1")
    return manifest


def _archive_row(folder, name):
    path = os.path.join(folder, name)
    stat = os.stat(path)
    try:
        manifest = _read_manifest(path)
        parts = [PartInfo(p["name"], int(p.get("size") or 0), list(p.get("files") or []), p["name"] in _PARTS)
                 for p in manifest["parts"] if isinstance(p, dict) and _PART.match(str(p.get("name", "")))]
        return Archive(name, path, str(manifest.get("taken_at") or ""), str(manifest.get("taken_by") or ""),
                       str(manifest.get("reason") or ""), str(manifest.get("version") or ""), parts, stat.st_size,
                       complete=set(_PARTS) <= {p.name for p in parts})
    except (OSError, tarfile.TarError, EOFError, ValueError, KeyError, TypeError):
        taken = datetime.fromtimestamp(stat.st_mtime, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        return Archive(name, path, taken, "", "", "", [], stat.st_size, readable=False)


def list_archives(directory=None):
    """The archives, newest first. One whose manifest can't be read still lists, `readable` false."""
    folder = directory or backups_dir()
    try:
        names = os.listdir(folder)
    except FileNotFoundError:
        return []
    except OSError as exc:
        raise BackupError(f"{folder} can't be listed: {exc.strerror or exc}") from exc
    rows = []
    for name in sorted((n for n in names if _NAME.match(n)), key=_order, reverse=True):
        try:
            rows.append(_archive_row(folder, name))
        except OSError:
            continue
    return rows


def _check_name(name):
    if not isinstance(name, str) or not _NAME.match(name):
        raise BackupError(f"There is no backup named {name!r}.")
    return name


def _safe_member(name):
    return bool(name) and not name.startswith("/") and ".." not in name.split("/") and "\\" not in name


def verify(name, directory=None):
    """Checks one archive: it opens, the manifest is format 1, every listed file is there with the recorded size and
    SHA-256, and nothing unlisted is in it."""
    folder = directory or backups_dir()
    path = os.path.join(folder, _check_name(name))
    problems = []
    try:
        manifest = _read_manifest(path)
    except FileNotFoundError as exc:
        raise BackupError(f"There is no backup named {name!r}.") from exc
    except (OSError, tarfile.TarError, EOFError, ValueError) as exc:
        return Verification(False, [f"{name} can't be read ({exc})"])
    expected = {}
    for part in manifest["parts"]:
        if not isinstance(part, dict) or not _PART.match(str(part.get("name", ""))):
            problems.append("the manifest lists a part with no usable name")
            continue
        for item in part.get("files") or []:
            rel = str(item.get("path", ""))
            if not _FILE.match(rel):
                problems.append(f"{part['name']}: {rel!r} is not a usable file name")
                continue
            expected[f"{part['name']}/{rel}"] = item
    seen = set()
    try:
        with tarfile.open(path, "r:bz2") as tar:
            for member in tar:
                if member.name == MANIFEST:
                    continue
                if not _safe_member(member.name) or not member.isfile():
                    problems.append(f"{member.name!r} is not a plain file inside the archive")
                    continue
                item = expected.get(member.name)
                if item is None:
                    problems.append(f"{member.name} is in the archive but not in its manifest")
                    continue
                digest, size = hashlib.sha256(), 0
                handle = tar.extractfile(member)
                for block in iter(lambda: handle.read(1 << 20), b""):
                    digest.update(block)
                    size += len(block)
                if size != item.get("size") or digest.hexdigest() != item.get("sha256"):
                    problems.append(f"{member.name} does not match its recorded hash; the backup is damaged")
                seen.add(member.name)
    except (OSError, tarfile.TarError, EOFError) as exc:
        problems.append(f"{name} can't be read to the end ({exc})")
    problems.extend(f"{missing} is in the manifest but not in the archive" for missing in sorted(set(expected) - seen))
    return Verification(not problems, problems, manifest)


def read_file(name, part, filename, directory=None):
    """The bytes of one file of one part, checked against the manifest's size and hash (for a page's diff).

    Raises:
        BackupError: No such archive, part or file, or the file doesn't match the manifest.
    """
    folder = directory or backups_dir()
    path = os.path.join(folder, _check_name(name))
    if not _PART.match(part or "") or not _FILE.match(filename or ""):
        raise BackupError(f"There is no file {part}/{filename} in {name}.")
    try:
        manifest = _read_manifest(path)
        item = next((f for p in manifest["parts"] if p.get("name") == part for f in p.get("files") or []
                     if f.get("path") == filename), None)
        if item is None:
            raise BackupError(f"{name} has no file {part}/{filename}.")
        with tarfile.open(path, "r:bz2") as tar:
            member = tar.getmember(f"{part}/{filename}")
            if not member.isfile() or member.size > MAX_FILE_BYTES:
                raise BackupError(f"{name} has no usable file {part}/{filename}.")
            data = tar.extractfile(member).read()
    except FileNotFoundError as exc:
        raise BackupError(f"There is no backup named {name!r}.") from exc
    except KeyError as exc:
        raise BackupError(f"{name} has no file {part}/{filename}.") from exc
    except (OSError, tarfile.TarError, EOFError, ValueError) as exc:
        if isinstance(exc, BackupError):
            raise
        raise BackupError(f"{name} can't be read ({exc}).") from exc
    if len(data) != item.get("size") or hashlib.sha256(data).hexdigest() != item.get("sha256"):
        raise BackupError(f"{name}: {part}/{filename} does not match its recorded hash; the backup is damaged.")
    return data


def _prune(folder, keep=None):
    if keep is None:
        keep = kept()
    rows = list_archives(folder)
    protected = next((a.name for a in rows if a.complete and a.readable), None)  # the newest full archive stays
    gone = []
    for old in rows[max(int(keep), 1):]:
        if old.name == protected:
            continue
        try:
            os.unlink(old.path)
            gone.append(old.name)
        except OSError:
            pass
    return gone


def prune(directory=None, keep=None):
    """Deletes all but the newest `keep` archives (default `settings.backups_kept`) and never the newest one that
    holds every registered part; returns the names deleted."""
    folder = directory or backups_dir()
    with _Lock(folder):
        return _prune(folder, keep)


def kept():
    """`settings.backups_kept`, or the default when the settings can't be read."""
    from planetgen.util import settings
    try:
        return settings.get_settings().settings.backups_kept
    except settings.SettingsError:
        return DEFAULT_KEPT


def restore(name, parts=None, context=None, directory=None):
    """
    Restores the whole archive `name`, or only `parts`, after verifying it.

    Each part to restore is first backed up in a safety archive (`RestoreReport.safety`), then its provider's
    `restore` is given a hash-checked copy of the part's files. A part the archive holds but this build has no
    provider for is skipped and listed in `skipped`; naming a part the archive doesn't hold is an error.

    Raises:
        BackupError: The archive is missing or fails `verify`, or a named part isn't in it. Nothing is changed.
    """
    context = context or RestoreContext()
    folder = directory or backups_dir()
    checked = verify(name, folder)
    if not checked.ok:
        raise BackupError(f"{name} is not a sound backup and was not restored: " + "; ".join(checked.problems))
    held = {p["name"]: p for p in checked.manifest["parts"] if isinstance(p, dict)}
    wanted = list(held) if parts is None else list(dict.fromkeys(parts))
    missing = [p for p in wanted if p not in held]
    if missing:
        raise BackupError(f"{name} holds no " + ", ".join(repr(p) for p in missing) + " part.")
    report = RestoreReport(archive=name)
    report.skipped = [p for p in wanted if p not in _PARTS]
    todo = [p for p in wanted if p in _PARTS]
    if not todo:
        return report
    report.safety = take(todo, who=context.who, reason=f"before restoring {name}", directory=folder).name
    path = os.path.join(folder, name)
    work = tempfile.mkdtemp(prefix=".restore-", dir=folder)
    try:
        with tarfile.open(path, "r:bz2") as tar:
            for part in todo:
                part_dir = os.path.join(work, part)
                os.mkdir(part_dir)
                for item in held[part].get("files") or []:
                    target = os.path.join(part_dir, *item["path"].split("/"))
                    os.makedirs(os.path.dirname(target), exist_ok=True)
                    with open(target, "wb") as out:
                        shutil.copyfileobj(tar.extractfile(f"{part}/{item['path']}"), out)
        for part in todo:
            report.results[part] = _PARTS[part].restore(os.path.join(work, part), context)
            report.restored.append(part)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    activity_log.event("DB", "backup.restore", user=context.who, archive=name, parts=",".join(report.restored),
                       skipped=",".join(report.skipped))
    if context.audit is not None:
        try:
            context.audit("backup.restore", f"backup:{name}", json.dumps({"parts": report.restored,
                                                                          "skipped": report.skipped}))
        except Exception as exc:  # noqa: BLE001 -- the restore is done; say so, don't undo it
            activity_log.event("DB", "backup.audit_failed", user=context.who, error=str(exc)[:200])
    return report
