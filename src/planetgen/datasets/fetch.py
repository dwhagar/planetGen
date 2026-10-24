# planetgen/datasets/fetch.py

"""
The shared setup-time download step (GEN.548): fetch a dataset's files once, check them, keep them outside the
checkout, and hand a reader the files or say that they are missing.

Boss's rule (2026-10-18): a static file is downloaded at setup, never at run time, unless the data needs constant
updating. So the two halves are separate. `fetch` is the setup half: `python -m planetgen.cli.datasets fetch`, run by
`update.sh` and `install.sh`, sends the standard User-Agent (`planetgen.util.user_agent`, API.31), reads at most the
size the catalogue states, compares the sha256 and puts the file in the cache folder by an atomic rename. Any failure
(host down, timeout, HTTP error, too many bytes, wrong hash) keeps what was there before, gives a reason and leaves the
dataset missing; it never raises, so a host that is down cannot stop an update. `locate` is the run-time half:
it opens no connection, only looks in the folder (hashing each file once per process), and answers with the folder or
that the dataset is missing. The generator requires every dataset (GEN.590, Boss's answer on the card): a reader raises
`DataMissing` and bulk generation refuses to start (`run_galaxy.require_datasets`), so no value comes from a fallback unless
the run was told, in so many words, to disregard missing data (GEN.599, `set_disregard_missing_data`: `--disregard-missing-data`
on the command line, a checkbox on the Generate page). Then a reader asks `Located.use_stand_in` and takes the dataset's own
labelled `fallback`, and `stand_ins_used` lists what was taken. What each fill replaced is in the stand-in registry
(docs/design/research/gap-filling-plan.md section 6).

Some sources are large archives of which only a few columns are used (the MIST isochrones are 220 MB). A dataset can
then name a `derive` function and the compact files it writes (`derived`); its archive piece is `transient`: the step
downloads and checks it, calls `derive(folder)` once, records the sha256 of every piece and of every derived file in
`.derived.json`, and deletes the archive. `locate` trusts the derived files only while that record matches the
catalogue and the files, so a new catalogue sha256 or `derive_version` makes the next setup download and derive again.

Tests never reach the network: `fetch` takes the function that opens a request (`opener`), and a test passes one that
serves bytes it built itself, with the sha256 computed from them.
"""

import contextlib
import hashlib
import http.client
import json
import os
import re
import tempfile
import time
import urllib.error
import urllib.request
from dataclasses import dataclass

from planetgen.util import log
from planetgen.util.user_agent import user_agent

CACHE_DIR_ENV = "PLANETGEN_DATA_DIR"
"""str: The environment variable that moves the cache folder."""

DEFAULT_CACHE_DIR = "/var/cache/planetgen/datasets"
"""str: Where the downloaded files are kept (beside the tile cache, outside the checkout); `update.sh` fills it as root
and makes it world-readable, since the web app only reads."""

TIMEOUT_SECONDS = 30
"""int: How long one request may wait for the host before the dataset is given up as missing."""

ATTEMPTS = 4
"""int: How many times a file is requested when the host does not answer or stalls (a wrong size or hash is never retried). A
retry asks for the rest of the file (a `Range` request) and starts over only when the host does not honour it."""

STALL_WINDOW_SECONDS = 60
"""int: A download that brings fewer than `STALL_MIN_BYTES` in this long is given up as stalled (the socket timeout does not catch a
host that trickles a few bytes a minute) and retried from where it stopped."""

STALL_MIN_BYTES = 1 << 17
"""int: The least a download must bring in `STALL_WINDOW_SECONDS` (about 2 kB/s)."""

PROGRESS_SECONDS = 10
"""int: How often the plain progress (`PlainProgress`, output that is not a terminal) says how far a download has come, so a slow one is
not mistaken for a hang."""

RETRY_WAIT_SECONDS = 3
"""int: The wait between two attempts."""

STAMP = ".derived.json"
"""str: The record `fetch` keeps in a dataset's folder when it derives files."""

CHUNK_BYTES = 1 << 16

PROGRESS_MIN_BYTES = 1 << 22
"""int: A file at least this big gets a line when its download starts (and one every `PROGRESS_SECONDS`) from `PlainProgress`."""

_SHA256 = re.compile(r"[0-9a-f]{64}")
_FILE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")


@dataclass(frozen=True)
class Piece:
    """One file of a dataset: its name in the cache folder, where it comes from, its sha256 and its exact size."""

    name: str
    url: str
    sha256: str
    size: int
    transient: bool = False

    def __post_init__(self):
        if not _FILE_NAME.fullmatch(self.name) or self.name.startswith(".."):
            raise ValueError(f"a dataset file name is a plain file name, not {self.name!r}")
        if not self.url.startswith("https://"):
            raise ValueError(f"{self.name}: a dataset is fetched over https, not {self.url!r}")
        if not _SHA256.fullmatch(self.sha256):
            raise ValueError(f"{self.name}: sha256 must be 64 lower-case hex digits")
        if self.size <= 0:
            raise ValueError(f"{self.name}: the size in bytes must be positive")


@dataclass(frozen=True)
class Dataset:
    """A set of files from one published source and what stands in for it until they are present.

    `key` names the cache sub-folder. `source` is the citation.
    `replaces` says in a few words what the data replaced (the stand-in of the gap plan's registry, named in the message
    when the files are missing), and `fill` is the registry id (gap plan section 6) of that fill. A dataset that keeps only compact files derived from a large
    archive names them in `derived`, the function that writes them in `derive` (`derive(folder)` returns
    `{file name: bytes}` from the verified pieces in `folder`) and a `derive_version` to raise when that function
    changes; the archive pieces are marked `transient` and are removed once derived. `fallback` says in a few words what a
    reader answers with when the files are missing and the run disregards missing data (GEN.599); empty when there is
    none, and then a missing dataset is refused whatever the run was told.
    """

    key: str
    title: str
    source: str
    pieces: tuple
    replaces: str
    fill: str
    derived: tuple = ()
    derive: object = None
    derive_version: int = 1
    fallback: str = ""

    def __post_init__(self):
        if not _FILE_NAME.fullmatch(self.key):
            raise ValueError(f"a dataset key is a plain folder name, not {self.key!r}")
        if not self.pieces:
            raise ValueError(f"{self.key}: a dataset has at least one file")
        names = [piece.name for piece in self.pieces]
        names += list(self.derived)
        if len(set(names)) != len(names):
            raise ValueError(f"{self.key}: two files share a name")
        if bool(self.derived) != (self.derive is not None):
            raise ValueError(f"{self.key}: `derived` and `derive` go together")
        if any(piece.transient for piece in self.pieces) and not self.derive:
            raise ValueError(f"{self.key}: a transient file is removed after deriving, so a derive function is needed")


@dataclass(frozen=True)
class PieceResult:
    """What `fetch` did with one file (GEN.600): `state` is `fetched` (requested and kept), `kept` (already there with
    the right size and sha256, or an archive whose derived files are current: not requested) or `failed` (requested or
    not writable; `reason` says why)."""

    piece: Piece
    state: str
    reason: str = ""


@dataclass(frozen=True)
class Result:
    """What `fetch` did: `status` is `fetched` (something was downloaded), `present` (everything was already there)
    or `missing` (a file could not be had; `reasons` say why for each). `pieces` is one `PieceResult` per file the
    dataset lists, in order (empty when the folder could not be made)."""

    dataset: Dataset
    status: str
    reasons: tuple = ()
    pieces: tuple = ()


@dataclass(frozen=True)
class Located:
    """What a reader gets: the folder holding every verified file, or `folder` None and the reason it is missing.
    `lacking` is, for a missing dataset, one `(Piece, reason)` for each file that has to be downloaded again (GEN.600),
    so a message can name every URL."""

    dataset: Dataset
    folder: str
    reason: str = ""
    lacking: tuple = ()

    @property
    def missing(self):
        """True when the files are missing and the reader cannot answer."""
        return self.folder is None

    @property
    def label(self):
        """A few words saying where the numbers come from, for a reader to pass on."""
        if self.missing:
            return f"missing ({self.dataset.fill}, replaces {self.dataset.replaces})"
        return f"downloaded: {self.dataset.source}"

    def describe(self):
        """The refusal text for a missing dataset: its name, why, every file to download with its URL and the setup
        command (GEN.599)."""
        urls = ", ".join(piece.url for piece, _reason in self.lacking) or "none"
        return (f"{self.dataset.key} ({self.dataset.fill}): {self.reason}; files to download: {urls}; "
                f"run 'python3 -m planetgen.cli.datasets fetch {self.dataset.key}' (update.sh does it at setup)")

    def require(self):
        """Self, or `DataMissing` when the files are missing, for a reader that has no value of its own to give."""
        if self.missing:
            raise DataMissing(f"{self.label}: {self.describe()}")
        return self

    def use_stand_in(self, what):
        """
        Whether a reader answers from the dataset's labelled fallback instead of the files: False when the files are there,
        True when they are missing and the run disregards missing data (`what` names what was taken, for the run's record),
        and `DataMissing` when they are missing and it does not (or the dataset has no fallback).
        """
        if not self.missing:
            return False
        if not _DISREGARD or not self.dataset.fallback:
            self.require()
        _STAND_INS.add((self.dataset.key, what))
        return True

    def path(self, name):
        """The path of the file `name`; not available while the files are missing."""
        if self.missing:
            raise FileNotFoundError(f"{self.dataset.key} is missing ({self.reason}); it has no file {name}")
        if name not in {piece.name for piece in self.dataset.pieces if not piece.transient} | set(self.dataset.derived):
            raise KeyError(f"{self.dataset.key} has no file {name!r}")
        return os.path.join(self.folder, name)


_DISREGARD = False
_STAND_INS = set()


def set_disregard_missing_data(on):
    """Tells the process whether a reader may answer from a dataset's fallback while its files are missing (GEN.599). Only
    a run that was told to (`--disregard-missing-data`) turns it on; the default is to refuse."""
    global _DISREGARD
    _DISREGARD = bool(on)


def disregarding_missing_data():
    """Whether this process was told to disregard missing data."""
    return _DISREGARD


@contextlib.contextmanager
def disregard_missing_data():
    """Disregards missing data inside the block and forgets what was used when it ends (for a test or a one-off call)."""
    before = _DISREGARD
    set_disregard_missing_data(True)
    try:
        yield
    finally:
        set_disregard_missing_data(before)
        _STAND_INS.clear()


def stand_ins_used():
    """The fallbacks taken in this process so far, as sorted `(dataset key, what)` pairs."""
    return sorted(_STAND_INS)


def cache_dir():
    """The folder the datasets are kept in: `PLANETGEN_DATA_DIR`, else `/var/cache/planetgen/datasets`."""
    return os.environ.get(CACHE_DIR_ENV) or DEFAULT_CACHE_DIR


def folder(dataset):
    """The dataset's own folder inside the cache."""
    return os.path.join(cache_dir(), dataset.key)


def _sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


_VERIFIED = {}
"""dict: `(path, size, mtime_ns, sha256)` of every file already hashed in this process, so a reader that asks for a
dataset again does not read its files again."""


def _is_good(path, piece):
    """True when `path` is the file `piece` describes: right size, right sha256."""
    try:
        stat = os.stat(path)
    except OSError:
        return False
    if stat.st_size != piece.size:
        return False
    key = (path, stat.st_size, stat.st_mtime_ns, piece.sha256)
    if key in _VERIFIED:
        return True
    try:
        good = _sha256_of(path) == piece.sha256
    except OSError:
        return False
    if good:
        _VERIFIED[key] = True
    return good


def _stamp_path(dataset):
    return os.path.join(folder(dataset), STAMP)


def _expected_stamp(dataset):
    return {"version": dataset.derive_version, "pieces": {piece.name: piece.sha256 for piece in dataset.pieces}}


def _derived_ok(dataset):
    """True when the dataset's derived files are all there and match the record made when they were derived, and the
    record was made from the catalogue's current pieces and derive version."""
    try:
        with open(_stamp_path(dataset), encoding="utf-8") as handle:
            stamp = json.load(handle)
    except (OSError, ValueError):
        return False
    expected = _expected_stamp(dataset)
    if not isinstance(stamp, dict) or any(stamp.get(key) != value for key, value in expected.items()):
        return False
    files = stamp.get("files")
    if not isinstance(files, dict) or set(files) != set(dataset.derived):
        return False
    where = folder(dataset)
    return all(_is_good(os.path.join(where, name), _derived_piece(name, files[name]["sha256"], files[name]["size"]))
               for name in dataset.derived)


def _derived_piece(name, sha256, size):
    """A stand-in `Piece` for a derived file, so the same hash check serves both."""
    return Piece(name, "https://derived.invalid/" + name, sha256, size)


def locate(dataset):
    """Where the dataset's files are, or that they are missing. Opens no connection (setup downloads, run time
    only reads): a missing, short or altered file gives a `Located` with the reason and every file to download again,
    never an error."""
    where = folder(dataset)
    lacking = []
    for piece in dataset.pieces:
        if piece.transient:
            continue
        path = os.path.join(where, piece.name)
        if not os.path.exists(path):
            lacking.append((piece, f"{piece.name} has not been downloaded"))
        elif not _is_good(path, piece):
            log.debug("dataset %s: %s does not match the catalogue (size or sha256)", dataset.key, piece.name)
            lacking.append((piece, f"{piece.name} does not match the catalogue's size or sha256"))
    if lacking:
        return Located(dataset, None, lacking[0][1], tuple(lacking))
    if dataset.derived and not _derived_ok(dataset):
        # the derived files are made again from the archive, so the archive is what a setup run downloads
        reason = "its derived files have not been made from the current catalogue"
        archive = [piece for piece in dataset.pieces if piece.transient] or list(dataset.pieces)
        return Located(dataset, None, reason, tuple((piece, reason) for piece in archive))
    return Located(dataset, where)


class DataMissing(LookupError):
    """A reader asked for a dataset whose files were not downloaded; the message is the label and the reason."""


class DownloadError(Exception):
    """One file could not be had; the message says why."""


def _now():
    return time.monotonic()


class PlainProgress:
    """
    The progress of a download as plain lines, for output that is not a terminal (a log file, a pipe): a line when a file of at least
    `PROGRESS_MIN_BYTES` starts and one every `PROGRESS_SECONDS` with the bytes, the speed and the time left. `fetch(progress=)` takes this or
    any object with the same four methods (the update command's bars, `planetgen.cli.datasets`).
    """

    def __init__(self, write):
        self._write = write
        self._began = {}
        self._said = {}

    def begin(self, piece):
        if piece.size >= PROGRESS_MIN_BYTES:
            self._write(f"    downloading {piece.name} ({_megabytes(piece.size)}) from {piece.url}")
        self._began[piece.name] = self._said[piece.name] = _now()

    def update(self, piece, count):
        now = _now()
        if piece.size < PROGRESS_MIN_BYTES or now - self._said.get(piece.name, now) < PROGRESS_SECONDS:
            return
        self._said[piece.name] = now
        speed = count / max(now - self._began[piece.name], 1e-9)
        left = f", about {int((piece.size - count) / speed)} s left" if speed > 0 and count < piece.size else ""
        self._write(f"    {piece.name}: {_megabytes(count)} of {_megabytes(piece.size)} ({100 * count // max(piece.size, 1)}%), "
                    f"{speed / 1e6:.2f} MB/s{left}")

    def note(self, text):
        self._write(f"    {text}")

    def end(self, piece):
        self._began.pop(piece.name, None)
        self._said.pop(piece.name, None)


def _megabytes(count):
    return f"{count / 1e6:.1f} MB"


def _request(piece, opener, out, progress=None):
    """Streams the piece into the open file `out` from where it ends (a `Range` request when it already holds a start); returns the
    byte count the file holds. Raises `DownloadError`, also when the host stalls (`STALL_MIN_BYTES` in `STALL_WINDOW_SECONDS`)."""
    out.seek(0, os.SEEK_END)
    count = out.tell()
    headers = {"User-Agent": user_agent()}
    if count:
        headers["Range"] = f"bytes={count}-"
    request = urllib.request.Request(piece.url, headers=headers)
    limit = piece.size + 1
    try:
        with opener(request, timeout=TIMEOUT_SECONDS) as response:
            if count and getattr(response, "status", 200) != 206:
                out.seek(0)
                out.truncate()
                count = 0
            window_start = _now()
            window_count = count
            while count < limit:
                block = response.read(min(CHUNK_BYTES, limit - count))
                if not block:
                    break
                out.write(block)
                count += len(block)
                now = _now()
                if now - window_start >= STALL_WINDOW_SECONDS:
                    if count - window_count < STALL_MIN_BYTES:
                        raise DownloadError(f"{piece.url}: stalled at {_megabytes(count)} of {_megabytes(piece.size)}")
                    window_start, window_count = now, count
                if progress:
                    progress.update(piece, count)
    except urllib.error.HTTPError as exc:
        raise DownloadError(f"{piece.url}: the host answered {exc.code} {exc.reason}") from exc
    except (urllib.error.URLError, OSError, http.client.HTTPException) as exc:
        raise DownloadError(f"{piece.url}: {getattr(exc, 'reason', None) or exc}") from exc
    return count


class _WrongBytes(DownloadError):
    """The host answered, but with the wrong file: never retried."""


def _sha256_of_open(out):
    out.seek(0)
    digest = hashlib.sha256()
    for block in iter(lambda: out.read(CHUNK_BYTES), b""):
        digest.update(block)
    return digest.hexdigest()


def _download(piece, path, opener, progress=None):
    handle, temporary = tempfile.mkstemp(prefix=f".{piece.name}.", dir=os.path.dirname(path))
    try:
        with os.fdopen(handle, "w+b") as out:
            if progress:
                progress.begin(piece)
            for attempt in range(ATTEMPTS):
                try:
                    count = _request(piece, opener, out, progress)
                    break
                except DownloadError as exc:
                    if attempt == ATTEMPTS - 1:
                        raise
                    if progress:
                        progress.note(f"{exc}; trying again ({attempt + 2} of {ATTEMPTS})")
                    time.sleep(RETRY_WAIT_SECONDS)
            if count != piece.size:
                raise _WrongBytes(f"{piece.url}: expected {piece.size} bytes, got {'more' if count > piece.size else count}")
            found = _sha256_of_open(out)
        if found != piece.sha256:
            raise _WrongBytes(f"{piece.url}: sha256 is {found}, the catalogue says {piece.sha256}; not kept")
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    except OSError as exc:
        _unlink(temporary)
        raise DownloadError(f"{path}: cannot be written ({exc})") from exc
    except BaseException:
        _unlink(temporary)
        raise
    finally:
        if progress:
            progress.end(piece)


def _unlink(path):
    try:
        os.unlink(path)
    except OSError:
        pass


def _write_atomically(path, data):
    handle, temporary = tempfile.mkstemp(prefix=f".{os.path.basename(path)}.", dir=os.path.dirname(path))
    try:
        with os.fdopen(handle, "wb") as out:
            out.write(data)
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    except BaseException:
        _unlink(temporary)
        raise


def _derive(dataset):
    """Runs the dataset's derive function over its verified pieces and writes the files and the record."""
    where = folder(dataset)
    made = dataset.derive(where)
    if set(made) != set(dataset.derived):
        raise ValueError(f"derive made {sorted(made)}, the catalogue says {sorted(dataset.derived)}")
    files = {}
    for name in dataset.derived:
        _write_atomically(os.path.join(where, name), made[name])
        files[name] = {"sha256": hashlib.sha256(made[name]).hexdigest(), "size": len(made[name])}
    stamp = dict(_expected_stamp(dataset), files=files)
    _write_atomically(_stamp_path(dataset), json.dumps(stamp, sort_keys=True).encode("utf-8"))
    for piece in dataset.pieces:
        if piece.transient:
            _unlink(os.path.join(where, piece.name))


def fetch(dataset, *, opener=None, progress=None):
    """Brings the dataset's files into the cache folder. A file already there with the right size and sha256 is left
    alone; a different one is replaced. Never raises for a network or disk problem: the result is `missing` with a
    reason per file that could not be had, and the files that were fine stay. A dataset with a derive function whose
    derived files are current does not download its archive again. `Result.pieces` says what happened to each file. `progress`, when
    given (`PlainProgress`, or an object with its `begin`, `update`, `note` and `end`), is told how each download goes."""
    opener = opener or urllib.request.urlopen
    where = folder(dataset)
    try:
        os.makedirs(where, mode=0o755, exist_ok=True)
    except OSError as exc:
        return Result(dataset, "missing", (f"{where}: cannot be made ({exc})",))
    needs_derive = bool(dataset.derived) and not _derived_ok(dataset)
    changed, reasons, handled = False, [], []
    for piece in dataset.pieces:
        if piece.transient and not needs_derive:
            handled.append(PieceResult(piece, "kept"))
            continue
        path = os.path.join(where, piece.name)
        if _is_good(path, piece):
            handled.append(PieceResult(piece, "kept"))
            continue
        try:
            _download(piece, path, opener, progress)
            changed = True
            handled.append(PieceResult(piece, "fetched"))
        except DownloadError as exc:
            reasons.append(str(exc))
            handled.append(PieceResult(piece, "failed", str(exc)))
    if reasons:
        return Result(dataset, "missing", tuple(reasons), tuple(handled))
    if needs_derive:
        try:
            _derive(dataset)
            changed = True
        except Exception as exc:  # a parser meeting an unexpected file must not stop an update
            return Result(dataset, "missing", (f"{dataset.key}: deriving its files failed ({type(exc).__name__}: {exc})",),
                          tuple(handled))
    return Result(dataset, "fetched" if changed else "present", (), tuple(handled))
