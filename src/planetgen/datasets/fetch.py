# planetgen/datasets/fetch.py

"""
The shared setup-time download step (GEN.548): fetch a dataset's files once, check them, keep them outside the
checkout, and hand a reader either the files or a labelled stand-in.

Boss's rule (2026-10-18): a static file is downloaded at setup, never at run time, unless the data needs constant
updating. So the two halves are separate. `fetch` is the setup half: `python -m planetgen.cli.datasets fetch`, run by
`update.sh` and `install.sh`, sends the standard User-Agent (`planetgen.util.user_agent`, API.31), reads at most the
size the catalogue states, compares the sha256 and puts the file in the cache folder by an atomic rename. Any failure
(host down, timeout, HTTP error, too many bytes, wrong hash) keeps what was there before, gives a reason and leaves the
dataset on its stand-in; it never raises, so a host that is down cannot stop an update. `locate` is the run-time half:
it opens no connection, only looks in the folder (hashing each file once per process), and answers with the folder or
with the stand-in label. A reader that gets a stand-in must say so in what it returns (the stand-in registry,
docs/design/gap-filling-plan.md section 6).

Tests never reach the network: `fetch` takes the function that opens a request (`opener`), and a test passes one that
serves bytes it built itself, with the sha256 computed from them.
"""

import hashlib
import http.client
import os
import re
import tempfile
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
"""int: How long one request may wait for the host before the dataset falls back to its stand-in."""

CHUNK_BYTES = 1 << 16

_SHA256 = re.compile(r"[0-9a-f]{64}")
_FILE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")


@dataclass(frozen=True)
class Piece:
    """One file of a dataset: its name in the cache folder, where it comes from, its sha256 and its exact size."""

    name: str
    url: str
    sha256: str
    size: int

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

    `key` names the cache sub-folder. `source` is the citation and `licence` the licence of the data (a record only).
    `standin` says in a few words what the readers use when the files are missing, and `fill` is the stand-in
    registry id (gap plan section 6) of that fill.
    """

    key: str
    title: str
    source: str
    licence: str
    pieces: tuple
    standin: str
    fill: str

    def __post_init__(self):
        if not _FILE_NAME.fullmatch(self.key):
            raise ValueError(f"a dataset key is a plain folder name, not {self.key!r}")
        if not self.pieces:
            raise ValueError(f"{self.key}: a dataset has at least one file")
        names = [piece.name for piece in self.pieces]
        if len(set(names)) != len(names):
            raise ValueError(f"{self.key}: two files share a name")


@dataclass(frozen=True)
class Result:
    """What `fetch` did: `status` is `fetched` (something was downloaded), `present` (everything was already there)
    or `stand-in` (a file could not be had; `reasons` say why for each)."""

    dataset: Dataset
    status: str
    reasons: tuple = ()


@dataclass(frozen=True)
class Located:
    """What a reader gets: the folder holding every verified file, or `folder` None and the reason for the stand-in."""

    dataset: Dataset
    folder: str
    reason: str = ""

    @property
    def standin(self):
        """True when the files are missing and the reader must use the stand-in."""
        return self.folder is None

    @property
    def label(self):
        """A few words saying where the numbers come from, for a reader to pass on."""
        if self.standin:
            return f"stand-in ({self.dataset.fill}): {self.dataset.standin}"
        return f"downloaded: {self.dataset.source}"

    def path(self, name):
        """The path of the file `name`; not available on a stand-in."""
        if self.standin:
            raise FileNotFoundError(f"{self.dataset.key} is on its stand-in ({self.reason}); it has no file {name}")
        if name not in {piece.name for piece in self.dataset.pieces}:
            raise KeyError(f"{self.dataset.key} has no file {name!r}")
        return os.path.join(self.folder, name)


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


def locate(dataset):
    """Where the dataset's files are, or that it is on its stand-in. Opens no connection (setup downloads, run time
    only reads): a missing, short or altered file is a stand-in with the reason, never an error."""
    where = folder(dataset)
    for piece in dataset.pieces:
        path = os.path.join(where, piece.name)
        if not os.path.exists(path):
            return Located(dataset, None, f"{piece.name} has not been downloaded")
        if not _is_good(path, piece):
            log.debug("dataset %s: %s does not match the catalogue (size or sha256)", dataset.key, piece.name)
            return Located(dataset, None, f"{piece.name} does not match the catalogue's size or sha256")
    return Located(dataset, where)


class DownloadError(Exception):
    """One file could not be had; the message says why."""


def _download(piece, path, opener):
    request = urllib.request.Request(piece.url, headers={"User-Agent": user_agent()})
    limit = piece.size + 1
    data = bytearray()
    try:
        with opener(request, timeout=TIMEOUT_SECONDS) as response:
            while len(data) < limit:
                block = response.read(min(CHUNK_BYTES, limit - len(data)))
                if not block:
                    break
                data.extend(block)
    except urllib.error.HTTPError as exc:
        raise DownloadError(f"{piece.url}: the host answered {exc.code} {exc.reason}") from exc
    except (urllib.error.URLError, OSError, http.client.HTTPException) as exc:
        raise DownloadError(f"{piece.url}: {getattr(exc, 'reason', None) or exc}") from exc
    if len(data) != piece.size:
        raise DownloadError(f"{piece.url}: expected {piece.size} bytes, got {'more' if len(data) > piece.size else len(data)}")
    found = hashlib.sha256(data).hexdigest()
    if found != piece.sha256:
        raise DownloadError(f"{piece.url}: sha256 is {found}, the catalogue says {piece.sha256}; not kept")
    handle, temporary = tempfile.mkstemp(prefix=f".{piece.name}.", dir=os.path.dirname(path))
    try:
        with os.fdopen(handle, "wb") as out:
            out.write(data)
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    except OSError as exc:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise DownloadError(f"{path}: cannot be written ({exc})") from exc


def fetch(dataset, *, opener=None):
    """Brings the dataset's files into the cache folder. A file already there with the right size and sha256 is left
    alone; a different one is replaced. Never raises for a network or disk problem: the result is `stand-in` with a
    reason per file that could not be had, and the files that were fine stay."""
    opener = opener or urllib.request.urlopen
    where = folder(dataset)
    try:
        os.makedirs(where, mode=0o755, exist_ok=True)
    except OSError as exc:
        return Result(dataset, "stand-in", (f"{where}: cannot be made ({exc})",))
    downloaded, reasons = 0, []
    for piece in dataset.pieces:
        path = os.path.join(where, piece.name)
        if _is_good(path, piece):
            continue
        try:
            _download(piece, path, opener)
            downloaded += 1
        except DownloadError as exc:
            reasons.append(str(exc))
    if reasons:
        return Result(dataset, "stand-in", tuple(reasons))
    return Result(dataset, "fetched" if downloaded else "present")
