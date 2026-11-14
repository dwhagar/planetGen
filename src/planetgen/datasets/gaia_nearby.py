# planetgen/datasets/gaia_nearby.py

"""
The Gaia DR3 stars within 25 parsecs (GEN.1534): 5,323 sources with a parallax over 40 milliarcseconds measured to better than 10 sigma, each with its position, parallax, proper motion,
G magnitude and BP-RP colour. Gaia is complete to the faintest M dwarfs at this distance, which the HYG catalogue (`hyg_stars`: incomplete, and mixing wide pairs) is not, so it is the
sample the spacing of stars is read from: how near a system's nearest neighbour is (`nearest_neighbour_distances_pc`). The galaxy's random placement keeps two systems apart by a share of
their galactic Hill radii (`tuning.SECTOR_HILL_SEPARATION_SHARE`, stand-in SEP1), and this is the statistic that share is fitted to (`scripts/fit_sector_separation.py`,
`tests/test_sector_separation_fit.py`). Nothing in the generator reads the data: it is a fit target and a test oracle, so its dataset is `required` False, as HYG's is.

The file is the result of one query of the Gaia archive's TAP service (`QUERY`), sorted by `source_id` so that the same query returns the same bytes (Gaia DR3 does not change). The data are
Gaia Collaboration, DR3, released under CC BY-SA 3.0 IGO (ESA/Gaia/DPAC); the file is never put in the repository, it goes to the cache folder like the other datasets.

The reader raises `fetch.DataMissing` when the file has not been downloaded.
"""

import csv
import functools
import http.cookiejar
import io
import math
import os
import time
import urllib.error
import urllib.parse
import urllib.request

from planetgen.datasets import fetch
from planetgen.util import secrets_mask
from planetgen.util.user_agent import user_agent

FILE = "gaia_dr3_25pc.csv"

QUERY = ("SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source "
         "WHERE parallax > 40 AND parallax_over_error > 10 ORDER BY source_id")
"""str: The ADQL query that makes the file: every DR3 source nearer than 25 pc (parallax over 40 mas) whose parallax is measured to better than 10 sigma."""

URL = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": QUERY})

ESA_TAP = "https://gea.esac.esa.int/tap-server/tap"
"""str: The Gaia archive's TAP service (ESA). Its anonymous synchronous query stops at 2,000 rows and 429s when the queue is full (OPS.183), which is why the sources below exist."""

ARI_TAP = "https://gaia.ari.uni-heidelberg.de/tap/sync"
"""str: The Heidelberg ARI mirror of the Gaia archive; the same ADQL gives the same file byte for byte, with the row limit raised by `MAXREC`."""

VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
"""str: VizieR's TAP service (CDS), which holds Gaia DR3 as `I/355/gaiadr3` and the Gaia Catalogue of Nearby Stars as `I/360/gcnsmain`."""

ESA_USER_ENV, ESA_PASSWORD_ENV = "PLANETGEN_ESA_USER", "PLANETGEN_ESA_PASSWORD"
"""str: Environment variables for an optional registered Gaia archive account (OPS.183). With both set the asynchronous query logs in first, which lifts the anonymous queue limit.
The password is registered with the log masking (SEC.48), so it can never reach a log."""

ASYNC_WAIT_SECONDS = 600
"""int: How long the asynchronous query may run before it is given up."""

ASYNC_POLL_SECONDS = 3

MAX_ROWS = 5700
MIN_ROWS = 5000
"""int: The rows a file of this sample holds, between these (DR3 gives 5,323; another catalogue's cuts may differ a little, which is why a source that is not byte for byte the same must stay inside)."""

COLUMNS = ("source_id", "ra", "dec", "parallax", "parallax_error", "pmra", "pmdec", "phot_g_mean_mag", "bp_rp")

LIMIT_PC = 25.0
"""float: The distance the sample reaches, parsecs (parallax 40 mas)."""

MERGE_PC = 0.1
"""float: Stars closer than this to a neighbour are one system (wide pairs and their companions, whose separations run to 20,000 AU), as the galaxy places a system once however many stars it holds."""

CENTER_PC = 15.0
"""float: Only systems nearer than this are asked for their nearest neighbour, so that the neighbour is always inside the complete sample (`LIMIT_PC`) and none is missed at its edge: the nearest neighbour of a system nearer than `CENTER_PC` is nearer than 10 pc."""

def _post(opener, url, fields, timeout, cookies=None):
    """One request to `url` (a POST with `fields`, or a GET when there are none) through `opener`; returns `(bytes, final url)`."""
    data = urllib.parse.urlencode(fields).encode("ascii") if fields is not None else None
    request = urllib.request.Request(url, data=data, headers={"User-Agent": user_agent()})
    try:
        with opener(request, timeout=timeout) as response:
            return response.read(), (response.geturl() if hasattr(response, "geturl") else url)
    except urllib.error.HTTPError as exc:
        raise fetch.DownloadError(f"{url}: the host answered {exc.code} {exc.reason}") from exc
    except (urllib.error.URLError, OSError) as exc:
        raise fetch.DownloadError(f"{url}: {getattr(exc, 'reason', None) or exc}") from exc


def esa_async(opener, timeout, wait=ASYNC_WAIT_SECONDS, poll=ASYNC_POLL_SECONDS, environ=None):
    """
    The query as an asynchronous job on the ESA archive (OPS.183): it has no 2,000-row cap, and the file is the one `URL` gives. Starts the job, waits for it to complete, reads its result and
    deletes the job (the anonymous queue holds 50). With `PLANETGEN_ESA_USER` and `PLANETGEN_ESA_PASSWORD` set it logs in first. Raises `fetch.DownloadError` when the job fails or takes too long.
    """
    environ = os.environ if environ is None else environ
    user, password = environ.get(ESA_USER_ENV, ""), environ.get(ESA_PASSWORD_ENV, "")
    if user and password:
        secrets_mask.register(password)
        jar = http.cookiejar.CookieJar()
        opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar)).open
        _post(opener, "https://gea.esac.esa.int/tap-server/login", {"username": user, "password": password}, timeout)
    _body, job = _post(opener, ESA_TAP + "/async", {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "PHASE": "RUN", "QUERY": QUERY}, timeout)
    if "/async/" not in job:
        raise fetch.DownloadError(f"{ESA_TAP}/async: it gave no job address")
    deadline = time.monotonic() + wait
    try:
        while True:
            phase = _post(opener, job + "/phase", None, timeout)[0].decode("utf-8", "replace").strip().upper()
            if phase == "COMPLETED":
                return _post(opener, job + "/results/result", None, timeout)[0]
            if phase in ("ERROR", "ABORTED"):
                raise fetch.DownloadError(f"{job}: the job ended {phase.lower()}")
            if time.monotonic() > deadline:
                raise fetch.DownloadError(f"{job}: still {phase.lower() or 'waiting'} after {wait} seconds")
            time.sleep(poll)
    finally:
        try:
            _post(opener, job, {"ACTION": "DELETE"}, timeout)
        except fetch.DownloadError:
            pass


def _vizier_columns(table, names, extra=""):
    """The ADQL of the same sample asked of a VizieR table whose columns are `names` (id, ra, dec, parallax, its error, pmra, pmdec, G, BP-RP)."""
    ident, ra, dec, plx, e_plx, pmra, pmdec, g, bp_rp = names
    return (f'SELECT {ident}, {ra}, {dec}, {plx}, {e_plx}, {pmra}, {pmdec}, {g}, "{bp_rp}" FROM "{table}" WHERE {plx} > 40{extra} ORDER BY {ident}')


def _vizier(query):
    def get(opener, timeout):
        fields = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": "100000", "QUERY": query}
        return _post(opener, VIZIER_TAP, fields, timeout)[0]
    return get


_VIZIER_NAMES = {"source_id": ("source", "gaiaedr3", "gaiadr3"), "ra": ("ra_icrs", "raj2000", "ra"), "dec": ("de_icrs", "dej2000", "dec"),
                 "parallax": ("plx", "parallax"), "parallax_error": ("e_plx", "parallax_error"), "pmra": ("pmra",), "pmdec": ("pmde", "pmdec"),
                 "phot_g_mean_mag": ("gmag", "phot_g_mean_mag"), "bp_rp": ("bp-rp", "bp_rp")}


def vizier_to_file(data):
    """
    A VizieR csv answer as this sample's file: the columns under their Gaia archive names, in `COLUMNS` order, the rows sorted by `source_id` and kept only where the parallax is
    measured to better than 10 sigma (the cut a catalogue that only cuts on parallax leaves to be made here). Raises `ValueError` when a column cannot be found.
    """
    rows = list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
    if not rows:
        raise ValueError("no rows")
    names = {key.strip().lower(): key for key in rows[0]}
    found = {}
    for column, candidates in _VIZIER_NAMES.items():
        match = next((names[c] for c in candidates if c in names), None)
        if match is None:
            raise ValueError(f"no column for {column} among {sorted(names)}")
        found[column] = match
    kept = [row for row in rows if row[found["parallax"]].strip() and row[found["parallax_error"]].strip()
            and float(row[found["parallax"]]) > 40 and float(row[found["parallax"]]) / float(row[found["parallax_error"]]) > 10]
    kept.sort(key=lambda row: int(float(row[found["source_id"]])) if "e" in row[found["source_id"]].lower() else int(row[found["source_id"]]))
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(COLUMNS)
    for row in kept:
        writer.writerow([row[found[column]].strip() for column in COLUMNS])
    return out.getvalue().encode("utf-8")


def accept_rows(data):
    """`None` when `data` is a file of this sample (the header and between `MIN_ROWS` and `MAX_ROWS` rows), else the reason it is not (a capped answer is the usual one)."""
    lines = data.decode("utf-8", "replace").splitlines()
    if not lines or [c.strip() for c in lines[0].split(",")] != list(COLUMNS):
        return "its header is not the sample's columns"
    rows = len(lines) - 1
    if not MIN_ROWS <= rows <= MAX_ROWS:
        return f"it has {rows} rows where the sample has about 5,300 (a capped or an emptied answer)"
    return None


def describe_rows(path):
    """How many sources the file at `path` holds, for the line that says which source served it."""
    with open(path, "rb") as handle:
        return f"{sum(1 for _ in handle) - 1:,} sources"


MIRRORS = (
    fetch.Mirror("the ESA archive's asynchronous query", get=esa_async),
    fetch.Mirror("the Heidelberg ARI mirror", url=ARI_TAP + "?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": "100000", "QUERY": QUERY})),
    fetch.Mirror("VizieR (Gaia DR3, I/355)", get=_vizier(_vizier_columns(
        "I/355/gaiadr3", ("Source", "RA_ICRS", "DE_ICRS", "Plx", "e_Plx", "pmRA", "pmDE", "Gmag", "BP-RP"), " AND Plx/e_Plx > 10")), convert=vizier_to_file, accept=accept_rows),
    fetch.Mirror("VizieR (the Gaia nearby-star catalog GCNS, I/360)", get=_vizier(_vizier_columns(
        "I/360/gcnsmain", ("GaiaEDR3", "RA_ICRS", "DE_ICRS", "Plx", "e_Plx", "pmRA", "pmDE", "Gmag", "BP-RP"))), convert=vizier_to_file, accept=accept_rows),
)
"""tuple: What is tried, in order, when the archive's own synchronous query fails or is capped (OPS.183, Boss: have all four and if one fails go to the next). The first two give the file byte for byte
(held to the sha256); the VizieR ones are the same stars from another catalogue, taken when their file keeps to `accept_rows`."""

DATASET = fetch.Dataset(
    key="gaia_nearby",
    title="Gaia DR3 stars within 25 parsecs (the spacing of stars)",
    source="Gaia DR3, ESA/Gaia/DPAC, CC BY-SA 3.0 IGO, via the Gaia archive TAP service (or its Heidelberg mirror, or VizieR)",
    pieces=(fetch.Piece(FILE, URL, "c2ec659c7bd10ebf4e8b5a3a5d76271919626064f9b959176f61b04cf6c6cb7c", 766879, mirrors=MIRRORS, describe=describe_rows),),
    replaces="SEP1: the share of the galactic Hill spheres that random placement keeps between two systems",
    fill="SEP1",
    required=False,
)


@functools.lru_cache(maxsize=2)
def _positions(path):
    """The `(x, y, z)` of every source in parsecs, in the equatorial frame with the Sun at the origin, from the file at `path`."""
    out = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            distance = 1000.0 / float(row["parallax"])
            ra, dec = math.radians(float(row["ra"])), math.radians(float(row["dec"]))
            out.append((distance * math.cos(dec) * math.cos(ra), distance * math.cos(dec) * math.sin(ra), distance * math.sin(dec)))
    return tuple(out)


def positions_pc():
    """The `(x, y, z)` in parsecs of every Gaia source within 25 pc, the Sun at the origin."""
    return _positions(fetch.locate(DATASET).require().path(FILE))


def systems(positions=None, merge_pc=MERGE_PC):
    """
    The centers (numpy array, parsecs) of the systems the sources make: sources closer than `merge_pc` to one another, directly or through a chain, are one system at the mean of their places.
    `positions` is the list of `(x, y, z)` to use instead of the downloaded sample (a test's).
    """
    import numpy as np
    from scipy.spatial import cKDTree

    points = np.asarray(positions_pc() if positions is None else positions, dtype=float).reshape(-1, 3)
    parent = list(range(len(points)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in cKDTree(points).query_pairs(merge_pc):
        parent[find(i)] = find(j)
    groups = {}
    for i in range(len(points)):
        groups.setdefault(find(i), []).append(i)
    return np.array([points[members].mean(axis=0) for members in groups.values()])


def nearest_neighbour_distances_pc(centers=None, center_pc=CENTER_PC):
    """
    The distance to the nearest other system, parsecs, of every system within `center_pc` of the Sun (sorted), from the system centers `centers` (`systems()` when omitted). The neighbours are
    searched among all the systems of the sample, which is complete out to 25 pc, so a system within 15 pc never misses its nearest neighbour at the edge.
    """
    import numpy as np
    from scipy.spatial import cKDTree

    centers = systems() if centers is None else np.asarray(centers, dtype=float).reshape(-1, 3)
    inner = centers[np.linalg.norm(centers, axis=1) < center_pc]
    distances, _ = cKDTree(centers).query(inner, k=2)
    return np.sort(distances[:, 1])


def system_density_per_pc3(centers=None):
    """The number of systems per cubic parsec inside `LIMIT_PC`."""
    centers = systems() if centers is None else centers
    return len(centers) / (4.0 / 3.0 * math.pi * LIMIT_PC ** 3)
