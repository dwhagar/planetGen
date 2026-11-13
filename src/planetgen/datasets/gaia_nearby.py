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
import math
import urllib.parse

from planetgen.datasets import fetch

FILE = "gaia_dr3_25pc.csv"

QUERY = ("SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source "
         "WHERE parallax > 40 AND parallax_over_error > 10 ORDER BY source_id")
"""str: The ADQL query that makes the file: every DR3 source nearer than 25 pc (parallax over 40 mas) whose parallax is measured to better than 10 sigma."""

URL = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": QUERY})

LIMIT_PC = 25.0
"""float: The distance the sample reaches, parsecs (parallax 40 mas)."""

MERGE_PC = 0.1
"""float: Stars closer than this to a neighbour are one system (wide pairs and their companions, whose separations run to 20,000 AU), as the galaxy places a system once however many stars it holds."""

CENTER_PC = 15.0
"""float: Only systems nearer than this are asked for their nearest neighbour, so that the neighbour is always inside the complete sample (`LIMIT_PC`) and none is missed at its edge: the nearest neighbour of a system nearer than `CENTER_PC` is nearer than 10 pc."""

DATASET = fetch.Dataset(
    key="gaia_nearby",
    title="Gaia DR3 stars within 25 parsecs (the spacing of stars)",
    source="Gaia DR3, ESA/Gaia/DPAC, CC BY-SA 3.0 IGO, via the Gaia archive TAP service",
    pieces=(fetch.Piece(FILE, URL, "c2ec659c7bd10ebf4e8b5a3a5d76271919626064f9b959176f61b04cf6c6cb7c", 766879),),
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
