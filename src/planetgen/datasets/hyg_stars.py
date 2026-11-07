# planetgen/datasets/hyg_stars.py

"""
The HYG star database v4.1 (VIEW.30, S12 of the system-sky plan): 119,626 rows of the Hipparcos, Yale Bright Star and Gliese catalogs, each with its apparent V magnitude, distance and
position. It is a TEST ORACLE only: the real sky the generated sky is held against (`tests/test_sky_hyg_calibration.py`: about 8,900 stars at V <= 6.5 and 25,800 at V <= 7.5). Nothing in
the generator reads it, so its dataset is `required` False: setup downloads it with the rest but bulk generation does not require it.

The file is `hygdata_v41.csv`, from the archive of the HYG-Database repository (the project moved to Codeberg and the archive, which holds this file, takes no more updates, so the hash
does not change). The licence is CC BY-SA 4.0, share-alike, so the file is never put in the repository: it goes to the cache folder, like the other datasets.

The reader raises `fetch.DataMissing` when the file has not been downloaded.
"""

import csv
import functools

from planetgen.datasets import fetch

FILE = "hygdata_v41.csv"

DATASET = fetch.Dataset(
    key="hyg_stars",
    title="HYG star database v4.1 (test oracle for the generated sky)",
    source="HYG database v4.1 (Hipparcos, Yale Bright Star and Gliese catalogs), CC BY-SA 4.0, astronexus",
    pieces=(fetch.Piece(FILE, "https://raw.githubusercontent.com/astronexus/HYG-Database/master/hyg/CURRENT/hygdata_v41.csv",
                        "d9f69fd86bbf90a4e4d52b4c5c53eacfa6dfc0bfdef85bfd94f095e0bebe4ebd", 33932548),),
    replaces="nothing: the real sky the generated sky is checked against in a test",
    fill="the generated-sky check",
    required=False,
)


@functools.lru_cache(maxsize=2)
def _magnitudes(path):
    """`(apparent V magnitude, distance in parsecs or None)` of every star but the Sun, from the file at `path`."""
    stars = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["id"] == "0" or not row["mag"]:
                continue
            distance = float(row["dist"]) if row["dist"] else None
            stars.append((float(row["mag"]), distance if distance is not None and distance < 100000.0 else None))
    return tuple(stars)


def stars():
    """The `(apparent V magnitude, distance in parsecs or None)` of every star of the catalogue except the Sun (a distance of 100,000 pc is the catalogue's "unknown")."""
    return _magnitudes(fetch.locate(DATASET).require().path(FILE))


def count_at_or_brighter(limit_mag, known_distance=False):
    """How many stars the catalogue lists at an apparent V magnitude of `limit_mag` or brighter (only those with a known distance when `known_distance`)."""
    return sum(1 for magnitude, distance in stars() if magnitude <= limit_mag and (distance is not None or not known_distance))
