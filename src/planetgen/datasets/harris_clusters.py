# planetgen/datasets/harris_clusters.py

"""
The Harris catalog of parameters for Milky Way globular clusters (GEN.297): 157 clusters with positions, distances,
metallicities, integrated magnitudes and colours, radial velocities, velocity dispersions and structural parameters.

The file is `mwgc.dat`, the December 2010 revision of the catalog, from Harris's page at McMaster
(https://physwww.mcmaster.ca/~harris/mwgc.dat; [R1827]; cite it as "[R1827] (2010 edition)"). It is three
fixed-width tables, one per part (positions; metallicity and photometry; velocities and structure), each row a cluster. The reader
joins them by cluster name into one `Cluster` per row; a blank cell is None. Until now the cluster laws of `tuning.py` were
"recalled from [R1827]"; this is the table they are checked against and fitted to.

The reader raises `fetch.DataMissing` when the file has not been downloaded.
"""

import functools
import re
from dataclasses import dataclass
from typing import Optional

from planetgen.datasets import fetch

FILE = "mwgc.dat"

DATASET = fetch.Dataset(
    key="harris_clusters",
    title="Harris catalog of parameters for Milky Way globular clusters (2010 edition)",
    source="[R1827] (2010 edition), McMaster University",
    pieces=(fetch.Piece(FILE, "https://physwww.mcmaster.ca/~harris/mwgc.dat",
                        "b457a690adc422fb2e8e05a23e1de50b62d1d5ac10dce99963c956c1b17ca658", 49364),),
    replaces="the globular cluster laws of tuning.py, recalled from [R1827]",
    fill="SC4",
    fallback="the recalled medians of the cluster laws in tuning.py",
)

PARTS = ("Part I:", "Part II:", "Part III:")
_NAME_WIDTH = 12
_PART_I = (("alias", 12, 25), ("ra", 25, 37), ("dec", 37, 50), ("l_deg", 50, 58), ("b_deg", 58, 66), ("r_sun_kpc", 66, 73),
           ("r_gc_kpc", 73, 79), ("x_kpc", 79, 85), ("y_kpc", 85, 91), ("z_kpc", 91, 97))
_PART_II = (("fe_h", 12, 18), ("fe_h_weight", 18, 21), ("e_bv", 21, 28), ("v_hb", 28, 34), ("distance_modulus_v", 34, 40),
            ("v_total", 40, 46), ("m_v_total", 46, 53), ("u_b", 53, 60), ("b_v", 60, 66), ("v_r_colour", 66, 72),
            ("v_i", 72, 78), ("spectral_type", 78, 85), ("ellipticity", 85, 91))
_PART_III = (("v_r_km_s", 12, 19), ("v_r_error_km_s", 19, 25), ("v_lsr_km_s", 25, 33), ("sigma_v_km_s", 33, 41),
             ("sigma_v_error_km_s", 41, 47), ("concentration", 47, 55), ("collapse_flag", 55, 59), ("core_radius_arcmin", 59, 64),
             ("half_light_radius_arcmin", 64, 70), ("mu_v", 70, 78), ("log_rho_0", 78, 85), ("log_core_relaxation_yr", 85, 92),
             ("log_median_relaxation_yr", 92, 98))
_TEXT = {"alias", "ra", "dec", "spectral_type", "collapse_flag"}


@dataclass(frozen=True)
class Cluster:
    """One cluster. Distances are in kiloparsecs (the Galactic centre at 8.0 kpc), radii in arcminutes, magnitudes as printed;
    `ra_deg` and `dec_deg` are J2000. A field the catalog leaves blank is None."""

    name: str
    alias: str
    ra_deg: Optional[float]
    dec_deg: Optional[float]
    l_deg: Optional[float]
    b_deg: Optional[float]
    r_sun_kpc: Optional[float]
    r_gc_kpc: Optional[float]
    x_kpc: Optional[float]
    y_kpc: Optional[float]
    z_kpc: Optional[float]
    fe_h: Optional[float]
    fe_h_weight: Optional[float]
    e_bv: Optional[float]
    v_hb: Optional[float]
    distance_modulus_v: Optional[float]
    v_total: Optional[float]
    m_v_total: Optional[float]
    u_b: Optional[float]
    b_v: Optional[float]
    v_r_colour: Optional[float]
    v_i: Optional[float]
    spectral_type: str
    ellipticity: Optional[float]
    v_r_km_s: Optional[float]
    v_r_error_km_s: Optional[float]
    v_lsr_km_s: Optional[float]
    sigma_v_km_s: Optional[float]
    sigma_v_error_km_s: Optional[float]
    concentration: Optional[float]
    core_collapsed: bool
    core_radius_arcmin: Optional[float]
    half_light_radius_arcmin: Optional[float]
    mu_v: Optional[float]
    log_rho_0: Optional[float]
    log_core_relaxation_yr: Optional[float]
    log_median_relaxation_yr: Optional[float]


def _number(text):
    text = text.strip()
    return float(text) if text else None


def _sexagesimal(text, hours):
    """Degrees from `"00 24 05.67"` (hours when `hours`) or `"-72 04 52.6"`; None when blank."""
    parts = text.split()
    if len(parts) != 3:
        return None
    sign = -1.0 if parts[0].startswith("-") else 1.0
    first, minutes, seconds = abs(float(parts[0])), float(parts[1]), float(parts[2])
    value = sign * (first + minutes / 60.0 + seconds / 3600.0)
    return value * 15.0 if hours else value


def _rows(lines, part):
    """The data rows of one part: every line after the `ID` header line up to the next rule, as `(name, line)`."""
    start = next(i for i, line in enumerate(lines) if line.strip().startswith(part))
    header = next(i for i in range(start, len(lines)) if re.match(r"\s+ID\s", lines[i]))
    rows = []
    for line in lines[header + 1:]:
        if line.startswith("_"):
            break
        if line.strip():
            rows.append((line[:_NAME_WIDTH].strip(), line))
    return rows


def _fields(columns, line):
    return {name: line[start:end].strip() if name in _TEXT else _number(line[start:end]) for name, start, end in columns}


def parse(text):
    """The clusters of a `mwgc.dat`, in the file's order, each joined from its three parts by name."""
    lines = text.splitlines()
    one, two, three = ({name: _fields(columns, line) for name, line in _rows(lines, part)}
                       for part, columns in zip(PARTS, (_PART_I, _PART_II, _PART_III)))
    if not one or set(one) != set(two) or set(one) != set(three):
        raise ValueError("the three parts of the Harris catalog do not list the same clusters")
    clusters = []
    for name, first in one.items():
        flag = three[name]["collapse_flag"]
        clusters.append(Cluster(
            name=name, alias=first["alias"], ra_deg=_sexagesimal(first["ra"], True), dec_deg=_sexagesimal(first["dec"], False),
            **{key: first[key] for key in ("l_deg", "b_deg", "r_sun_kpc", "r_gc_kpc", "x_kpc", "y_kpc", "z_kpc")},
            **{key: value for key, value in two[name].items()},
            **{key: value for key, value in three[name].items() if key != "collapse_flag"},
            core_collapsed=flag.startswith("c")))
    return tuple(clusters)


@functools.lru_cache(maxsize=2)
def _load(path):
    with open(path, encoding="utf-8", errors="replace") as handle:
        return parse(handle.read())


def clusters():
    """The 157 clusters of the catalog, in the file's order."""
    return _load(fetch.locate(DATASET).require().path(FILE))
