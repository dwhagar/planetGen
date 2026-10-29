# planetgen/datasets/hhe_eos.py

"""
The hydrogen-helium equation of state of giant planets and brown dwarfs (GEN.733, registry fill GS1).

Source: Chabrier and Debras 2021 (ApJ, "A new equation of state for dense hydrogen-helium mixtures II: taking into account
hydrogen-helium interactions") and Chabrier, Mazevet and Soubiran 2019 (ApJ 872, 51), the (T, P) tables of `DirEOS2021.tar.gz`
(16 MB, 12 text files) on Gilles Chabrier's page at the ENS Lyon (perso.ens-lyon.fr/gilles.chabrier/DirEOS). Three tables are
read, each a rectangular grid of 121 isotherms (log T = 2.0 to 8.0, step 0.05) of 441 pressures (log P = -9.0 to 13.0, P in GPa,
step 0.05), with columns log T, log P, log rho [g/cm3], log U [MJ/kg], log S [MJ/kg/K], the response functions and the adiabatic
gradient:

- `TABLEEOS_2021_TP_Y0275_v1`, the mixture of helium mass fraction Y = 0.275 with the hydrogen-helium interactions;
- `TABLE_H_TP_v1` and `TABLE_HE_TP_v1`, pure hydrogen and pure helium.

Setup downloads the archive once, checks its sha256, keeps the density and entropy of the mixture and the share its volume moves by
for a change of helium (`ivl_slope`, below) in `hhe_eos.npz` (about 0.5 MB) and removes the archive. The tables' own warnings hold: values
in the solid phase, beyond the first-order quantum correction of the ions and at very low density are unphysical, and the
interactions are uncertain outside the temperature-density range of the simulations they come from (the README of the archive).

The reader gives the density at any helium mass fraction by the ideal volume law, the way Chabrier et al. 2019 mix hydrogen
and helium: the mixture table at 0.275 with its specific volume moved by `(Y - 0.275)` times `ivl_slope`, the volume the two pure
fluids differ by over the mixture's, which keeps the interaction term of the 0.275 table (a stand-in of the structure model, GS3).
The slope is held to `IVL_SLOPE_RANGE`, since the pure tables hold unphysical values where a fluid is solid or the density is
tiny. An adiabat is a line of the mixture table's entropy, not of its adiabatic gradient, which the tables' own README says carries
numerical noise (`adiabat_log_t`). A point off the grid takes the value on its nearest edge (no extrapolation); the grid starts at
100 K, and the first scale height of a colder atmosphere is read at 100 K.

The readers raise `fetch.DataMissing` when the file has not been downloaded; the structure model (`physics/giant_structure.py`)
asks `Located.use_stand_in` and a run told to disregard missing data takes the polytrope profile instead (GEN.599).
"""

import functools
import io
import math
import tarfile
from dataclasses import dataclass, field

import numpy as np

from planetgen.datasets import fetch

ARCHIVE = "DirEOS2021.tar.gz"
DERIVED = "hhe_eos.npz"

MIXTURE_MEMBER = "TABLEEOS_2021_TP_Y0275_v1"
HYDROGEN_MEMBER = "TABLE_H_TP_v1"
HELIUM_MEMBER = "TABLE_HE_TP_v1"
MEMBERS = (MIXTURE_MEMBER, HYDROGEN_MEMBER, HELIUM_MEMBER)
"""tuple: The archive members read."""

IVL_SLOPE_RANGE = (-1.2, 3.0)
"""tuple: The bounds of `ivl_slope`, the pure fluids' volume difference over the mixture's volume (a bound of the ideal volume law is
-1/(1 - Y) and 1/Y; the real helium volume is about half of hydrogen's, a slope of about -0.5)."""

ONE_BAR_LOG_P_GPA = -4.0
"""float: log10 of the one-bar level in the tables' units (GPa)."""

DIATOMIC_ADIABATIC_GRADIENT = 2.0 / 7.0
"""float: d ln T / d ln P of an ideal diatomic gas (molecular hydrogen above its rotational freeze-out), for the part of a cold adiabat
under the table's 100 K edge."""

TABLE_HELIUM_SHARE = 0.275
"""float: The helium mass fraction Y of the mixture table (Chabrier and Debras 2021)."""

DATASET = fetch.Dataset(
    key="hhe_eos",
    title="Chabrier and Debras 2021 hydrogen-helium equation of state, (T, P) tables",
    source="Chabrier and Debras 2021 (ApJ) and Chabrier, Mazevet and Soubiran 2019 (ApJ 872, 51), DirEOS2021, ENS Lyon",
    pieces=(fetch.Piece(ARCHIVE, "https://perso.ens-lyon.fr/gilles.chabrier/DirEOS/" + ARCHIVE,
                        "45f316790ce20d5d1ce0abee4db308521b5bfdc5526d0997141a2784834feeff", 16367726, transient=True),),
    replaces="the n = 1 polytrope of a gas giant, the n = 3/2 polytrope of a brown dwarf and the constant-density layers of an ice giant",
    fill="GS1",
    fallback="the polytrope profile of the layer pressures (n = 1 gas giant, n = 3/2 brown dwarf, constant-density ice giant)",
    derived=(DERIVED,),
    derive=lambda folder: derive(folder),
    derive_version=2,
)

_COLUMNS = {"log_t": 0, "log_p": 1, "log_rho": 2, "log_s": 4}
_WIDTH = 10


def _read_table(data, name):
    """`(log_t axis, log_p axis, columns)` of one table's text (`data` bytes): isotherms of ten-column rows, `columns` a dict of
    `(n_t, n_p)` float64 arrays by `_COLUMNS` name."""
    rows = np.loadtxt(io.BytesIO(data), comments="#", dtype=np.float64, ndmin=2)
    if rows.shape[1] != _WIDTH:
        raise ValueError(f"{name}: {rows.shape[1]} columns, expected {_WIDTH}")
    log_t = np.unique(rows[:, _COLUMNS["log_t"]])
    log_p = np.unique(rows[:, _COLUMNS["log_p"]])
    if len(log_t) < 2 or len(log_p) < 2 or len(rows) != len(log_t) * len(log_p):
        raise ValueError(f"{name}: {len(rows)} rows do not fill a {len(log_t)} by {len(log_p)} grid")
    grid = rows.reshape(len(log_t), len(log_p), _WIDTH)
    if not (np.allclose(grid[:, 0, 0], log_t) and np.allclose(grid[0, :, 1], log_p)):
        raise ValueError(f"{name}: the rows are not isotherm by isotherm with rising pressure")
    return log_t, log_p, {key: grid[:, :, column] for key, column in _COLUMNS.items() if key not in ("log_t", "log_p")}


def derive(folder):
    """The compact file from the downloaded archive: `{DERIVED: npz bytes}` with the two axes, the mixture's `log_rho` and `log_s`, and
    `ivl_slope`, the specific volume of pure helium less that of pure hydrogen over the ideal-volume-law mixture's at Y = 0.275,
    held to `IVL_SLOPE_RANGE`."""
    tables = {}
    with tarfile.open(f"{folder}/{ARCHIVE}") as archive:
        for member in archive:
            name = member.name.rsplit("/", 1)[-1]
            if member.isfile() and name in MEMBERS:
                tables[name] = _read_table(archive.extractfile(member).read(), name)
    missing = [name for name in MEMBERS if name not in tables]
    if missing:
        raise ValueError(f"the archive lacks {', '.join(missing)}")
    log_t, log_p, mixture = tables[MIXTURE_MEMBER]
    for name in (HYDROGEN_MEMBER, HELIUM_MEMBER):
        if not (np.allclose(tables[name][0], log_t) and np.allclose(tables[name][1], log_p)):
            raise ValueError(f"{name}: not on the mixture table's grid")
    v_h, v_he = 10.0 ** -tables[HYDROGEN_MEMBER][2]["log_rho"], 10.0 ** -tables[HELIUM_MEMBER][2]["log_rho"]
    ivl_slope = np.clip((v_he - v_h) / ((1.0 - TABLE_HELIUM_SHARE) * v_h + TABLE_HELIUM_SHARE * v_he), *IVL_SLOPE_RANGE)
    if not np.all(np.isfinite(ivl_slope)) or not all(np.all(np.isfinite(mixture[key])) for key in mixture):
        raise ValueError("the tables hold a value that is not a number")
    buffer = io.BytesIO()
    np.savez_compressed(buffer, log_t=log_t, log_p=log_p, log_rho=mixture["log_rho"].astype(np.float32),
                        log_s=mixture["log_s"].astype(np.float32), ivl_slope=ivl_slope.astype(np.float32))
    return {DERIVED: buffer.getvalue()}


@dataclass(frozen=True)
class Grid:
    """
    The tables on their grid. The integrator of `physics/giant_structure.py` reads them thousands of times a body, and plain floats in
    flat tuples are faster than numpy scalars, so the density and the slope are kept as flat row-major tuples, and the entropy
    also as an array for `adiabat_log_t`.

    Attributes:
        log_t0, d_log_t, n_t (float, float, int): The first log10 temperature (K), the step and the number of isotherms.
        log_p0, d_log_p, n_p (float, float, int): The first log10 pressure (GPa), the step and the number of pressures.
        log_rho, ivl_slope (tuple): Each `n_t * n_p` values, isotherm by isotherm: log10 density of the Y = 0.275 mixture (g/cm3)
            and the ideal-volume-law slope of its specific volume with helium share.
        log_s (ndarray): log10 specific entropy of the mixture (MJ/kg/K), `(n_t, n_p)`.
    """

    log_t0: float
    d_log_t: float
    n_t: int
    log_p0: float
    d_log_p: float
    n_p: int
    log_rho: tuple
    ivl_slope: tuple
    log_s: np.ndarray = field(compare=False, repr=False)

    @property
    def log_t_range(self):
        return self.log_t0, self.log_t0 + self.d_log_t * (self.n_t - 1)

    @property
    def log_p_range_gpa(self):
        return self.log_p0, self.log_p0 + self.d_log_p * (self.n_p - 1)

    def sample(self, log_t, log_p_gpa):
        """`(log10 density of the mixture in g/cm3, ivl_slope)` at the point, bilinear in (log T, log P) and held to the grid's edge
        outside it."""
        x = (log_t - self.log_t0) / self.d_log_t
        y = (log_p_gpa - self.log_p0) / self.d_log_p
        x = 0.0 if x < 0.0 else (self.n_t - 1.0 if x > self.n_t - 1.0 else x)
        y = 0.0 if y < 0.0 else (self.n_p - 1.0 if y > self.n_p - 1.0 else y)
        i = int(x)
        j = int(y)
        if i > self.n_t - 2:
            i = self.n_t - 2
        if j > self.n_p - 2:
            j = self.n_p - 2
        fx = x - i
        fy = y - j
        n_p = self.n_p
        a = i * n_p + j
        b = a + n_p
        w00 = (1.0 - fx) * (1.0 - fy)
        w10 = fx * (1.0 - fy)
        w01 = (1.0 - fx) * fy
        w11 = fx * fy
        rho, slope = self.log_rho, self.ivl_slope
        return (rho[a] * w00 + rho[b] * w10 + rho[a + 1] * w01 + rho[b + 1] * w11,
                slope[a] * w00 + slope[b] * w10 + slope[a + 1] * w01 + slope[b + 1] * w11)

    def log_entropy(self, log_t, log_p_gpa):
        """log10 of the mixture's specific entropy (MJ/kg/K) at the point, bilinear and held to the edge outside the grid."""
        x = min(max((log_t - self.log_t0) / self.d_log_t, 0.0), self.n_t - 1.0)
        y = min(max((log_p_gpa - self.log_p0) / self.d_log_p, 0.0), self.n_p - 1.0)
        i = min(int(x), self.n_t - 2)
        j = min(int(y), self.n_p - 2)
        fx, fy = x - i, y - j
        s = self.log_s
        return float(s[i, j] * (1.0 - fx) * (1.0 - fy) + s[i + 1, j] * fx * (1.0 - fy) + s[i, j + 1] * (1.0 - fx) * fy
                     + s[i + 1, j + 1] * fx * fy)

    def log_entropy_of_surface_temperature(self, temperature_k):
        """
        log10 entropy (MJ/kg/K) of the adiabat that is at `temperature_k` at the one-bar level. Under the table's coldest isotherm
        (100 K) the line is carried up the molecular gas's own adiabat, T proportional to P^(2/7), to where it reaches 100 K, and the
        entropy there is read (the table holds nothing colder, GS5).
        """
        cold = 10.0 ** self.log_t0
        if temperature_k >= cold:
            return self.log_entropy(math.log10(temperature_k), ONE_BAR_LOG_P_GPA)
        return self.log_entropy(self.log_t0, ONE_BAR_LOG_P_GPA + math.log10(cold / temperature_k) / DIATOMIC_ADIABATIC_GRADIENT)

    def surface_temperature_k(self, log_s0):
        """The one-bar temperature (K) of the adiabat of entropy `log_s0`, the inverse of `log_entropy_of_surface_temperature`: the
        temperature of the one-bar column where it is hotter than the table's coldest isotherm, else the pressure along that isotherm
        where the entropy is `log_s0`, taken down to one bar along the molecular gas's adiabat."""
        t = self.adiabat_log_t(log_s0)
        j = int(round((ONE_BAR_LOG_P_GPA - self.log_p0) / self.d_log_p))
        if t[j] > self.log_t0 + 1e-9:
            return 10.0 ** t[j]
        for k in range(j, self.n_p):                      # the first pressure where the adiabat is on the table
            if t[k] > self.log_t0 + 1e-9:
                log_p = self.log_p0 + k * self.d_log_p
                return 10.0 ** self.log_t0 * 10.0 ** (DIATOMIC_ADIABATIC_GRADIENT * (ONE_BAR_LOG_P_GPA - log_p))
        return 10.0 ** self.log_t0

    def adiabat_log_t(self, log_s0):
        """
        log10 of the temperature (K) at each pressure node of the table (`log_p0 + k * d_log_p`, GPa) along the line of constant
        entropy `log_s0` (log10 MJ/kg/K): in each column the first isotherm whose entropy passes it, interpolated; held to the grid's
        coldest and hottest isotherm where the entropy is outside the column's range.
        """
        s = self.log_s
        above = s > log_s0
        none = ~above.any(axis=0)
        upper = np.where(none, self.n_t - 1, above.argmax(axis=0))        # the first isotherm hotter than the adiabat
        lower = np.maximum(upper - 1, 0)
        columns = np.arange(self.n_p)
        s_lo, s_hi = s[lower, columns], s[upper, columns]
        span = s_hi - s_lo
        fraction = np.where((span > 0.0) & (upper > lower), np.clip((log_s0 - s_lo) / np.where(span > 0.0, span, 1.0), 0.0, 1.0), 0.0)
        return tuple((self.log_t0 + (lower + np.where(none, 1.0, fraction)) * self.d_log_t).tolist())


@functools.lru_cache(maxsize=2)
def _load(path):
    with np.load(path) as data:
        log_t, log_p = data["log_t"], data["log_p"]
        arrays = {name: tuple(float(v) for v in data[name].reshape(-1)) for name in ("log_rho", "ivl_slope")}
        log_s = data["log_s"].astype(np.float64)
    if len(log_t) < 2 or len(log_p) < 2:
        raise ValueError("the hydrogen-helium table has no grid")
    return Grid(float(log_t[0]), float(log_t[1] - log_t[0]), len(log_t), float(log_p[0]), float(log_p[1] - log_p[0]), len(log_p),
                log_s=log_s, **arrays)


def grid():
    """The tables, loaded once. Raises `fetch.DataMissing` when they were not downloaded."""
    return _load(fetch.locate(DATASET).require().path(DERIVED))


def density_kg_m3(temperature_k, pressure_pa, helium_share=TABLE_HELIUM_SHARE):
    """The density of a hydrogen-helium fluid of helium mass fraction `helium_share` at the temperature and pressure: the
    Y = 0.275 interacting mixture with its volume moved by the ideal volume law for the difference in Y."""
    log_rho, slope = grid().sample(math.log10(temperature_k), math.log10(pressure_pa) - 9.0)
    return 1000.0 / (10.0 ** -log_rho * (1.0 + (helium_share - TABLE_HELIUM_SHARE) * slope))


def entropy_mj_kg_k(temperature_k, pressure_pa):
    """The specific entropy of the Y = 0.275 mixture at the temperature and pressure (MJ/kg/K)."""
    return 10.0 ** grid().log_entropy(math.log10(temperature_k), math.log10(pressure_pa) - 9.0)
