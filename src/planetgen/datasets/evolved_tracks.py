# planetgen/datasets/evolved_tracks.py

"""
The post-main-sequence tracks of the MIST grid (GEN.935): what a star does between the end of its main sequence and its
white dwarf or its collapse, read from the same isochrone archive as the main sequence (`isochrones.py`, GEN.548).

An isochrone file lists, for each of 107 ages, the stars of every initial mass in one list of equivalent evolutionary
points (EEPs, integers: 454 is the end of the main sequence, 605 the tip of the red giant branch, 707 the end of core
helium burning, 808 the start of the thermal pulses, 1409 the end of the asymptotic giant branch, 1710 the white dwarf
cooling sequence; [R0773, R1150]). At one EEP the initial mass falls as the age rises, so for a mass `m`
the age it reaches that EEP is found by interpolating across the isochrone ages in log mass. Doing that for every EEP
gives the star's own track: it resolves the Hertzsprung gap, the blue loops and the thermal-pulse stages, which the
isochrones (0.05 dex apart in age) cannot show as ages, because the EEPs are the isochrones' rows.

Setup keeps, for 130 log-spaced masses from 0.8 to 300 solar masses and every [Fe/H] of the grid, every sixth EEP from
the end of the main sequence on (and the landmark EEPs): the time since the main sequence ended as a share of the
main-sequence lifetime, the present mass, log L, log Teff, log R, the surface hydrogen and carbon fractions and the
phase (MIST's own: 2 red giant branch, 3 core helium burning, 4 early asymptotic giant branch, 5 thermal pulses, 6
post-asymptotic giant branch, 9 Wolf-Rayet). The readers blend the four nearest tracks (two masses, two metallicities)
at each EEP; a mass below the lightest track that reaches the thermal pulses takes that track, and a track that ends
earlier (a massive star ends at EEP 808, core helium exhaustion) is cut where it ends.
"""

import bisect
import math
from dataclasses import dataclass

import numpy as np

FIRST_EEP = 454
"""int: The end of the main sequence (terminal-age main sequence), the first EEP kept."""
RGB_TIP_EEP = 605
"""int: The tip of the red giant branch."""
CORE_HELIUM_END_EEP = 707
"""int: The end of core helium burning."""
THERMAL_PULSE_EEP = 808
"""int: The start of the thermal-pulsing asymptotic giant branch, and where a massive star's track ends."""
AGB_END_EEP = 1409
"""int: The end of the asymptotic giant branch (the post-asymptotic giant branch follows)."""
LAST_EEP = 1710
"""int: The last EEP of the grid (a white dwarf cooling sequence)."""
EEP_STRIDE = 6
"""int: Every sixth EEP is kept (the landmarks too)."""
SAMPLE_EEPS = tuple(sorted(set(range(FIRST_EEP, LAST_EEP + 1, EEP_STRIDE))
                           | {RGB_TIP_EEP, CORE_HELIUM_END_EEP, THERMAL_PULSE_EEP, AGB_END_EEP, LAST_EEP}))
"""tuple: The EEPs kept, ascending."""
MASS_NODES_SOL = tuple(float(m) for m in np.logspace(np.log10(0.8), np.log10(300.0), 130))
"""tuple: The initial masses (solar) the tracks are made for, log-spaced (0.02 dex)."""

LANDMARK_EEPS = (RGB_TIP_EEP, CORE_HELIUM_END_EEP, AGB_END_EEP, LAST_EEP)
"""tuple: The equivalent evolutionary points whose times `Grid.landmark_fracs` gives: the ends of the first ascent, of core helium
burning, of the asymptotic giant branch and of the track."""

COLUMNS = ("frac", "log_l", "star_mass", "log_teff", "log_r", "h1", "c12", "phase")
"""tuple: What each sample holds, in order."""
_FRAC, _LOG_L, _MASS, _LOG_TEFF, _LOG_R, _H1, _C12, _PHASE = range(8)

# The columns of a post-main-sequence row as `isochrones.parse_isochrone_file` returns it.
ROW_COLUMNS = ("eep", "log_age", "initial_mass", "star_mass", "log_l", "log_teff", "log_r", "h1", "c12", "phase")
_R_EEP, _R_LOG_AGE, _R_INITIAL, _R_MASS, _R_LOG_L, _R_LOG_TEFF, _R_LOG_R, _R_H1, _R_C12, _R_PHASE = range(10)


def sample_index(eep, count=None):
    """The index of the first sample at or after `eep`, held to the last (of the first `count` samples when given)."""
    last = len(SAMPLE_EEPS) if count is None else count
    return min(bisect.bisect_left(SAMPLE_EEPS, eep), last - 1)


def build(post_rows):
    """
    The tracks of one [Fe/H] from its post-main-sequence rows (an array with the columns of `ROW_COLUMNS`):
    `(tracks, ms_age_gy)` where `tracks` is a float32 array over the mass nodes, the sampled EEPs and `COLUMNS`, NaN
    where the star does not reach that EEP inside the grid's ages, and `ms_age_gy` the age (Gy, float64) each mass
    node leaves the main sequence at (NaN where the grid does not hold it).
    """
    nodes = np.log10(np.asarray(MASS_NODES_SOL))
    tracks = np.full((len(nodes), len(SAMPLE_EEPS), len(COLUMNS)), np.nan)
    ms_age = np.full(len(nodes), np.nan)
    if len(post_rows) == 0:
        return tracks.astype(np.float32), ms_age
    eeps = post_rows[:, _R_EEP].astype(int)
    order = np.argsort(eeps, kind="stable")
    rows, eeps = post_rows[order], eeps[order]
    unique, begin = np.unique(eeps, return_index=True)
    finish = np.append(begin[1:], len(eeps))
    spans = {int(e): (int(a), int(b)) for e, a, b in zip(unique, begin, finish)}
    at = {}

    def along(eep):
        """`(log10 initial mass, rows)` of one EEP in ascending initial mass, or None when there are fewer than two."""
        if eep not in at:
            a, b = spans.get(eep, (0, 0))
            if b - a < 2:
                at[eep] = None
            else:
                chosen = rows[a:b][np.argsort(rows[a:b, _R_INITIAL])]
                at[eep] = (np.log10(chosen[:, _R_INITIAL]), chosen)
        return at[eep]

    first = along(FIRST_EEP)
    if first is None:
        return tracks.astype(np.float32), ms_age
    ms_age = np.interp(nodes, first[0], first[1][:, _R_LOG_AGE], left=np.nan, right=np.nan)
    ms_age = 10.0 ** ms_age
    for j, eep in enumerate(SAMPLE_EEPS):
        found = along(eep)
        if found is None:
            continue
        log_m, chosen = found

        def read(column, log_m=log_m, chosen=chosen):
            return np.interp(nodes, log_m, chosen[:, column], left=np.nan, right=np.nan)

        age = 10.0 ** read(_R_LOG_AGE)
        tracks[:, j, _FRAC] = (age - ms_age) / ms_age
        tracks[:, j, _MASS] = read(_R_MASS)
        tracks[:, j, _LOG_L] = read(_R_LOG_L)
        tracks[:, j, _LOG_TEFF] = read(_R_LOG_TEFF)
        tracks[:, j, _LOG_R] = read(_R_LOG_R)
        tracks[:, j, _H1] = read(_R_H1)
        tracks[:, j, _C12] = read(_R_C12)
        position = np.interp(nodes, log_m, np.arange(len(log_m)), left=np.nan, right=np.nan)
        inside = ~np.isnan(position)
        phase = np.full(len(nodes), np.nan)
        phase[inside] = chosen[np.rint(position[inside]).astype(int), _R_PHASE]
        tracks[:, j, _PHASE] = phase
    _make_prefix_and_monotone(tracks)
    return tracks.astype(np.float32), ms_age


def _make_prefix_and_monotone(tracks):
    """Each node's track is cut where it first has a gap, and its time since the main sequence ended is made
    non-decreasing (the interpolation across isochrone ages can wobble by a rounding error)."""
    for node in range(tracks.shape[0]):
        valid = ~np.isnan(tracks[node, :, _LOG_L])
        if not valid.any():
            tracks[node] = np.nan
            continue
        stop = int(np.argmin(valid)) if not valid.all() else len(valid)
        tracks[node, stop:] = np.nan
        tracks[node, :stop, _FRAC] = np.maximum.accumulate(tracks[node, :stop, _FRAC])


@dataclass(frozen=True)
class Track:
    """
    The track of one initial mass and [Fe/H], blended from the four nearest of the grid: the arrays run over the
    sampled EEPs `eep` from the end of the main sequence to the end of the track. `frac` is the time since the main
    sequence ended as a share of its (MIST) lifetime `ms_age_gy`; `phase` is MIST's phase number.
    """

    mass_sol: float
    fe_h: float
    ms_age_gy: float
    eep: np.ndarray
    frac: np.ndarray
    star_mass: np.ndarray
    log_l: np.ndarray
    log_teff: np.ndarray
    log_r: np.ndarray
    h1: np.ndarray
    c12: np.ndarray
    phase: np.ndarray

    def __len__(self):
        return len(self.eep)

    def index_of(self, eep):
        """The index of the first sample at or after `eep`, held to the last."""
        return sample_index(eep, len(self.eep))

    def frac_at(self, eep):
        """The time since the main sequence ended, as a share of its lifetime, when the star reaches `eep` (held to the track)."""
        return float(self.frac[self.index_of(eep)])

    def reaches(self, eep):
        """Whether the track runs to `eep`."""
        return bool(self.eep[-1] >= eep)

    def at_frac(self, frac, first_eep, last_eep):
        """
        The state at `frac` of the main-sequence lifetime after it ended, between `first_eep` and `last_eep` (the
        window's EEPs, held to the track): `(star_mass, log_l, log_teff, log_r, h1, c12, phase)`. Values are linear
        in `frac` between samples; the phase is that of the sample at or before.
        """
        low, high = self.index_of(first_eep), self.index_of(last_eep)
        if high < low:
            high = low
        x = self.frac[low:high + 1]
        frac = min(max(frac, x[0]), x[-1])
        k = min(int(np.searchsorted(x, frac, side="right")) - 1, len(x) - 2) if len(x) > 1 else 0
        k = max(k, 0)
        if len(x) == 1:
            weight = 0.0
        else:
            span = x[k + 1] - x[k]
            weight = 0.0 if span <= 0.0 else (frac - x[k]) / span
        i = low + k

        def mix(values):
            if len(x) == 1:
                return float(values[i])
            return float(values[i] * (1.0 - weight) + values[i + 1] * weight)

        return (mix(self.star_mass), mix(self.log_l), mix(self.log_teff), mix(self.log_r), mix(self.h1), mix(self.c12),
                int(self.phase[i]))


class Grid:
    """The tracks of every [Fe/H] of the grid, with the mass nodes and the main-sequence age of each, as the derived file holds them."""

    def __init__(self, fe_h, tracks, ms_age_gy):
        self.fe_h = tuple(float(v) for v in fe_h)
        self.tracks = tracks                                  # [fe_h][node, sample, column] float32
        self.ms_age_gy = ms_age_gy                            # [fe_h][node] float64
        self.nodes = np.log10(np.asarray(MASS_NODES_SOL))
        self._nodes = [float(v) for v in self.nodes]
        eep = np.asarray(SAMPLE_EEPS)
        stop = int(np.searchsorted(eep, THERMAL_PULSE_EEP))
        self._counts = [[int(c) for c in (~np.isnan(found[:, :, _LOG_L])).sum(axis=1)] for found in tracks]
        self._lightest = []
        for found in tracks:
            reaching = np.nonzero(~np.isnan(found[:, stop, _LOG_L]))[0]
            self._lightest.append(int(reaching[0]) if len(reaching) else None)
        # The share of the main-sequence lifetime at the landmarks of every node (held to the node's last sample), which
        # is all a window needs.
        self._helium = False
        self._landmarks = []
        for found, counts in zip(tracks, self._counts):
            table = np.full((found.shape[0], len(LANDMARK_EEPS)), np.nan)
            for node, count in enumerate(counts):
                if count:
                    for k, landmark in enumerate(LANDMARK_EEPS):
                        table[node, k] = found[node, min(SAMPLE_EEPS.index(landmark), count - 1), _FRAC]
            self._landmarks.append(table)

    def has_tracks(self):
        return any(index is not None for index in self._lightest)

    def helium_star_relation(self):
        """
        `(intercept, slope, sigma)` of `log10 L/Lsun = intercept + slope log10 (M/Msun)` for the stars MIST calls Wolf-Rayet
        (phase 9: hydrogen all but gone from the surface), fitted over every sample of every track of the grid, with the
        standard deviation of the residual (dex); None when no track has one. The luminosity of a stripped helium star follows
        its mass, whatever the star it came from (GEN.939).
        """
        if self._helium is False:
            masses, log_ls = [], []
            for tracks in self.tracks:
                stripped = (tracks[:, :, _PHASE] == 9.0) & np.isfinite(tracks[:, :, _MASS]) & (tracks[:, :, _MASS] > 0.0)
                masses.append(tracks[:, :, _MASS][stripped].astype(float))
                log_ls.append(tracks[:, :, _LOG_L][stripped].astype(float))
            mass, log_l = np.concatenate(masses), np.concatenate(log_ls)
            if len(mass) < 10:
                self._helium = None
            else:
                slope, intercept = np.polyfit(np.log10(mass), log_l, 1)
                sigma = float(np.std(log_l - (intercept + slope * np.log10(mass))))
                self._helium = (float(intercept), float(slope), sigma)
        return self._helium

    def _bracket(self, fe_h):
        values = self.fe_h
        fe_h = min(max(fe_h, values[0]), values[-1])
        upper = min(bisect.bisect_left(values, fe_h), len(values) - 1)
        lower = upper if values[upper] == fe_h else upper - 1
        if lower == upper:
            return lower, upper, 0.0
        return lower, upper, (fe_h - values[lower]) / (values[upper] - values[lower])

    def _corners(self, mass_sol, fe_h):
        """The tracks that make this one, `(index, low node, high node, share of the high one, weight)`, and how many
        samples all of them have; None when no track of the grid reaches the thermal pulses."""
        lower, upper, weight_fe = self._bracket(fe_h)
        log_m = math.log10(mass_sol)
        nodes = self._nodes
        corners = []
        for index, weight in ((lower, 1.0 - weight_fe), (upper, weight_fe)):
            if weight <= 0.0 or self._lightest[index] is None:
                continue
            high = min(max(bisect.bisect_left(nodes, log_m), 1), len(nodes) - 1)
            low = high - 1
            along = min(max((log_m - nodes[low]) / (nodes[high] - nodes[low]), 0.0), 1.0)
            if low < self._lightest[index]:      # lighter than the lightest track that runs on: that track
                low = high = self._lightest[index]
                along = 0.0
            corners.append((index, low, high, along, weight))
        if not corners:
            return None, 0
        length = min(self._counts[index][node] for index, low, high, _, _ in corners for node in (low, high))
        return corners, length

    def track(self, mass_sol, fe_h):
        """The `Track` of this mass and [Fe/H]; None when no track of the grid reaches the thermal pulses."""
        corners, length = self._corners(mass_sol, fe_h)
        if corners is None:
            return None
        total = sum(weight for *_, weight in corners)
        blended = np.zeros((length, len(COLUMNS) - 1))
        age = 0.0
        best = max(corners, key=lambda corner: corner[4])
        for index, low, high, along, weight in corners:
            data = self.tracks[index]
            mix = data[low, :length, :_PHASE] * (1.0 - along) + data[high, :length, :_PHASE] * along
            if low == high and mass_sol < 10.0 ** self._nodes[low]:
                mix[:, _MASS] *= mass_sol / 10.0 ** self._nodes[low]      # the lightest track, for a star lighter than it: its present mass in proportion
            blended += mix * (weight / total)
            age += (self.ms_age_gy[index][low] * (1.0 - along) + self.ms_age_gy[index][high] * along) * (weight / total)
        index, low, high, along, _ = best
        pick = high if along >= 0.5 else low
        phase = self.tracks[index][pick, :length, _PHASE].astype(int)
        return Track(mass_sol, fe_h, float(age), np.asarray(SAMPLE_EEPS[:length]), blended[:, _FRAC], blended[:, _MASS],
                     blended[:, _LOG_L], blended[:, _LOG_TEFF], blended[:, _LOG_R], blended[:, _H1], blended[:, _C12], phase)

    def light(self, mass_sol, fe_h):
        """`(frac, log_l)` of the track of this mass and [Fe/H]: the two columns the bright windows read, blended as `track` does
        and cheaper; None when no track of the grid reaches the thermal pulses."""
        corners, length = self._corners(mass_sol, fe_h)
        if corners is None:
            return None
        total = sum(weight for *_, weight in corners)
        blended = np.zeros((length, 2))
        for index, low, high, along, weight in corners:
            data = self.tracks[index]
            blended += (data[low, :length, :2] * (1.0 - along) + data[high, :length, :2] * along) * (weight / total)
        return blended[:, 0], blended[:, 1]

    def landmark_fracs(self, mass_sol, fe_h):
        """The time since the main sequence ended, as a share of its lifetime, when a star of this mass and [Fe/H] reaches each
        of `LANDMARK_EEPS` (held to the track's last sample), blended as `track` does; None when no track of the grid
        reaches the thermal pulses."""
        corners, _ = self._corners(mass_sol, fe_h)
        if corners is None:
            return None
        total = sum(weight for *_, weight in corners)
        found = np.zeros(len(LANDMARK_EEPS))
        for index, low, high, along, weight in corners:
            table = self._landmarks[index]
            found += (table[low] * (1.0 - along) + table[high] * along) * (weight / total)
        return tuple(float(v) for v in found)
