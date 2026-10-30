# planetgen/physics/post_main_sequence.py

"""
A star after its main sequence (GEN.935): the age windows of its subgiant, giant, helium-burning and asymptotic giant
or supergiant phases, and its luminosity, temperature and radius in each, read from the MIST tracks
(`datasets/evolved_tracks.py`) instead of drawn.

The windows and the tracks come from two sources, as the main sequence's do (`main_sequence`): the generator's own phase
times say when a star is in a phase, and MIST says what the star looks like there. The main-sequence lifetime is
Hurley, Pols and Tout 2000 (`ms_lifetime`); from it on, for every star lighter than `SUPERGIANT_MIN_MASS_SOL`, the end of
the Hertzsprung gap (the base of the giant branch) and of the first ascent (helium ignition) are theirs too
(`giant_branch`, GEN.502, whose durations reproduce the paper's Table 1). A star lighter than `HELIUM_FLASH_MASS_SOL`
then spends `HORIZONTAL_BRANCH_DURATION_GY` on the horizontal branch and `AGB_MS_FRACTION` of its main-sequence
lifetime on the asymptotic giant branch (GEN.286, GEN.287); a star of `HELIUM_FLASH_MASS_SOL` to `SUPERGIANT_MIN_MASS_SOL`
burns helium in its core for as long, as a share of the main-sequence lifetime, as MIST's track does, and then climbs the
asymptotic giant branch the same way; a heavier star has one window, from the end of its main sequence to the end of
its core helium burning, MIST's whole track. Each window is mapped onto its own stretch of MIST's equivalent evolutionary
points, linearly in time: a star halfway through the window is where MIST's track is halfway through the same stretch
(so the two sources need not agree on how long a phase lasts, only on what it is).

The state is a function of the age (a star in a window is at one point of the track), not a draw. A horizontal-branch star
of a star lighter than `HELIUM_FLASH_MASS_SOL` is the exception: its colour and luminosity come from its [Fe/H] and age
(`stellar_evolution.sample_horizontal_branch_colour`), because the track does not hold the spread of a horizontal branch
(it is one mass), and the track gives only the present mass there.
"""

import functools
import math
from collections import namedtuple

import numpy as np

from planetgen import tuning
from planetgen.datasets import evolved_tracks, isochrones
from planetgen.physics import giant_branch, ms_lifetime, remnant_mass

SUBGIANT = "subgiant"
GIANT = "giant"
HORIZONTAL_BRANCH = "horizontal-branch"
CORE_HELIUM_BURNING = "core-helium-burning"
ASYMPTOTIC_GIANT = "asymptotic-giant"
SUPERGIANT = "supergiant"

Window = namedtuple("Window", "name start_gy end_gy f_start f_end first_eep last_eep")
"""A phase: its name, the ages it begins and ends (Gy), the share of the main-sequence lifetime after the main sequence ended
that MIST's track is at each end (`evolved_tracks.Track.frac`) and the equivalent evolutionary points it reads (the stretch
of the track the window is mapped onto)."""

State = namedtuple("State", "window mass_sol luminosity_sol temperature_k radius_sol surface_h1 surface_c12 mist_phase")
"""A star's place on the track: the window it is in, its present mass (Msun, `Timeline.present_mass`), luminosity (Lsun), temperature (K) and
radius (Rsun), its surface hydrogen and carbon mass fractions and MIST's phase number."""

_FIRST_ASCENT_FIRST = evolved_tracks.FIRST_EEP
_FIRST_ASCENT_LAST = evolved_tracks.RGB_TIP_EEP
_HELIUM_FIRST = evolved_tracks.RGB_TIP_EEP
_HELIUM_LAST = evolved_tracks.CORE_HELIUM_END_EEP
_ASYMPTOTIC_FIRST = evolved_tracks.CORE_HELIUM_END_EEP
_ASYMPTOTIC_LAST = evolved_tracks.AGB_END_EEP


def has_core_helium_phases(mass_sol):
    """Whether a star of initial mass `mass_sol` has a horizontal branch (its helium ignites in a flash): below
    `HELIUM_FLASH_MASS_SOL`."""
    return mass_sol < tuning.HELIUM_FLASH_MASS_SOL


class Timeline:
    """The windows of one initial mass and birth [Fe/H], in order, and the state in each (see the module note)."""

    __slots__ = ("mass_sol", "fe_h", "main_sequence_gy", "windows", "_track", "_light")

    def __init__(self, mass_sol, fe_h, main_sequence_gy, windows):
        self.mass_sol, self.fe_h, self.main_sequence_gy, self.windows = mass_sol, fe_h, main_sequence_gy, windows
        self._track = self._light = None

    @property
    def end_gy(self):
        """The age the star's last window ends: it becomes a white dwarf (a lighter star) or collapses (a heavier one)."""
        return self.windows[-1].end_gy

    @property
    def track(self):
        if self._track is None:
            self._track = isochrones.post_main_sequence_track(self.mass_sol, self.fe_h)
        return self._track

    def window(self, name):
        """The window of this name, or None when the star has none."""
        for window in self.windows:
            if window.name == name:
                return window
        return None

    def window_at(self, age_gy):
        """The window `age_gy` is in (a window includes its start), or None before the first or from the end of the last."""
        for window in self.windows:
            if window.start_gy <= age_gy < window.end_gy:
                return window
        return None

    def state(self, age_gy):
        """The `State` at `age_gy`, which must be inside a window."""
        window = self.window_at(age_gy)
        if window is None:
            raise ValueError(f"{age_gy} Gy is outside the post-main-sequence windows of a {self.mass_sol} solar mass star")
        span = window.end_gy - window.start_gy
        return self._state_in(window, min(max((age_gy - window.start_gy) / span, 0.0), 1.0) if span > 0.0 else 0.0)

    def _state_in(self, window, share):
        """The `State` `share` of the way through `window` (0 at its start, 1 at its end)."""
        frac = window.f_start + share * (window.f_end - window.f_start)
        mass, log_l, log_teff, log_r, h1, c12, phase = self.track.at_frac(frac, window.first_eep, window.last_eep)
        frac = min(max(frac, self.track.frac[self.track.index_of(window.first_eep)]),
                   self.track.frac[self.track.index_of(window.last_eep)])
        return State(window.name, self.present_mass(frac, mass, window), 10.0 ** log_l, 10.0 ** log_teff, 10.0 ** log_r, h1, c12, phase)

    def radius_max_sol(self, age_gy):
        """The largest radius, in solar radii, the star has had after its main sequence by `age_gy` (GEN.1009): the largest of the
        track's radius at the end of every window the star has finished and at its place in the window it is in, which is
        where the first ascent and the asymptotic branch end (the red-giant and AGB tips). 0.0 before the first window."""
        largest = 0.0
        for window in self.windows:
            if age_gy <= window.start_gy:
                break
            span = window.end_gy - window.start_gy
            share = min(max((age_gy - window.start_gy) / span, 0.0), 1.0) if span > 0.0 else 1.0
            largest = max(largest, self._state_in(window, share).radius_sol)
        return largest

    def present_mass(self, frac, track_mass, window):
        """
        The mass the star has at `frac` of its main-sequence lifetime after the main sequence ended, where MIST's track has
        `track_mass` (GEN.940). A star of `SUPERGIANT_MIN_MASS_SOL` and more keeps the track's, winds included. A lighter
        one loses, on the first ascent, the track's loss times `RED_GIANT_BRANCH_REIMERS_ETA` over `MIST_REIMERS_ETA`, keeps
        the track's small loss through core helium burning, and on the asymptotic branch goes from the mass it has to the
        white dwarf it makes (`remnant_mass`, so a giant and the dwarf it becomes are one mass) as the track loses its own
        mass; never under that white dwarf's. A star of 6.5 solar masses and more has a track that stops before the thermal
        pulses, so its asymptotic branch goes to the dwarf in a straight line in time (a stand-in).
        """
        initial = self.mass_sol
        if initial >= tuning.SUPERGIANT_MIN_MASS_SOL:
            return track_mass
        track = self.track
        scale = tuning.RED_GIANT_BRANCH_REIMERS_ETA / tuning.MIST_REIMERS_ETA
        tip, helium_end = track.index_of(evolved_tracks.RGB_TIP_EEP), track.index_of(evolved_tracks.CORE_HELIUM_END_EEP)
        dwarf = remnant_mass.white_dwarf_mass_sol(initial)
        if frac <= track.frac[tip]:
            return max(initial - scale * (initial - track_mass), dwarf)
        extra = (scale - 1.0) * (initial - track.star_mass[tip])
        if frac <= track.frac[helium_end]:
            return max(track_mass - extra, dwarf)
        helium_mass = max(track.star_mass[helium_end] - extra, dwarf)
        if not track.reaches(evolved_tracks.AGB_END_EEP):
            span = window.f_end - window.f_start
            share = min(max((frac - window.f_start) / span, 0.0), 1.0) if span > 0.0 else 0.0
            return helium_mass - (helium_mass - dwarf) * share if window.name == ASYMPTOTIC_GIANT else helium_mass
        lost = track.star_mass[helium_end] - track.star_mass[track.index_of(evolved_tracks.AGB_END_EEP)]
        share = 1.0 if lost <= 0.0 else min(max((track.star_mass[helium_end] - track_mass) / lost, 0.0), 1.0)
        return helium_mass - (helium_mass - dwarf) * share

    def bright_pieces(self, min_luminosity_sol, skip=()):
        """
        Every age interval, as `(window name, start_gy, end_gy)`, in which the star is at least `min_luminosity_sol` bright
        on the track (log L linear between the track's samples), except the windows named in `skip` (the horizontal
        branch of a lighter star, whose luminosity is not the track's).
        """
        if self._light is None:
            self._light = isochrones.post_main_sequence_light(self.mass_sol, self.fe_h)
        frac, log_l = self._light
        excess = log_l - math.log10(min_luminosity_sol)
        if excess.max() < 0.0:
            return []
        above = excess >= 0.0
        pieces = []
        # Where, in the share of the main-sequence lifetime, each stretch between two samples is bright: the whole of it when both
        # ends are, the part after the crossing when only the later one is, the part before it when only the earlier one is.
        delta = np.diff(excess)
        cross = frac[:-1] + (-excess[:-1] / np.where(delta == 0.0, np.inf, delta)) * np.diff(frac)
        first = np.where(above[:-1], frac[:-1], cross)
        last = np.where(above[1:], frac[1:], cross)
        bright = above[:-1] | above[1:]
        for window in self.windows:
            if window.name in skip:
                continue
            low, high = evolved_tracks.sample_index(window.first_eep, len(frac)), evolved_tracks.sample_index(window.last_eep, len(frac))
            span = window.f_end - window.f_start
            if span <= 0.0 or high <= low:
                if above[low]:
                    pieces.append((window.name, window.start_gy, window.end_gy))
                continue
            keep = bright[low:high]
            if not keep.any():
                continue
            starts = np.maximum(first[low:high][keep], window.f_start)
            ends = np.minimum(last[low:high][keep], window.f_end)
            inside = ends > starts
            if not inside.any():
                continue
            starts, ends = starts[inside], ends[inside]
            joined = np.concatenate(([True], starts[1:] > ends[:-1] + 1.0e-15))
            starts = starts[joined]
            ends = ends[np.append(joined[1:], True)]
            scale = (window.end_gy - window.start_gy) / span
            for start, end in zip(starts.tolist(), ends.tolist()):
                pieces.append((window.name, window.start_gy + (start - window.f_start) * scale,
                               window.start_gy + (end - window.f_start) * scale))
        return pieces


@functools.lru_cache(maxsize=1 << 17)
def timeline(mass_sol, fe_h=0.0):
    """The `Timeline` of a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` (cached: the sampling tables ask for the
    same masses again and again). Raises `fetch.DataMissing` when the isochrones were not downloaded."""
    landmarks = isochrones.post_main_sequence_grid().landmark_fracs(mass_sol, fe_h)
    if landmarks is None:
        raise ValueError(f"the isochrone grid holds no track after the main sequence for {mass_sol} solar masses")
    tip, helium_end, asymptotic_end, last = landmarks
    t_ms = ms_lifetime.lifetime_gy(mass_sol, fe_h)
    if mass_sol >= tuning.SUPERGIANT_MIN_MASS_SOL:
        windows = (Window(SUPERGIANT, t_ms, t_ms * (1.0 + last), 0.0, last, _FIRST_ASCENT_FIRST, evolved_tracks.LAST_EEP),)
        return Timeline(mass_sol, fe_h, t_ms, windows)

    base_gy = t_ms * giant_branch.subgiant_end_fraction(mass_sol, fe_h)
    ignition_gy = t_ms * giant_branch.giant_end_fraction(mass_sol, fe_h)
    gap_share = (base_gy - t_ms) / (ignition_gy - t_ms)
    windows = [Window(SUBGIANT, t_ms, base_gy, 0.0, gap_share * tip, _FIRST_ASCENT_FIRST, _FIRST_ASCENT_LAST),
               Window(GIANT, base_gy, ignition_gy, gap_share * tip, tip, _FIRST_ASCENT_FIRST, _FIRST_ASCENT_LAST)]
    if has_core_helium_phases(mass_sol):
        horizontal_end = ignition_gy + tuning.HORIZONTAL_BRANCH_DURATION_GY
        asymptotic_gy = t_ms * tuning.AGB_MS_FRACTION
        windows.append(Window(HORIZONTAL_BRANCH, ignition_gy, horizontal_end, tip, helium_end, _HELIUM_FIRST, _HELIUM_LAST))
        windows.append(Window(ASYMPTOTIC_GIANT, horizontal_end, horizontal_end + asymptotic_gy, helium_end,
                              asymptotic_end, _ASYMPTOTIC_FIRST, _ASYMPTOTIC_LAST))
    else:
        helium_gy = t_ms * (helium_end - tip)
        asymptotic_gy = t_ms * (asymptotic_end - helium_end)
        windows.append(Window(CORE_HELIUM_BURNING, ignition_gy, ignition_gy + helium_gy, tip, helium_end, _HELIUM_FIRST,
                              _HELIUM_LAST))
        if asymptotic_gy > 0.0:
            windows.append(Window(ASYMPTOTIC_GIANT, ignition_gy + helium_gy, ignition_gy + helium_gy + asymptotic_gy,
                                  helium_end, asymptotic_end, _ASYMPTOTIC_FIRST, _ASYMPTOTIC_LAST))
    return Timeline(mass_sol, fe_h, t_ms, tuple(windows))
