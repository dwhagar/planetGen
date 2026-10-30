# Stars after the main sequence, from the MIST tracks (GEN.935)

Status: built (Star 1, 2026-10-29 and 30), part of the star mass function epoch (GEN.932); the present mass is GEN.940. Read with the star class audit (`reports/star-class-audit-2026-10-29.md`, findings H1 to H3 and section 6) which found the faults this change removes.

## What was wrong

The model drew a post-main-sequence star's luminosity and temperature from hand-set ranges that had no link to its mass or to each other:

- A supergiant's luminosity was its end-of-main-sequence value times a uniform 1 to 1.5, a third of them red (3,500 to 4,500 K) and the rest blue (9,000 to 35,000 K), whatever the mass. Thirty percent were smaller than the star's own radius leaving the main sequence, a third of those brighter than 300,000 Lsun were red (100 and 150 solar mass "red hypergiants" of 3,000 to 4,500 Rsun against 1,500 Rsun for the largest real ones), and no star was yellow.
- A giant or asymptotic giant had a uniform 3,500 to 5,000 K at any luminosity, so a 1 solar mass tip was at most 54 Rsun where the models give 65 to 123.
- A star from 2.2 solar masses up had no core helium burning, so a 3 solar mass giant was 970 Lsun where the models give 113 to 144, and a fixed tenth of the main-sequence life as a subgiant where the Hurley, Pols and Tout formulae give under 1 percent.

## What it is now

**Data (setup, `datasets/isochrones.py`, `datasets/evolved_tracks.py`).** The MIST v1.2 isochrone archive that setup already downloads (v/vcrit 0.4, solar-scaled, 15 values of [Fe/H], Choi et al. 2016) lists, for each of 107 ages, the stars of every initial mass as rows of equivalent evolutionary points (EEPs): the same integer means the same stage in every star's life (454 the end of the main sequence, 605 the tip of the red giant branch, 707 the end of core helium burning, 808 the start of the thermal pulses and the end of a massive star's track, 1409 the end of the asymptotic giant branch, 1710 a white dwarf). At one EEP the initial mass falls as the age rises, so the age a given mass reaches that EEP is found by interpolating across the isochrone ages in log mass. Doing that for 130 log-spaced masses from 0.8 to 300 solar masses and every sixth EEP (and the landmarks) gives each mass its own track. This resolves the Hertzsprung gap and the blue loops, which the isochrones cannot show as ages because they are 0.05 dex apart: the rows are the EEPs, not the ages. The derived file `post_main_sequence.npz` (5 MB) holds, per mass and EEP, the time since the main sequence ended as a share of its lifetime, the present mass, log L, log Teff, log R, the surface hydrogen and carbon fractions and MIST's phase number. The main-sequence file gains the present mass (used by GEN.940).

**Windows (`physics/post_main_sequence.py`).** The generator's own phase times say when a star is in a phase; MIST says what it looks like there.

| Initial mass | Windows |
|---|---|
| under `HELIUM_FLASH_MASS_SOL` (2.2) | subgiant and giant (the Hurley, Pols and Tout base of the giant branch and helium ignition of `giant_branch`, GEN.502), horizontal branch (`HORIZONTAL_BRANCH_DURATION_GY`, with the colour and luminosity of `stellar_evolution`), asymptotic giant (`AGB_MS_FRACTION` of the main-sequence life) |
| 2.2 to 8 | subgiant and giant (the same formulae), core-helium-burning giant and asymptotic giant, each as long, as a share of the main-sequence lifetime, as MIST's track |
| 8 and up | one supergiant window, MIST's whole track from the end of the main sequence to the end of core helium burning |

Each window is mapped linearly in time onto its own stretch of EEPs, so a star halfway through is where MIST's track is halfway through the same stretch. The state is a function of the age, not a draw, except on the horizontal branch of a light star (the track holds one mass, not the spread of a branch).

**Classes.** Subgiant is IV, the giant, core-helium-burning and asymptotic windows are III (II from `BRIGHT_GIANT_MIN_MASS_SOL`), a supergiant is IB, IAB, IA or 0 by luminosity as before. A giant or supergiant born with at most 13 solar masses, at least 500 Lsun and between 5,000 and 7,000 K is a classical Cepheid (the type text says so); the yellow stars of 20 to 45 solar masses are supergiants and hypergiants, not Cepheids.

**Bright-star tables (`generation/star_population.py`).** `_bright_windows` is the list of ages in which the track is at least as bright as the threshold, found from the track (`Timeline.bright_pieces`); the horizontal branch and a Wolf-Rayet star keep a window with the chance their drawn luminosity clears it. The cell bound takes the earliest start and latest end over five masses across the cell at the threshold less `BRIGHT_TABLE_MAIN_SEQUENCE_SLACK`.

## What it removed

`GIANT_TEMPERATURE_RANGE_K`, `GIANT_BRANCH_TIP_LUMINOSITY_SOL`, `GIANT_BRANCH_TIME_SLOPE`, `BRIGHT_GIANT_LUMINOSITY_RANGE_SOL`, `AGB_LUMINOSITY_RANGE_SOL`, `SUBGIANT_MAX_BRIGHTENING`, `SUBGIANT_END_TEMPERATURE_K`, `SUPERGIANT_LUMINOSITY_GROWTH_RANGE`, `RED_SUPERGIANT_*`, `BLUE_SUPERGIANT_*`, `SUBGIANT_PHASE_END_MS_FRACTION` and `GIANT_PHASE_END_MS_FRACTION`, and the functions that drew from them (`giant_luminosity_law`, `sample_giant_luminosity_sol`, `giant_bright_chance`, `agb_bright_chance`, `subgiant_end_luminosity_sol`, `subgiant_end_fraction`, `giant_end_fraction` of `stellar_evolution`). No shim.

## Checks (all from the audit's section 6, read in the tests)

- 1 solar mass first ascent: 28 Lsun at 4,507 K and 8.7 Rsun; 2 solar masses: 164 Lsun at 4,383 K and 22 Rsun. Tip (Serenelli et al. 2017, Table 2): Z = 0.0001 gives 65 Rsun at 4,550 K, Z = 0.004 gives 123 Rsun at 3,640 K; the tracks are within 12%.
- Core helium burning: 3 solar masses 113 to 144 Lsun, 5 solar masses about 1,500 Lsun, 8 solar masses 6,477 Lsun at 4,097 K.
- No red supergiant above 320,000 Lsun (Humphreys-Davidson limit, Davies, Crowther and Beasor 2018); the tracks themselves stop being red above about 25 solar masses.

## The present mass (GEN.940)

A star's mass is no longer its birth mass. `Timeline.present_mass` and `main_sequence.State.mass_sol` give:

- **Main sequence.** MIST's present mass (`star_mass`, which the setup now keeps for every row): a Sun-like star is its birth mass to a thousandth, a 40 solar mass star has lost 5 percent by the end of its hydrogen burning.
- **Heavier than `SUPERGIANT_MIN_MASS_SOL`.** The track's mass all through: the Dutch wind scheme of MESA (Vink et al. 2000 and 2001, de Jager et al. 1988, Nugis and Lamers 2000) is in the rows. A 20 solar mass star ends its red supergiant life at 13 and a 40 solar mass one at 23.
- **Lighter, first ascent.** MIST took the Reimers (1975) wind with eta_R = 0.1 (Choi et al. 2016, section 3.7.1); McDonald and Zijlstra 2015 fit 0.477 +- 0.070 to the horizontal-branch morphology of 22 globular clusters. The loss of the first ascent is MIST's times 0.477 over 0.1 (`RED_GIANT_BRANCH_REIMERS_ETA`, `MIST_REIMERS_ETA`): a one solar mass star is 0.83 at the tip against MIST's 0.954, its horizontal branch 0.78, and a 0.85 to 0.9 solar mass star of an old cluster has a horizontal branch of 0.58 to 0.64, inside the observed 0.55 to 0.7.
- **Core helium burning:** the track's small loss. **Asymptotic branch:** the mass goes from the horizontal-branch mass to the white dwarf the star becomes (`remnant_mass.white_dwarf_mass_sol`) as the track loses its own, so a giant and the dwarf it turns into are one mass. A star of 6.5 solar masses and more has a track that stops before the thermal pulses, so its asymptotic branch goes to the dwarf in a straight line in time (a stand-in, PS7).
- A star lighter than the lightest track of the grid takes that track and its present mass in proportion.

The star's record carries it as `mass_kg`, so the surface gravity, the orbital periods and the habitable-zone and Hill-sphere sizes of its planets follow. The hot-subdwarf mass of the audit's L8 (0.47 solar masses) and the Wolf-Rayet mass (GEN.939) are not in it.

## Wolf-Rayet stars (GEN.939)

In the Wolf-Rayet window a star's present mass is the track's (GEN.940, 13 to 25 solar masses at the Sun's metallicity for initial masses of 30 to 150) and its luminosity follows that mass by the helium-star relation of MIST's own Wolf-Rayet rows (phase 9, 24,610 samples over every [Fe/H]): log L/Lsun = 4.31 + 1.16 log M/Msun, 0.064 dex scatter, fitted when the tracks are read (`isochrones.helium_star_relation`). The luminosity is drawn from that normal law in log luminosity and cut at `WR_LUMINOSITY_RANGE_SOL` (1e5 to 3e6; it was 1e5 to 1e6 log-uniform for every star), and the chance a star clears a bright-star floor follows its mass (`wolf_rayet.bright_chance`). The temperature and subtype are as before. The audit's finding that the code was too low above 60 solar masses does not show in MIST's tracks: their Wolf-Rayet stars of 60 to 150 solar masses end at 4e5 to 9e5 Lsun because the winds have taken them down to 13 to 23 solar masses.

## Known limits

- The tracks of 30 to 45 solar masses spend about half their post-main-sequence life as yellow hypergiants (5,000 to 9,000 K, 0.3 to 1 million Lsun). MIST has no eruptive mass loss, so these stars do not become luminous blue variables or Wolf-Rayet stars as the real ones do. They are rare in the sky (the mass range is a few percent of the supergiants) and the rows are tier T1 data used as the table gives them.
- The Wolf-Rayet window is still the last 0.3 Myr of the supergiant window (GEN.243) and the threshold mass the old one; MIST's own Wolf-Rayet stars of 100 solar masses and more are Wolf-Rayet from the main sequence (WNh), which the generator does not model.
- At the 2.2 solar mass boundary the horizontal branch of the lighter star (a fixed 100 Myr) and the core helium burning of the heavier (MIST, about 17 percent of the main-sequence life, 130 Myr at 2.2) differ, so the end of the last giant phase steps by 10 percent of the main-sequence life there.
- A mass below the lightest track of the grid that reaches the thermal pulses (0.85 to 1.0 solar masses, by metallicity) takes that track.
