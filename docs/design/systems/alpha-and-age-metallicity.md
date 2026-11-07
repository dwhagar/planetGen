# Alpha-enhanced stars and the age-metallicity relation (GEN.744, lock-free step)

Status: built and wired. Age follows the star's [Fe/H] inside the disc's populations; every track lookup reads the alpha-enhanced equivalent [Fe/H] of the star's branch (GEN.1024); the thick disc spans 8 to 13 Gyr (GEN.1025). Code: `physics/alpha_tracks.py`, `physics/age_metallicity.py`. Tests: `tests/test_alpha_and_age_metallicity.py`. Registry 6.34 (AL1, AR1 to AR4).

## Alpha-enhanced tracks

The MIST grid is solar-scaled. A star of the high-alpha branch (`stellar_abundances.BRANCH_HIGH_ALPHA`) behaves like a solar-scaled star of higher [Fe/H] [R3878]: `[M/H] = [Fe/H] + log10(0.638 f + 0.362)`, `f = 10 ** [alpha/Fe]`. `alpha_tracks.equivalent_fe_h(fe_h, branch)` gives the [Fe/H] the isochrone lookup takes: +0.27 dex for a high-alpha star at [Fe/H] -0.8, about 0.02 for a thin-disc star at solar. The star's [alpha/Fe] is the mean of its branch's [O/Fe], [Mg/Fe], [Si/Fe] and [Ca/Fe] (the hinge fits `composition` already uses), so the track and the composition of one star agree.

Wiring (as built, GEN.1024): the branch is a star's own uniform number (`stellar_abundances.branch_of`, from its [Fe/H] and age), so the track (`alpha_tracks.track_fe_h`), the composition and the radiogenic offsets agree. `evolve_star`, `post_main_sequence.timeline`, `ms_lifetime` and the scatter's windows (`alpha_tracks.track_segments`, `possible_segments`) read the equivalent [Fe/H]. Halo, bulge and cluster stars are older than the thick-disc gate, so their branch is the logistic in [Fe/H]: a metal-poor halo star reads the high-alpha track.

## Age-metallicity relation

[R1863]: the thick disc has a tight relation rising 0.15 dex per Gyr over its 4 to 5 Gyr of formation; the thin disc gains 0.025 dex per Gyr. `age_metallicity.mean_fe_h(age, branch)` is a line at each slope through the anchors AR1 (thin, +0.05 at birth) and AR2 (thick, -0.10 at 8 Gyr, held to 13 Gyr). `sample_fe_h(age, branch, normal)` draws about it with scatters 0.15 (thin) and 0.10 (thick).

It agrees with today's populations (young and intermediate disc 0, old disc -0.2) and [R0349]'s thin-disc floor near -0.7 sits at the thin disc's 3 sigma.

Wiring (as built): the order of the draws is kept, [Fe/H] first and the age after it. The scatter's tables, the places' own [Fe/H] distributions and the bins of GEN.651 all stand on [Fe/H] being drawn from the place before the star is, so reversing the order would rebuild them. Instead the age is chosen to suit the [Fe/H]: `age_metallicity.tilted_age` draws eight ages the way the caller draws one (the population's window and age law, or the track's bright window) and keeps one in proportion to the density of the relation at the star's [Fe/H] (`log_density`: the thin disc's normal, with the thick disc's beside it past 8 Gyr, `thick_share` = AR4). The callers are `sample_living_star` (both of its age draws), `star_population`'s bright-star and mass-band draws and the pinned-mass path of `Star`. A metal-poor star of the old disc comes out older than a metal-rich one (a test fixes the sign and size); the population's [Fe/H] and its age window are unchanged. The young population's window is 0.1 Gyr wide, so the relation is flat across it and nothing moves. Halo, bulge, cluster and satellite stars have histories of their own and are not tilted; nor is a star of no population (the whole disc, 0 to 13 Gyr, drawn with the intermediate disc's [Fe/H]), whose metals are not the relation's mix over those ages and which a tilt would crowd toward the young. Over each of the three disc populations the tilt moves the mean age by less than 0.05 Gyr (a test fixes it).

It changes the age of every disc star, so it needs `--affects population` and `fill` and a census re-baseline (ages shift inside each window, the counts by mass do not).

The 10 to 13 Gyr thick-disc window (GEN.944 part 1) was built as GEN.1025: the thin disc's own history ends at 10 Gyr (`tuning.THIN_DISK_OLD_WINDOW_GY`) and the old population reaches 13 (`galaxy.thick_ages`).

## Not here

- Variable-star flags (Cepheids, Mira, carbon and S stars) of GEN.944: they need Star 1's post-main-sequence phases (GEN.935).
- The alpha shift in the tracks was built as GEN.1024 (see Wiring above).
- Alpha-enhanced tracks themselves (Dartmouth, BaSTI-alpha) as a download: the Salaris shift is the accepted secondary-source stand-in until the equivalence is shown to fail.
