# Alpha-enhanced stars and the age-metallicity relation (GEN.744, lock-free step)

Status: modules built, nothing in generation calls them (wiring is GEN.744's second step and rides train 3 with the star epoch). Code: `physics/alpha_tracks.py`, `physics/age_metallicity.py`. Tests: `tests/test_alpha_and_age_metallicity.py`. Registry 6.34 (AL1, AR1 to AR3).

## Alpha-enhanced tracks

The MIST grid is solar-scaled. A star of the high-alpha branch (`stellar_abundances.BRANCH_HIGH_ALPHA`) behaves like a solar-scaled star of higher [Fe/H] (Salaris, Chieffi and Straniero 1993): `[M/H] = [Fe/H] + log10(0.638 f + 0.362)`, `f = 10 ** [alpha/Fe]`. `alpha_tracks.equivalent_fe_h(fe_h, branch)` gives the [Fe/H] the isochrone lookup takes: +0.27 dex for a high-alpha star at [Fe/H] -0.8, about 0.02 for a thin-disc star at solar. The star's [alpha/Fe] is the mean of its branch's [O/Fe], [Mg/Fe], [Si/Fe] and [Ca/Fe] (the hinge fits `composition` already uses), so the track and the composition of one star agree.

Wiring: the lookup needs the branch before the track is read, and the branch is today decided after [Fe/H] and age in the composition step (`stellar_abundances.branch_from_uniform`); the wiring moves the branch draw ahead of `evolve_star` and passes the equivalent [Fe/H] to `isochrones.main_sequence_point` and the lifetime. Both belong to Bug 1's and Star 1's star-epoch files (`stellar_evolution.py`), so the wiring is coordinated with them through TODO.

## Age-metallicity relation

Haywood et al. 2013 (Section 3.3, Fig. 9): the thick disc has a tight relation rising 0.15 dex per Gyr over its 4 to 5 Gyr of formation; the thin disc gains 0.025 dex per Gyr. `age_metallicity.mean_fe_h(age, branch)` is a line at each slope through the anchors AR1 (thin, +0.05 at birth) and AR2 (thick, -0.10 at 8 Gyr, held to 13 Gyr). `sample_fe_h(age, branch, normal)` draws about it with scatters 0.15 (thin) and 0.10 (thick).

It agrees with today's populations (young and intermediate disc 0, old disc -0.2) and Bensby et al. 2014's thin-disc floor near -0.7 sits at the thin disc's 3 sigma.

Wiring: the age is drawn first and [Fe/H] follows from it and the branch (it reverses the order of GEN.553 for the disc, and the thick disc's window moves from the 3 to 10 Gyr "old" population to 8 to 13 Gyr: GEN.944 part 1, the "10 to 13 Gyr thick-disk window"). It changes every disc star, so `--affects population` or `fill` and a census re-baseline.

## Not here

- Variable-star flags (Cepheids, Mira, carbon and S stars) of GEN.944: they need Star 1's post-main-sequence phases (GEN.935).
- Alpha-enhanced tracks themselves (Dartmouth, BaSTI-alpha) as a download: the Salaris shift is the accepted secondary-source stand-in until the equivalence is shown to fail.
