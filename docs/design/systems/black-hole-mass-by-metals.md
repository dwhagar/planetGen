# The mass of a stellar-mass black hole follows the metals (GEN.554)

**Status:** built (stand-in BHM1, `docs/design/research/gap-filling-plan.md` section 6.337). Code: `physics/fryer_remnants.py`,
`physics/remnant_mass.py` (`black_hole_stretch`, `stretched_mass_sol`, `draw_mass_sol(..., stretch=)`), `galaxy/remnant_density.py`
(`stretch_array`, `place_stretch`, `stretch_from_peaks`), the scatter (`generation/phenomenon_scatter.py`), the sector fill
(`generation/run_sector.py`), the member tree (`generation/multistar.py`) and the remnant composition (`physics/remnant_composition.py`).

## What was wrong

A stellar-mass black hole's mass was uniform over 5 to 20 solar masses everywhere, whatever the metals of the stars that made it. The
failed share already follows the metals (`black_hole_formation.failed_share`, GEN.563), and so does the count, but the mass did not:
the holes of a metal-poor halo were as heavy as the disk's, where [R1466] and the gravitational-wave holes of metal-poor places say
they are heavier (a lighter wind leaves more of the star).

## What was tried first, and why it is not used

[R1466]'s delayed engine (their Eqs. 5, 7 and 9, baryonic mass, 0.9 of it for a hole, their Eq. 14) gives a hole mass from a progenitor's
initial mass and metals. Over the stars the failure windows hold (`remnant_progenitors.kind_windows`) it gives a mean of 4.4 solar masses
at the sun's metals, 14 at half of them, 21 at a third and 31 at a tenth. At the sun, 67% of the stars of the top window (31 to 150
solar masses) come out under 5 solar masses (their Eq. 9, wind-stripped), which no observed hole is, so the kind (the windows) and the
mass (the fits) disagree there. The neutron stars the fits give (1.17 to 1.35) are under the verified observed Gaussian
(1.49 +- 0.19, [R3396]). Verified wins, so neither is taken as it stands: the neutron stars keep the Gaussian and the black holes keep the
registered law as the anchor.

## The rule

* The anchor is the registered law, uniform over 5 to 20 solar masses (`tuning.BLACK_HOLE_MASS_RANGE_SOLAR`), the widest published range
  ([R3397]'s 7.8 +- 1.2 are 16 holes of transient X-ray binaries chosen by flux; [R1328] has holes to 28; [R3348] 14.2 and 11.7).
  It stands for the holes of the disk around the sun.
* Fryer's delayed engine sets only the trend: `fryer_remnants.mass_trend(fe_h)`, the IMF-weighted mean hole mass over the failure
  windows at the failed share of that [Fe/H] (`failed_share` of the reference structure) and the metals `Z = 10^[Fe/H]` (held at 1
  from the sun up), over its mean for the spread of [Fe/H] of the disk at the reference radius (`intermediate`, 0 +- 0.15). It is
  1 for that disk, 0.73 at and above the sun's metals, 1.75 at -0.2 dex, 3.2 at -0.5 and 6.0 with no metals.
* `remnant_mass.black_hole_stretch(fe_h)` multiplies the span above the 5 solar mass edge by the trend, held so that the top of the
  span is at most `tuning.BLACK_HOLE_MASS_TOP_SOLAR` = 50 solar masses. The top is inside the lower edge of the pair-instability
  mass gap, 45 solar masses, 40 to 56 for the carbon-oxygen reaction rate's one-sigma range ([R5365]): the pulsations that gap comes
  from strip the heavier stars before they collapse. The hole's mass is `5 + (anchor mass - 5) x stretch`: the same rank in the law
  at every place, and a flat law up to the top, no pile at the top.
* A star of the member tree and a remnant's composition take the stretch of their own [Fe/H]. The scatter's rows (the disk, bulge and halo scatter, the reach groups and the rows of a component's own pass), the rows of a ring
  sector and the holes of the sector fill take the stretch of their place: the densities' mean of each population's mean stretch
  over its [Fe/H] distribution (`remnant_density.stretch_array`, the same tables as the failed share's mean).
* The scatter's mass cut (`tuning.PHENOMENON_MIN_MASS_SOLAR`, 14) and `remnant_mass.share_above`, `moments` stay in the anchor law's
  masses, so the counts above and below the cut, and every count, do not move: a hole of a metal-poor place above the cut is above it
  by rank, and its mass is larger.

## What moves (Milky Way-like shape, `reports/gen554-mass-from-progenitor-2026-11-24/`)

The mean hole mass of the galaxy goes from 12.5 to 17.7 solar masses; the share at or above 10 solar masses from 67% to 79%, 20 from 0 to
38%, 30 from 0 to 6% and 40 from 0 to 1.7%; the anchor's cut at 14 solar masses holds 40% of the holes by rank and 63% by mass. The
galaxy as a whole is poorer in metals than the disk around the sun (the old disk is at -0.325, the halo at -1.6), so its holes are
heavier than the anchor's.

## What is left as a stand-in

* The standalone remnants (scatter rows, field hosts) draw their mass before any progenitor, so the mass does not follow a single progenitor;
  only the place's trend reaches them. A mass from the progenitor itself waits for the progenitor to be drawn first (the age, `remnant_progenitors.collapse_age_gy`, draws one
  after the mass today).
* The field hosts (`generation/field_host.py`) have no metallicity at hand and draw at the anchor.
* The trend uses the reference structure's failed shares, not the galaxy's own, as the windows do (`reference_failed_share`).
