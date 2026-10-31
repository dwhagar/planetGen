# Carbon stars and S stars (GEN.1026)

Code: `physics/agb_carbon.py`. Sheet: the derived "Carbon or S star" row of a star (`describe.sheet.CARBON_KEY`), in the identity group. Tests: `tests/test_agb_carbon.py`, `tests/test_property_sheet.py`. Registry 6.139 (CS1 to CS11).

Nothing is stored and nothing in generation changes: the sheet reads the star's stored state, lifespan (the end of its last giant phase, from which a bisection on the track gives the initial mass), mass, [Fe/H] and age. The post-main-sequence track has no thermal pulses (its asymptotic branch is one stretch), so the module counts them from the star's place in that stretch and lets the third dredge-up bring carbon up, as the sheet is made.

## The model

1. **Pulses.** A star of initial mass `M` runs `pulses_in_branch(M)` thermal pulses over the asymptotic window (17 at 1.5 solar masses to 25 at 2.75, Karakas and Lugaro 2016 Table 1), and has run that number times the share of the window it has used.
2. **Dredge-up.** A star at or above `tdu_minimum_mass_sol` (1.5 solar masses at solar [Fe/H], 0.2 lower per dex of [Fe/H] below it, never under 1.0) dredges up at pulse `n` the share `lambda(n) = lambda_max(M) (1 - exp(-n/4))` of the core growth of the pulse (0.004 solar masses); 22 percent of what is dredged is carbon.
3. **Surface C/O.** The dredged carbon over the oxygen of the envelope, the envelope being the mean of its mass at the start of the branch and now, less a 0.6 solar-mass core. The envelope's own C/O is the A21 photospheric one (0.59 at solar).
4. **Hot bottom burning.** A star of 4.0 solar masses or more burns the carbon again and is never one.

A star on the asymptotic branch with C/O at or above 1 is a **carbon star**; with C/O from 0.9 to 1 and some dredge-up, an **S star**. Because C/O climbs through the S band on its way past 1, an S star is a carbon star in the making: the row reads "S star" for a stretch of the branch and then "carbon star".

## What it gives

At solar [Fe/H] the carbon stars are the stars of 1.5 to 4 solar masses in the last half of their branch; a metal-poor star dredges at lower mass and turns carbon-rich after fewer pulses (the same carbon is a larger share of less oxygen), as the carbon-star counts of the metal-poor dwarf galaxies show. The C/O of a very metal-poor star comes out in the tens: it is the ratio the model's arithmetic gives, only the flag is used.

## Not here

- The oxygen-rich Mira-like stars of the branch's start are not marked (the Variability row handles the temperature).
- The track's [Fe/H] of a star with an alpha-enhanced branch (GEN.1024, PR #4397) is not used: the module reads the birth [Fe/H]; when #4397 lands the sheet's call passes `alpha_tracks.track_fe_h` instead.
- Hot bottom burning does not make a nitrogen-rich flag, and the s-process products are not tracked: an S star is named by its C/O alone.
