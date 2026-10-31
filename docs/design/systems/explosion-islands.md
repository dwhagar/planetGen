# Explosion islands: which massive stars leave a black hole (GEN.554, slice 2)

Status: built (GEN.554 slice 2). The kick step is the next slice. Code: `physics/remnant_progenitors.py`, `tuning.EXPLOSION_FAILURE_ISLANDS_SOL`, the progenitor draws of `generation/multistar.py` and `physics/remnant_progenitors.collapse_age_gy`, and the remnant moments of `generation/place_terms.py`.

## What changed

Before: a neutron star's progenitor was drawn over 8 to 30 solar masses and a black hole's over 20 to 100 (`NEUTRON_STAR_PROGENITOR_SOL`, `BLACK_HOLE_PROGENITOR_SOL`, "a design default"), the two ranges overlapped, and a black hole was a direct collapse when its progenitor was 40 solar masses or more. None of it knew the failed share of `black-hole-density.md`, so the neutron star and black hole counts and the progenitors behind them could disagree.

After: the star's initial mass decides, in windows. A population that fails a share `F` of its 8 to 150 solar mass stars (`black_hole_formation.failed_share`) fails them in the initial masses of `failure_windows(F)`, and the neutron stars are the massive stars outside them. A progenitor is drawn from the IMF over its kind's windows (`sample_mass_sol`), and a black hole whose progenitor lies in a window failed outright, which is a direct collapse (`multistar.DIRECT_COLLAPSE_PROGENITOR_SOL` is gone). The windows hold exactly `F` of the massive stars (a test), so the counts of `neutron-star-density.md` and the progenitors agree by construction.

## The windows

Only two windows are sourced: O'Connor and Ott 2011 give, for the Woosley and Heger 2007 progenitors at solar metallicity, black holes "only in ~23-25 Msun and ~35-45 Msun" (0 to 7% of the massive stars), the "islands of 'explodability' in a sea of black hole formation" of Sukhbold et al. 2016 (`black-hole-formation-sources.md`, sections 2.1 and 2.2). The generator's IMF puts 6.8% of the massive stars in those two windows, which is the top of O'Connor and Ott's solar band. How the windows grow with the share is not sourced and is a stand-in:

| Share `F` | Windows |
|---|---|
| 0 | none: every massive star leaves a neutron star |
| up to 6.8% | the two islands, shrunk about their centres in log mass until they hold `F` |
| above 6.8% | the two islands and the top of the mass range down to the mass that makes up `F` (9%: 89.7 to 150; 20%: 28.5 to 150 with 23 to 25; 32%: 18.6 to 150) |
| 45% (the widest verified share) | 14.5 to 150, where Ertl et al. 2016) put the lowest island at about 15 |

The windows nest: a larger share fails every mass a smaller one does (a test).

## Which share

A draw with a galaxy's structure around it could take that galaxy's share at its place, but the companion and age draws (`multistar._remnant`, `collapse_age_gy`) have no structure. They take the share of the middle draw at solar metallicity (`reference_failed_share`, 13.6%), the same for every galaxy; a caller that knows the place's share passes it (`failed_share`). The count of black holes and neutron stars at a place still follows the galaxy's own draws (slice 1).

## Stand-ins and open items

- **[stand-in]** The nesting of the windows above the islands (the top of the range first) is not sourced; it only has to hold `F`, and it makes the heaviest stars fail first, as at low metallicity (O'Connor and Ott: black holes for "MZAMS >~ 30 Msun").
- The direct-collapse flag of a black hole with no progenitor (a scatter row) keeps the hole mass rule (`DIRECT_COLLAPSE_HOLE_SOL`, 20 solar masses, **[stand-in]**).
- **[stand-in]** Remnant masses are unchanged: the black hole's mass is still uniform over 5 to 20 and the neutron star's over 1.1 to 2.2 solar masses. A mass that follows the progenitor (helium core 7.7 to 9.2, Sukhbold) is a later slice with the mass law (card 6 of `black-hole-density.md`).
- The magnetar rate keeps the tuning rate as a **[stand-in]** (GEN.554's later slice). The empirical placement of `remnant_distribution` is unchanged. Rotation, binaries and the metallicity shift of the windows (Maltsev et al. 2025's carbon-oxygen core windows) are not modelled.

## Tests

`tests/test_remnant_progenitors.py`: the windows hold the share across 0 to 100%, are ordered and inside 8 to 150, nest, equal the two islands at their own share and reach 14.5 solar masses at 45%; the two kinds share the range; progenitors lie in their windows and above the lowest mass asked, follow the IMF over the windows and repeat for a seed; a black hole of a failed progenitor is a direct collapse; the multi-star builder's remnants come from the windows.
