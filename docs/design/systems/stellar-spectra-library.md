# Stellar spectra library: model atmospheres plus the measured chromosphere (GEN.970)

Status: built (GEN.970, parent GEN.745) and wired (GEN.971, rides the star-epoch merge train 3). Code: `datasets/stellar_spectra.py` (download and reduction), `physics/stellar_library.py` (the lookup and the chromosphere). Tests: `tests/test_stellar_spectra.py` (blackbody fixtures), `tests/test_stellar_library_real.py` (the real grids against the Sun and MUSCLES).

## What it replaces

`physics/stellar_spectrum.py` (GEN.443) gives a star a blackbody of its effective temperature with a seeded UV-C share (0.1 to 1, log-uniform) and a UV-B and UV-A share on one family fitted to the Sun. The blackbody overpredicts the Sun's ultraviolet 3.9 times in UV-C, and the seeded share has no dependence on temperature or activity. The map colour of a star (`physics/photometry.py`) is a blackbody above 4500 K and a small BT-Settl offset table below.

## Sources

- ATLAS9 fluxes of Castelli and Kurucz 2003 (`fp00k2odfnew`, solar metallicity), 476 models, 3500 to 50000 K, log g 0 to 5, from Fiorella Castelli's page at INAF Trieste. This is the hot half (3500 K up).
- PHOENIX-ACES-AGSS-COND-2011 of Husser et al. 2013, from the Goettingen server (ftp.gwdg.de), 2300 to 3500 K in 200 K steps at log g 4.5, 5.0 and 5.5: the M dwarfs.
- The MUSCLES Treasury Survey (France et al. 2016, ApJ 820, 89), eleven M and K stars, for what the photosphere lacks. BT-Settl itself (phoenix.ens-lyon.fr) resets the connection from the cloud machines and the MUSCLES HLSP files (archive.stsci.edu) fail their TLS chain there, so PHOENIX-ACES stands for BT-Settl (same code, different dust treatment, the same family) and the MUSCLES numbers are the paper's Tables 1 and 3. Neither is a loss for the question asked: the tests below hold the result to the measured stars.

Setup downloads 130 MB (transient pieces, sha256 and size pinned), reduces each model to the share of its bolometric power in 1500 log-spaced bins from 9 nm to 100 micron (`stellar_spectra.npz`, 2.5 MB), and deletes the sources. Dataset id SL1.

## How it works

`stellar_spectra.band_fractions(T_eff, log g, bands)` interpolates the shares of the models around a star linearly in log T_eff and log g, in the logarithm of each bin's share (an ultraviolet tail that falls by a decade between two models is a straight step), and sums them over a band. `stellar_library.band_shares(T_eff, log g, xuv_over_lbol)` gives the generator's bands (`stellar_spectrum.BAND_EDGES_NM`) and adds the corona: far-ultraviolet (91.2 to 170 nm) power equal to 1.12 times the star's XUV power (5 to 91.2 nm), which `physics/activity.py` already computes as `l_xuv_w`, split between UV-C and the ionising band. `shares_of_blackbody` returns the three UV shares as multiples of the blackbody's: the quantities `stellar_spectrum.uv_shares` holds, so GEN.971 swaps one for the other.

## What it was checked against

| Check | Library | Measured |
|---|---|---|
| Sun, UV-C / UV-B / UV-A as a share of the blackbody's | 0.29 / 0.71 / 0.77 | 0.256 / 0.625 / 0.787 (ATLAS-1, Thuillier et al. 2004) |
| Near-ultraviolet (170 to 320 nm) of the eleven MUSCLES stars | within 0.33 dex for all eleven | France et al. 2016 Table 3 |
| Far-ultraviolet (91.2 to 170 nm) of the seven M dwarfs | photosphere 1e-8 to 1e-15 of bolometric, measured 1e-4.2 to 1e-4.8 | the chromosphere, not the model, supplies it |
| Far-ultraviolet with the corona (XUV times 1.12) | within 0.55 dex for all eleven | France et al. 2016 Table 3 |
| Closure (integrated flux over sigma T^4) | 0.97 to 1.00 | |

## Wired (GEN.971)

- `physics.activity.generate_activity` stores `uv_c_share` as the library's UV-C power over the blackbody's (`activity.uv_shares_of` -> `stellar_library.uv_shares`): the star's effective temperature, surface gravity (mass and radius) and its corona's XUV over luminosity (`log_lx_lbol`; the photosphere's own XUV is in the model, so only the corona is added). The seeded draw is gone except for a run that disregards a missing library (`fetch.stand_ins_used`). The value is a multiple of the blackbody's and exceeds 1 for a cool active star (an M5 dwarf at saturation: 100 to 1000, the chromosphere); the column unit is `uv-mult`, the sheet row "UV-C power against its blackbody's". A star outside 2300 to 50000 K takes 1. A white dwarf keeps none.
- `physics.body_radiation` takes all three shares (`uv_c`, `uv_b`, `uv_a`) from the library instead of the family built on `uv_c_share` (`stellar_spectrum.irradiance_w_m2(..., shares=)`, `atmospheric_radiation.surface_radiation(..., uv_shares=)`); a star that lacks mass, radius, temperature or luminosity keeps the family from its stored `uv_c_share`. Limitation: the extra UV-C of a cool star is placed in the 100 to 280 nm cells in proportion to the blackbody's shape there, so it sits near 280 nm, where the chromosphere's is really at 120 to 170 nm; the band total is right and the in-band shape is not (a later item: the spectrum on the radiation grid from the library bins).
- The colour of a star's disc (`bodylook.star_surface_color`) and the light of the sky colour (`atmosphere_colour.body_light` into `atmosphere_colour(light=)`, and the giant's) come from the library spectrum by `stellar_library.visible_spectrum` and `star_colour`; outside the library, or for a run that disregards a missing library, the blackbody with the BT-Settl offsets of `photometry.py` (still what a bare `atmosphere_colour(temperature_k=)` call takes). A cool star's seeded spread is no longer used with a library spectrum. Against the blackbody at the Sun's temperature the chromaticity is within 0.006; at 7500 K the library is bluer (line blanketing moves the ultraviolet's flux into the blue), and at 3300 K the red dwarf's sky is whiter.
- The grids are solar metallicity; a metal-poor star's ultraviolet differs. SL4 records it.
- Stars above 50000 K or below 2300 K get no spectrum (None); callers take the blackbody.
- This moves the UV of every star and the colour of every star and sky, so it is `--affects fill` (stars, and the body colour stage) and rides train 3.
