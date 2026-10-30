# Stellar spectra library: model atmospheres plus the measured chromosphere (GEN.970)

Status: lock-free step built (GEN.970, parent GEN.745). Nothing in generation calls it yet; the wiring that replaces the seeded UV-C share and the blackbody map colour is GEN.971 and rides the star-epoch merge train. Code: `datasets/stellar_spectra.py` (download and reduction), `physics/stellar_library.py` (the lookup and the chromosphere). Tests: `tests/test_stellar_spectra.py` (blackbody fixtures), `tests/test_stellar_library_real.py` (the real grids against the Sun and MUSCLES).

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

## Open for GEN.971

- Replace `activity.generate_activity`'s `uv_c_share` draw with `shares_of_blackbody(...)["uv_c"]` from the star's own T_eff, log g and `l_xuv_w / luminosity` (the draw goes away, the column stays; UV-B and UV-A shares come from the library too). This moves the UV of every star, so it rides the star-epoch train.
- Replace the BT-Settl offset table in `photometry.py` by the CIE chromaticity of the library spectrum (the library keeps the shares, not the spectra, so a star's colour needs the bins' wavelengths and shares, which it has: the bin centres and shares are enough for a 1500-point integral against the CIE tables).
- The grids are solar metallicity; a metal-poor star's ultraviolet differs. SL4 records it.
- Stars above 50000 K or below 2300 K get no answer (None); the caller takes the blackbody.
