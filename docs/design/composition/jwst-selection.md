# The JWST atmosphere target-selection model (GEN.506)

Status: built 2026-10-23 by lane M. Code `src/planetgen/physics/jwst_selection.py`, fit script `scripts/fit_jwst_selection.py`, data `src/planetgen/physics/data/jwst_selection.json`, tests `src/tests/test_jwst_selection.py`. It is the selection model `jwst-target-selection` of the tuning harness (GEN.503, `docs/design/composition/tunable-gaps.md` section 4) for retention and lava worlds (B3, GEN.511), giants (B8, GEN.516) and aerosols (B11, GEN.519).

## 1. What it answers

JWST observes planets chosen for signal: nearby, hot, large scale height, bright host. A tuned fraction of clear, retained or hazy atmospheres must be compared with the planets JWST has observed and not with the whole drawn population.

## 2. The model

A logistic regression of the label "has a published JWST transmission or emission spectrum" (the NASA Exoplanet Archive `spectra` table; direct images are not counted) on three things:

- log10 of the transmission spectroscopy metric `TSM = S Rp^3 Teq / (Mp R*^2) 10^(-mJ/5)` ([R2315], equation 1; `S` = 0.190, 1.26, 1.28, 1.15 for under 1.5, 1.5 to 2.75, 2.75 to 4 and over 4 Earth radii; [R2315] read);
- log10 of the emission metric `ESM = 4.29e6 B7.5(1.1 Teq) / B7.5(T*) (Rp/R*)^2 10^(-mK/5)` (equation 4);
- an offset for each of four radius classes (under 1.5, 1.5 to 4, 4 to 10, over 10 Earth radii).

The sample is the 4,618 archive planets that transit and have a radius, a mass (measured or from the archive's mass-radius relation), an equilibrium temperature, the host's radius and effective temperature and the host's J and K magnitudes; 78 have a JWST spectrum. The fit is by Newton's method (`scripts/fit_jwst_selection.py`, retrieved 2026-10-23). Only the fitted numbers are in the repository.

| Coefficient | Value |
|---|---|
| intercept (class under 1.5 Earth radii) | -6.39 |
| per dex of TSM | 3.88 |
| per dex of ESM | 0.31 |
| class offsets (1.5 to 4, 4 to 10, over 10) | -4.32, -5.16, -5.44 |
| area under the ROC curve | 0.937 |

The generator has no apparent magnitude, so `with_magnitudes` draws J from a normal truncated to the archive's range (mean 11.9, sd 1.94) and J-K from a line in log10 of the host's effective temperature (intercept 9.13, slope -2.32, scatter 0.065). `observed_sampler` adds them to a sampler's items; `probability(item)` is the chance the planet is on the list; `jwst_target` is the selection model.

## 3. What it shows (COMPUTED)

- The slope in TSM is steep, 3.9 per dex: the programs were chosen from a ranking by this metric, as Kempton et al. suggested. The emission metric adds little once the transmission metric and the class are known (0.31 per dex).
- The class offsets are what the programs chose, not physics: rocky planets (under 1.5 Earth radii) are selected 22 times from 906 planets, sub-Neptunes 19 from 2,525 (TSM rankings favour them, but few are bright enough), 4 to 10 Earth radii 10 from 442, and giants 27 from 745. The fit reproduces each count exactly by construction (`calibration` in the file, tested).
- Without the class offsets (one radius term) the fit under-counts rocky planets (9 against 22) and over-counts sub-Neptunes (36 against 19); that is why the offsets are in.
- Examples: a hot Jupiter (14 Earth radii, 1100 K, J 10.2) has a probability of 0.28; a GJ 1214 b-like sub-Neptune (2.7 Earth radii, 596 K, M dwarf, J 9.75) 0.53; an Earth around a Sun-like star at J 12, below 0.001.

## 4. What it does not do

- It conditions on a known transiting planet with a mass. A planet that was never found is not in the fit, so the host magnitudes are those of found planets; discovery and the choice of a JWST program are not separated.
- A planet is "observed" when its spectrum is published. Planets observed and not yet published, and programs not yet flown, are not in the list: the model is a statement about the published list on the retrieval date and changes when the archive grows.
- The class offsets and the steep slope are the choices of the time allocation committee for the first cycles; they are not physical constants.
- A host without magnitudes gets a J from the found-planet distribution; a faint M dwarf at J 15 can be a valid host that no program would observe.
- Validity: transiting planets, J between 3.98 and 16.2, the four radius classes.
