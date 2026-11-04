# A white dwarf as a star: how hot, how bright, how large, how crystallised (GEN.742)

**Status (2026-10-24):** built (`datasets/white_dwarf_models.py` `cooling_state`, `physics/stellar_evolution.py` `_white_dwarf_state`,
`tests/test_white_dwarf_cooling.py`, `tests/test_white_dwarf_cooling_real.py`). Row S3 of `reports/physics-gaps-2026-10-24.md`.

A white dwarf has no fusion and only cools. Its luminosity used to be a power law of the cooling age (`L = 1e-3 (M / 0.6) t^-1.4`,
clamped to 1e-5 to 50 solar luminosities, the Mestel law) and its radius one zero-temperature mass-radius relation, so a star of a given
mass had one radius whatever its age, a hydrogen layer made no difference, and the core never crystallised. Now all four come from the
cooling sequences the setup step already downloads for the star's makeup [R0308].

| Piece | Value | Source |
|---|---|---|
| Sequences | 23 masses, 0.20 to 1.30 solar masses, thick (q_H = 1e-4) and thin (q_H = 1e-10) hydrogen layer, carbon-oxygen cores; each lists Teff, log g, radius, cooling age, photon luminosity and the crystallised mass fraction from about 100,000 K down to about 1,500 K | [R0308] (`sequences.tar`, sha256 in the catalogue) |
| State at a cooling age | `cooling_state(mass, cooling_age_gy, envelope)`: the two grid masses around the mass, each read at the cooling age (linear in age, in the logarithm of Teff and of L, in the radius and in the crystallised share), then mixed linearly in mass | the files; interpolation is the reader's |
| Cooling age | the star's age minus the end of its asymptotic branch (`tuning.WD_MIN_COOLING_AGE_GY` at least). The sequences count from about 100,000 K; the tens of thousands of years from the end of the branch to there are not counted | stand-in (tier T5, tweak L0) |
| Hydrogen layer | `tuning.STELLAR_WD_DA_SHARE` (0.80) of white dwarfs are hydrogen atmospheres and read the thick sequences, the rest read the thin ones; the draw is kept in the state (`white_dwarf_atmosphere`) | the share is the existing design value of `stellar_composition` (tier T3, tweak L0); BC19, the share against temperature, stays open |
| Older than the last model | held at the last model (about 1,500 K and a few millionths of a solar luminosity); the answer's `extrapolated` is true | stand-in (tier T5, tweak L0): no model of a black dwarf is read |
| Heavier than 1.30 solar masses | read from the 1.30 sequence; `mass_clamped` is true. Oxygen-neon white dwarfs are not in the sequences | stand-in (tier T6, tweak L0) |
| Temperature, luminosity, radius | the sequence's three, so `L = 4 pi R^2 sigma Teff^4` holds to the interpolation (tested to 3 percent) | Stefan-Boltzmann |
| Crystallised share | the sequence's `Mx/M`, 0 to 1, carried in the star's state as `crystallised_fraction` | the files |

## When the sequences are missing

The dataset states its own fallback (`white_dwarf_models.DATASET.fallback`). A run that disregards missing data (GEN.599) takes the old
Mestel law and the Chandrasekhar-type radius, and `fetch.stand_ins_used()` lists "the Mestel power law of the white-dwarf cooling"; any
other run refuses (GEN.590). The law is the labelled fallback (`tuning.WD_COOLING_*`), not a second model beside the sequences.

## What it changes

Stored values: the temperature, luminosity and radius of every white dwarf the population model makes, so the generator epoch rises.
A young white dwarf is hotter than before (up to about 110,000 K) and a 10 Gyr one cooler (about 3,000 to 4,000 K for 0.6 solar masses, below
2,000 K with a thin hydrogen layer); the radius no longer depends on the mass alone. A white dwarf now draws one more uniform (the atmosphere).

## Not yet

- The hydrogen-to-helium atmosphere change with temperature (BC19: DA fraction against effective temperature from the large catalogues) is
  still a fixed share.
- Helium-core and oxygen-neon white dwarfs have no layers or sequences (the sequences are carbon-oxygen); GEN.743.
- A white dwarf placed by a user-given type (`-star_type D...`) and the companions of `multistar.py` still use the old mass-radius relation;
  the population model is the generator's path for a white dwarf primary.
- Star makeup (`stellar_composition.star_makeup`) draws its own atmosphere from the star's key; the star's drawn atmosphere should replace it
  (GEN.553).
