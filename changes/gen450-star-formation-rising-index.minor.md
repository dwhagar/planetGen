**Headline:** The disc's recent star formation now falls faster in its thin outer gas, as the observations of molecular and atomic gas say, instead of at one fixed index.

### Changed

- The share of a thin-disc star that is recent (`galaxy.star_formation.recent_factor`) now has the star-formation law's index rise where the disc's own gas thins (GEN.450): the rate goes as `gas^1.4 (gas / (gas + 9))^(2.1 - 1.4)`, the gas being the rings' gas of the disc's star-formation history. Molecular gas forms stars at a constant efficiency, atomic gas (under about 9 solar masses per pc2) at a falling one (Bigiel et al. 2008, Leroy et al. 2008, Schaye 2004). For a Milky Way-like disc the fall beyond the rate peak is -0.27 dex per kpc (was -0.23; the Herschel study observes -0.28, which is a test target and not an input); 91% of the star formation is inside the solar circle (was 86%, observed 84%) and 0.4% beyond 13.5 kpc (was 1.3%, observed 1%). The young and intermediate populations are lighter in the outer disc of every galaxy, so a placed star there is older on average.
- `tuning.STAR_FORMATION_LOW_DENSITY_EXPONENT` (2.1, the middle of the 1.5 to 2.7 that resolved total-gas fits give) and `tuning.STAR_FORMATION_SATURATION_GAS_MSUN_PC2` (9) are new and open: no source gives them for one galaxy.
- `tests/fixtures/density_identity.json` is recorded again for the new factor.

generation-output: changed
