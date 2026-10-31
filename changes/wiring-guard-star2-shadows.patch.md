**Headline:** Star 2's eighteen functions that only tests called now each have an owner: three deleted, two named to the photometry library, one to a star page item, twelve marked as test readers.

### Changed
- **Deleted three functions only tests called (TEST.523).** `stellar_evolution.main_sequence_end_luminosity_sol` (its tests read `main_sequence.end_state(mass, fe_h).luminosity_sol`, the one line it returned), `stellar_abundances.solar_metal_fraction` (the value stays as the private Z0 the metal ratio divides by; its test reads `composition(0.0, {}).z`) and `stellar_composition.cluster_key` (the cluster test spells the key `"cluster:33"`).
- **`docs/design/research/wiring-register.md`: eighteen lines no longer name TEST.520.** `photometry.blackbody_colour` and `photometry.bolometric_correction_v` name VIEW.11 (the photometry library reads them); `stellar_composition.star_makeup` names UX.254 (a star's chemistry on its page); twelve are `by design` with the reader named: the sampler, expected-share and published-table checks of their tests, the plausibility gate, and the neutron star table build script.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (TEST.523 closes when this merges; the TODO thread owns the file).
