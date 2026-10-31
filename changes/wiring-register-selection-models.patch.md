**Headline:** The wiring register says why the simulated survey samplers are only run by the tuning tests.

### Changed
- **Twelve test-only functions are named by design (TEST.527).** The observed samplers and expected counts of the five selection models (comet admission, JWST target selection, Kepler recovery, mass and radius sample, polluted white dwarf detection) are the simulated side of the tuning harness's checks, reached through `tuning_harness.SELECTION_MODEL_MODULES`. Each register line now names who runs it and why the generator never does. No code or generated output changes.

generation-output: unchanged

TODO items resolved: TEST.527.
