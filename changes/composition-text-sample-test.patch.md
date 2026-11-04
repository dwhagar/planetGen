**Headline:** The composition text test draws on until it has seen both a rocky and a giant body.

### Changed

- `test_every_generated_body_states_what_its_fractions_say` drew six systems and asserted it had seen a rocky and a giant body; under the current planet law those six hold no giant. It now draws at least six systems and keeps going (at most 60) until both kinds have been seen. The per-body check is unchanged. Tests only.

generation-output: unchanged
