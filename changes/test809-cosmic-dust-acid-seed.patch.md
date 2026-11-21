**Headline:** The cosmic dust tests find a default system with an acid cloud by searching seeds, instead of relying on one seed that has stopped holding one.

### Fixed

- `tests/test_cosmic_dust.py` built the default system of seed 32 and expected a hot rocky planet with an acid cloud (V1) in it; that seed no longer has one. The fixture now tries seed 137 (which holds one) and then seeds 0 to 199, stops at the first system with an acid cloud, and fails with the tried count when none has one. The cloud assertions are unchanged (TEST.809).

generation-output: unchanged
