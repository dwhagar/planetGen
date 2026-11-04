**Headline:** The forced belt and habitable world test checks the request's answer, not a requirement line that GEN.1130 stopped writing.

### Changed

- `test_asteroid_belt_and_habitable_world_forced_together_both_succeed` counted a system with no habitable world as a failure unless "a habitable world" stood in `unmet_requirements`. Since GEN.1130 a habitable request only raises the odds and answers through `StarSystem.habitable_request` ("aim", "fallback" or "unmet"), so G2V failed in 11 of 15 trials on main. The test now asks for that answer. The belt half of the test is unchanged. Tests only.

generation-output: unchanged
