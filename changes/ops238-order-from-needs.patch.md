**Headline:** the order of a run's stages is now checked and logged from what each command needs, not from a list of four requirements.

### Changed

- OPS.238, second part: `commands.REQUIRES` is deleted. `needs.requires_of` (the kinds that make what a command needs), `needs.required_stage_keys` and `needs.validate_order` take its place: a run whose steps are listed so that a need is met only after the step that has it is refused where the list is made; a stage's "released" line, a failed run's "did not run" line and a command's break-point dependencies (`after`) read the same edges. The four requirements the old table held are still there (a test says so); the run's stages still run in their listed order.
