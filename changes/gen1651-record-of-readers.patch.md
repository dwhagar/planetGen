**Headline:** the system scene and the pages that show a planet's air colour work again after the build record change.

### Fixed

- GEN.1651 follow-up: a planet restored for display is a plain namespace of its columns with no `build`, and the readers that take any object (`body_radiation._star_distance_au`, reached through the scene's air colour, and sixteen more in `physics/`) read `body.build.<field>` and raised `AttributeError`. They read `record_of(body)` now: the body's record, or an empty one for an object that built nothing, so they give what the old `getattr(body, "_x", None)` gave. Passes that write a field still use the real record of a real body.

generation-output: unchanged
