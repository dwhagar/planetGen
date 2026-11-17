**Headline:** A bad FORCE_CLASS request is refused whatever star comes up, since the check now runs before the first draw.

### Fixed
- `StarSystem`: `FORCE_CLASS` combined with `SLOTS`, or naming a letter that is not a planet class, raised only after the star was drawn and only for stars that reach the planet planner, so a pulsar, a post-AGB pair or a diskless black hole accepted it silently (found when `test_class_sample::test_the_switch_refuses_slots_and_unknown_letters` went red on main after a draw moved). The checks are the same, now made once at the start. The test's "unknown letter" was `Q`, which is a planet class now (every letter A to Z is one), so it names `1`.

generation-output: unchanged
