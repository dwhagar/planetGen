**Headline:** The planet-input fuzz test checks each class for a kind and a phrase, since the composition strings were retired.

### Fixed

- `test_every_class_has_a_composition_and_type` read `PLANET_CLASSES[letter]["composition"]`, which GEN.1619 removed, so it failed with `KeyError: 'composition'` on main (TEST.798). It is now `test_every_class_has_a_kind_and_a_phrase`: each class keeps its kind ('t' or 'g') in `tuning.PLANET_CLASSES`, and has a phrase in the class table (`planet_classes.describe`). The test was stale, not the code.

generation-output: unchanged
