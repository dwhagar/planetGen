**Headline:** the two sweep tests that break the life model by class letter no longer fail when the habitability code scores a body from its values.

### Changed

- The stand-ins that `test_sweep_life.py` puts in place of `life_classes.score_classes` read `body.planet_class`; the habitability code calls it with a column-named mapping (`_class_rater`), which has no letter, so they raised `AttributeError` as soon as a sweep body reached that call (train 24). They read the letter with `getattr(body, "planet_class", None)` now and pass a mapping to the real function. Test only.
