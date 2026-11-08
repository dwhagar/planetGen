**Headline:** The one-line description of a planet class and its class-page sentence now live in one row, with no number typed by hand and no Solar System feature stated as a fact of the class.

### Changed

- GEN.1432: each row of `planet_classes.json` now holds the class's `phrase` (what a planet or moon page says after "It is") beside its page `sentence`; `planet_class_text.json` is gone and `planet_class_text.py` reads the phrase from the class rows. The two had disagreed for about ten of the 25 letters (class A said "small, barren" though it asks for sulfur dioxide in the air; class F said "volcanic with shallow seas" though it asks for an air with almost no oxygen). The phrases are reconciled to what each class's conditions say.
- The class-page sentences hold no digits: the temperatures, pressures and masses they typed by hand are gone, and the page's Conditions line (read from `planet_class_fit.json`) carries the numbers. The loader refuses a row whose sentence has a digit or whose phrase is empty or starts in capitals.
- Class J no longer says every gas giant has rings, class I no longer says every ice giant has a tilted magnetic field, and classes N and U no longer say "Venus-like" and "Titan-like". A test checks that no class text names a Solar System feature.
- Nothing stored changes: the phrase is read when a page is built.

generation-output: unchanged
