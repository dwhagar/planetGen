**Headline:** The stale composition strings of the planet classes are gone; a class is its letter and its kind.

### Changed
- **GEN.1619: `tuning.PLANET_CLASSES` no longer carries a `composition` string for each class (they were out of date: class C anthracite, class V iridium and tungsten, class B potassium).** Nothing displayed them and no draw read them; the class pages read the words in `physics/data/planet_classes.json`. The mapping from those strings to material rows (`materials.class_composition`) goes in GEN.1660 (PR 6333), which this change follows. Nothing generated moves.

generation-affects: none
