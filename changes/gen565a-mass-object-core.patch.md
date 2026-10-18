**Headline:** Every body with a mass now shares one core class for its galactic motion and its save and load, the first step of the one object class for stars, planets, comets and wanderers.

### Changed
- New `generation/mass_object.py` (`MassObject`) holds what the nine body classes wrote out by hand: the four galactic-motion fields and the draw that sets them, and one `to_dict`/`from_dict` pair driven by each family's `SERIALIZABLE_FIELDS`, with `_encode` and `_decode` hooks for values that are not plain JSON. `Star`, `Planet`, `Comet`, `RoguePlanet`, `InterstellarComet` and `AsteroidField` use it; `BlackHole`, `NeutronStar` and `BinaryStarProxy` inherit it through `Star`, and `CompactRemnant`'s copy of `Star`'s save and load is gone.
- No stored value, field name or draw changes.

generation-output: unchanged

TODO items resolved: none (GEN.565 stays open) / Added: none
