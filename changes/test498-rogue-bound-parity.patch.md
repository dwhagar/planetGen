**Headline:** A test now checks that a rogue planet and a planet of one mass are made by one code path and stored the same way.

### Changed
- TEST.498: a rogue planet and a bound planet of the same mass are compared on the attributes they carry, the columns and kinds of value they store, the moon class and the class table; they differ only in the star's light, the liquid surface, the cloud layer and (through the class) the strength of a solid body. Nothing derives from `Planet`, so no rogue-only generator can exist unnoticed.
