**Headline:** Planets, moons and rogue planets now share one `PlanetaryObject` family class, the third family under the mass object core.

### Changed
- **`generation/planetary_object.py`** (GEN.565f): `PlanetaryObject(MassObject)` holds what a planet, a moon and a rogue planet share: the class letter and composition defaults, `is_giant`, and `parent_kind` (`"star"` for a planet, `"planet"` for a moon, `None` for a rogue planet). `Planet` and `RoguePlanet` subclass it; `RoguePlanet.kind_label` reads `is_giant`. No field, stored value or draw changes. GEN.565 stays open for folding the remnant classes into states and the table revision.
