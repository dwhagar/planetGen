**Headline:** What an asteroid field from a rogue-planet collision is made of now follows the two parents (the outer layers leave first: ice, then silicate mantle, with the iron staying in the remnant) as a tested calculation; no world uses it yet.

### Added
- **Collision field composition (GEN.1102).** A new module (`physics/field_composition`) takes two differentiated parents (mass, core share, ice share) and a collision speed, runs them through the collision law, and gives the mass and layers of the field that flies free, the remnant planet that stays as its own rogue planet with its core share, and the bound debris. A field is silicate with no iron, and ice shards when the parents were water-rich; a merger, hit-and-run or erosion makes no field. `small_body_composition` gains `makeup_from_components`. Three stand-in rows (IF23 to IF25) are registered. Nothing reads it yet: the field generation (GEN.1103) replaces the one-taxon draw `draw_field("collisional")` with it on a later generation train.

generation-output: unchanged
