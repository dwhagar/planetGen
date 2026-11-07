**Headline:** The habitability test galaxy uses a system that still has moons and scores only the bodies that carry domain ratings, so the habitability browser tests pass again.

### Changed

- **Habitability fixture follows the moon change (TEST.677).** The seeded system (seed 6) held five moons before train 14 and one after it, since captured moons moved to the irregular store (GEN.1082 to GEN.1084), so `test_the_lists_draw_each_bodys_stored_rating_as_the_symbol_of_its_tier` for moons found one scored moon where it needs more than three. The fixture now uses seed 9 (seven planets, six moons). That system has three gas giants, which have no domain ratings, and the fixture gave them a score, so the lists cell raised on a missing tier; the fixture now scores only bodies that carry domain ratings, as a real galaxy does. Test fixture only; no product code changes.

generation-output: unchanged
