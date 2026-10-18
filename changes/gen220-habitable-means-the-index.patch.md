**Headline:** A system's overview now counts potentially habitable worlds by their habitability rating, and the "Habitable" chips say they mean the planet's class.

### Changed
- The overview counts a rocky planet or moon as potentially habitable when its overall PHI-6 colour is Blue or Green (`habitability.POTENTIALLY_HABITABLE_TIERS`, one place to move the line), instead of counting the class flag. The sentence no longer talks about class M worlds, and the planet count is always listed.
- The row chips read "Habitable class" and "Habitable-class moon", the class pages label the flag "Habitable class", and the levels page explains that the equipment level beside the chip is the rating.
- `GET /api/systems/<id>/description` is cached under the system's tags in place of the removed `sections`.

TODO items resolved: GEN.220 / Added: none
