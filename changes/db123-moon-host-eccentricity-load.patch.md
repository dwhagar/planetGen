**Headline:** A moon loaded from the database keeps its host planet's eccentricity, so its air solved again gets the light it had.

### Changed
- **`host_eccentricity` of a moon survives a save and load (DB.123).** The field is not a column and came back as 0.0, so an edit that solved a loaded moon's air again (`physics.planets.calculate_atmospheric_conditions`, which takes the light of a moon from its parent's orbit, GEN.727) used a circular parent orbit and got another temperature than generation did. Generation tells each moon its parent's eccentricity (`moon.host_eccentricity = planet.eccentricity`); `store.load_star_system` now does the same for each loaded moon from the loaded planet, so no column and no migration are needed. A moon of a rogue planet has no star and keeps 0.0. Stored data is unchanged.

generation-output: unchanged
generation-affects: none
