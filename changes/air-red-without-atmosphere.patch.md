### Fixed

- **A world with no atmosphere now rates Red in Air, with a score of 0 (GEN.215).** It rated Green because the air had no toxins. Under 0.05 kPa (it prints as 0.0 kPa) there is no air to breathe, and the description says there is no atmosphere. Bodies generated before this release keep their stored Air rating until regenerated.
