### Added
- Schema v88: `species` gains the six tech indexes, `tech_level` (indexed for range search) and `tech_source`, ready for the tech-level pass.

### Fixed
- The harshness of a world (used for the tech level draw and facility ratings) read the air score under its old name and always fell back to the middle value; it reads `phi4_air` again.
