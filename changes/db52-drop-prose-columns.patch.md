**Headline:** The stored evolutionary paragraphs, flavour text and habitability note are gone from the database, replaced by the structured columns added earlier.

### Changed
- Migration 0111 (schema v111) drops `planet_evolutionary_paragraphs` and `moon_evolutionary_paragraphs`, `flavor_text` and `hab_note` of `planets` and `moons`, and `system_configs.markdown`. Nothing reads or writes them since DB.51; the life stage, pace and ages, the flavour key and the stored tiers carry what they said. `flavor_text_count` stays.
- The old-schema fixture for v110, the schema reference and the column coverage tests follow.

TODO items resolved: DB.52, DB.31 / Added: none
