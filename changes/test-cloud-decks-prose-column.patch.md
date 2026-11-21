**Headline:** The column check lists `planets.cloud_decks` as a JSON document read whole, as it does `air_history`.

### Fixed
- GEN.1577: `planets.cloud_decks` (a giant's cloud decks and trace-gas ceilings, JSON text) is in the prose-column list of `test_db_columns_and_indexes.py` with its reason, so the two tests that failed on the Train 51 head (`test_no_galaxy_table_keeps_a_json_or_serialized_blob`, `test_the_prose_list_has_no_stale_entries`) pass.

generation-output: unchanged
generation-affects: none
