**Headline:** Generation now keeps a world's life stage, pace, ages and flavour line as data instead of prose, and nothing reads the old prose back.

### Changed
- `life.apply_life_data` draws the timeline as a state (`evolution.get_evolutionary_state`: the pace, the latest milestone, the ages) and stores it in the DB.50 columns; `decide_flavor_text` stores `flavor_key`. The draws are unchanged. `Planet.evolutionary_data` and `Planet.flavor_text` are now read from those columns, and `life_stage_from_paragraphs`, the paragraph and keyword searches in `decide_flavor_text`, `population.parse_timeline` (now `timeline_of`, exact ages instead of rounded text) and the paragraph-reading habitability life stage (`habitability_world.native_life`) are gone.
- `hab_note` is no longer stored or parsed by regex: the tooltip note is worked out from the stored tiers (`habitability_world.note`) and a compact object's planet finds its host from the star (`store.compact_host_of_star`). `system_configs.markdown` is no longer read or written. The evolutionary paragraph tables and the prose columns are no longer touched (DB.52 drops them).
- The stored values change (exact ages, no paragraphs, no note), so this rides the pending generator epoch.

generation-output: changed

TODO items resolved: DB.51 / Added: none
