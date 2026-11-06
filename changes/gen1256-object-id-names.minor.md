**Headline:** Rogue planets, black holes, neutron stars, nebulae, remnants, quasars, comets and asteroid fields are now named from their Object ID, shown as a two-word local name with the sector beside it.

### Changed
- An object with no star stores its Object ID (`FE81000A2B-0000005-000`) as its `name`, not a 19-digit position ID. The API shows it as its local name under the naming key and adds `name_sector`, the birth sector. A key change still renames all of them at once, and nothing is stored in words.
- `GET /api/lookup?text=` also accepts a local name followed by its sector (`Banhunm Finkap, sector 00FA0FFFE9E240`).
- Migration `0214_object_id_names` rewrites the stored names of an older galaxy from each row's `uid`; a row with no `uid`, and a name given by hand, stay.

- The golden files change (epoch 134): the stored names of the rogue planets in the pinned galaxy change shape, so the generator epoch goes up.

### Removed
- The position ID of GEN.64 (`names/object_id.py`, the type codes, the 16-number collision chain, `store._claim_object_ids` and `_rank_chain`, the `object_ids` table order, `registry.id_code`) and the codec name of a position ID (`naming_key.codec_name`, `object_id_of`, `stored_name_for`). There is no alias.

generation-output: changed
