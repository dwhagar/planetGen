**Headline:** There is now a short two-word local name for any Object ID, spelling only its serial and body number, with the sector shown beside it.

### Changed

- **Local names (GEN.1255).** `names/naming_key.py` gains `local_name`, `local_name_with_sector` and `uid_of_local_name`: the serial and body number of an Object ID, 10 hex digits, spelled alone as two words in their own codec domain, so the sector's words never appear (`Zaidonr Jaisonn, sector 0000B0001803E8`). The 26 count bits of the serial are first scrambled by a keyed bijection (odd multiplier and offset drawn from the naming key and the sector), so serials 0, 1, 2 no longer share a first word: the first 4,096 serials of every sector have first words of their own. The name is exact inside a sector and reverses from the name, the sector and the key. `CODEC_VERSION` goes to 2. Nothing uses the local name yet: rogue planets and the other phenomena, scatter rows, the pages, the API and search follow in GEN.1256, GEN.1257, UX.259 and API.38. No stored data or generated output changes.

generation-output: unchanged
