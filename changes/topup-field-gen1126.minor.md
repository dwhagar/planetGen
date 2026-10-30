**Headline:** The mass the Rogues level cannot hold in a sector can now become thin clouds of small rubble, planned as a tested calculation that places exactly the mass that is left; no sector uses it yet.

### Added
- **Light top-up fields (GEN.1126).** A new module (`physics/topup_field`) takes the mass the rogue top-up leaves over and splits it into a few fields (a dozen at most) of smaller bodies, lower density and a larger volume than a field with a history, each placed in the sector with the radius and age of a field that has not yet retired, and with the same apparent origin and velocity a top-up rogue gets. Their origin is the registered stand-in "unbound remnant", with no collision event and no parent claimed. The placed mass is the remainder exactly. Four stand-in rows (TF1 to TF4) are registered. The rogue top-up shares its apparent-origin draw (`rogue_topup.apparent_history`, same numbers as before). Nothing reads it yet: the region run (GEN.1120) and the field items (GEN.1103 to GEN.1107) wire it on a later generation train.

generation-output: unchanged
