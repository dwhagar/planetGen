**Headline:** A glowing dark member of an orbit tree (a brown dwarf, a neutron star, a black hole with a disk) is now a source in a planet's light ledger, stored by its kind and member number.

### Added

- GEN.1528: the effect ledger stored only sources with a row of a kind table, so the light of a dark member (`dark_light`, GEN.1501) was in the planet's light and missing from its ledger share. `effect_ledger.dark_location` names such a source and `DARK_KIND_CODES` (brown dwarf 8, neutron star 9, black hole 10) stores it with its member number as `source_id`; `EffectLedger.rows` files it, and the refill job writes it again. `source_type` is a plain small integer, so no migration is needed (the item expected one). A dark member's light is a part in 1e11 of its host's when it is a brown dwarf, so in the stored galaxy it is mostly in the summary row under the floor; rows appear where it is near or hot (a young neutron star, a black hole disk).

generation-output: changed
generation-affects: scatter, fill
