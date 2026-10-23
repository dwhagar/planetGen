**Headline:** The collision number of an object ID now follows the sectors' addresses, so two sectors that share an ID get the same names whichever one saved first.

### Changed
- **The object-ID collision number is keyed on the lower address (GEN.307).** Objects named by their 76-bit position ID (rogue planets, black holes, neutron stars, nebulae, supernova remnants, quasars, comets, asteroid fields) that share an ID take the collision numbers 0, 1, 2 ... in the order of their sectors' addresses, and in one sector in the order they were generated. `store._claim_object_ids` took the next free number, so when two sectors held objects with one ID the sector that saved first kept number 0. It now gets the sector's address (`save_sector` passes it): when a stored object of a higher address holds a number the new object outranks, the stored object moves up to the next number (`_rank_chain`). An object saved on its own, with no sector address, ranks after every sector, as it did. Every quasar shares the galactic centre's ID, so the quasars of a galaxy are the case this settles.
- `tests/test_object_id.py` saves three sectors with one shared ID in three orders into three databases and finds the same numbers 0, 1, 2 by address in each.

generation-output: changed

TODO items resolved: GEN.307. Added: none.
