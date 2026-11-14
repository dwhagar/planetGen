**Headline:** A body keeps its Object ID when its system is drawn again, even when an earlier draw left its place empty: the number of each place in a system is remembered.

### Added

- `body_numbers` table (migration 0257): for a system, the path of a body's slot (`planet:primary:3`, `moon:planet:primary:3:1`, a star's role and place, a belt's or comet's place) and the body number it was given. Redrawing a stored system (`replace_star_system_content`) writes the numbers its bodies hold first, and a new body at a remembered slot takes that number. Before, a body kept its ID only if a body held the same slot at the moment of the draw; a planet that came back at a slot an earlier draw had emptied got a new number. Deleting a system deletes its rows. Left out of the fingerprint. Nothing is backfilled: a system's first redraw writes its rows.

### Changed

- DB.75 first part (`_slot_uids`/`_carry_uids`) is now `_remember_body_numbers`/`_carry_body_numbers`, reading the stored numbers instead of the rows alive at the moment of the draw. A fresh system numbers its bodies in the same order as before, so no ID changes in a galaxy made by this version. The key becomes Lane 05's `{kind}:{host_seed}:orbit:{slot}` with GEN.1548 (DB.75b); the slot is today's place in the star's orbits.

generation-output: unchanged
generation-affects: none
