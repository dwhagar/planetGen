**Headline:** A property test now keeps the habitability words and the habitability index saying the same thing, and the cooling-garment sentence reads "are needed" and the airless pressure reads like every other figure.

### Changed
- New `tests/test_wording_agrees_with_the_index.py`: over random worlds, the example worlds and generated systems, a row's colour is its score's band and any "inside the X band" in its sentence names that colour, the gear the sentences ask for is the stored kit, an airless body is Red in Air, a body whose text claims life has life in the index, a ratings block equals the stored scores, and no quantity in the describer or `protection` is written with a raw `:.Nf`.
- `protection` wrote the airless pressure as `0.0 kPa` and "cooling garments is needed"; it now uses the unit ladder (`0 kPa`) and "are needed".
- One finding stays a strict expected failure in that test: the words ask for an environment suit in a Yellow heat (45 to 122 C), a Yellow cold or a Yellow oxygen excess, where the stored kit (`equipment_tier`) is a mask or nothing. Changing the stored kit moves stored values (a generator epoch), so Boss decides which side moves.

TODO items resolved: TEST.155 / Added: the words and the stored kit disagree in a Yellow heat, cold or oxygen excess (strict expected failure in tests/test_wording_agrees_with_the_index.py)
