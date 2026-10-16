**Headline:** The cell key is now tested across a cube boundary on every table that holds a position.

### Changed

- New `test_cellkey_motion.py`: the orbit update carries a star system, a nebula, a facility, a bright star and a black-hole scatter point 0.04 pc over a cube face, and a hypervelocity star over a face on each axis, and the stored key is the next cube's, one cube away on that axis only; a step inside a cube leaves the key where it was.
- The same file writes awkward points (on a face, one step either side of it, off the box, `0.1 + 0.2`, the smallest double) to one row of every table with a cell key and reads back the same bits and the same key as `cellkey.cell_key`.

TODO items resolved: TEST.166 / Added: none
