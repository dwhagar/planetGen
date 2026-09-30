### Changed

- **Schema v35: hybrid master-wedge sector slots.** Each ring now holds
  the multiple of its master wedge count nearest `2*pi*(i + 1/2)`
  (3, 9, 15, 21, 27, 36, ...), with 3 master wedges at the center that
  double (6, 12, ... 1,536) once each would hold 8 slots. Slot
  boundaries line up on the master lines from the center to the edge,
  sector arcs stay within about 6% of 4 pc, and the total sector count
  is unchanged. `galaxyGeometry.ring_master_count` is new, and
  `galaxyprisms.js` mirrors both functions.
- The migration deletes every sector in a ring whose slot count changed
  (all but 15 rings) with its systems and phenomena. Run `update.sh`,
  then regenerate the galaxy; the skeleton (`generate.py plan`) is kept.
