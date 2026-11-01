**Headline:** A companion star stored with a system is rebuilt as the same star: the stored companions of massive stars no longer fall away when the system is filled in.

### Changed
- **A member's record is made from floats (GEN.951).** The companion mass law returned a numpy number, and the star model gives a slightly different star (about three parts in a hundred million) for a numpy number than for the float of the same value. The scatter stored the first star and the fill rebuilt the second, so the check that the built star is the stored member failed and the system was left a single star. The record and the companion mass are floats now, so both sides make one star.

generation-output: changed

TODO items resolved: GEN.951.
