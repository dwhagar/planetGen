**Headline:** The overlap verifier no longer reads a trojan planet or a trojan region as crossing its giant's orbit.

### Fixed
- **The overlap verifier did not know trojans (GEN.1885).** The overlap verifier (`physics/clearance.py`) now knows trojans (GEN.1885). A trojan planet and a trojan region stand on their giant's orbit by design (GEN.1150, GEN.1151), but the verifier read them as a planet crossing its giant (spacing 0 mutual Hill radii) and as a planet reaching into a belt, which made `test_clearance` pair-6 and triple-10 red on main, and the trojan regions sat in the ring `belt_order` counts, which could point a promoted planet at the wrong belt. A trojan planet is left out of the planet crossing and spacing rules (its giant stands for it), a region is not held against its own giant, a planet that is not the giant's is still held to the region, and a region is no belt of the ring. The limits are unchanged. No generated output changes.

generation-affects: none
generation-output: unchanged
