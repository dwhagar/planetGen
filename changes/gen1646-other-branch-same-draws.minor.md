**Headline:** The second climate branch of a bistable world now keeps the world's own volatile inventory instead of drawing a new one.

### Changed
- **GEN.1646, part 2 (the cause of the excess nitrogen air).** `air_branches.other_branch` (oddity W20) solved the air a second time on a stream of its own, which drew the inventory, the degassing efficiencies, the sink cycle and the late accretion again, so the two climate branches of one world could hold 20 kPa and 500 kPa of N2. `surface_air.solve` now records the ambient stream state it began from (`Air.entry_state`) and `other_branch` repeats it, so only the start temperature differs. On the Earth analog (varied key, 90 trials) the N2 median and spread move toward the Earth's; the numbers are in the PR. A world with W20 and no distinct second branch under the same inventory is no longer named for one; the pinned test seed moved (`w20-G2V-20`).
generation-output: changed
