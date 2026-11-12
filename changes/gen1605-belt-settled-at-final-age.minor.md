**Headline:** A belt's mass, largest body and density now match the star's final age and the belt's final place, and its largest body is its biggest member.

### Changed
- GEN.1605 (bug): a belt settled its mass, largest body and density once, when it was built, at the age the star had then and the place the belt had then. The star's age rises afterwards (`Star.adjust_age_for_planets`) and the spacing pass can move a belt, but the members were drawn later at the final age, so the stored mass and the members described two belts (5 of 25 belts stored a mass 2 to 4.5 times off), and late accretion and the debris disk read the stale mass. `AsteroidBelt.settle` now works the three out again at the final age and place in the settle pass (`StarSystem._settle_belts`, beside the secular rates), and a belt it does not change is left as it is.
- A belt's `largest_body_km` is its rank 1 member's diameter (the ledger's figure is the diameter at the grain density its collision law uses, the member's is the one its taxon's density gives it; they disagreed by 15 to 40 percent). A belt without members keeps the figure the ledger or the store gave it.
- `tests/test_belt_settle.py`: a belt settled again at a later age holds what a fresh belt of that age holds; the finished belts of seeded systems carry the mass of a fresh settle and a largest body equal to their rank 1 member's (red before, green after).

generation-output: changed
generation-affects: fill
TODO items resolved: GEN.1605
