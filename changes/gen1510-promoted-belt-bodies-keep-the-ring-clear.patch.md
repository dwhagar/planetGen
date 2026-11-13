**Headline:** A belt body made a planet no longer crosses the planets and belts beside it, so no generated system holds an overlap the clearance check finds.

### Changed
- **A belt member becomes a planet only when its orbit stays in its belt's ring and clear of the other planets (GEN.1510 follow-up).** The members chosen for promotion had orbits drawn and then forced up to an eccentricity of 0.5 whatever the ring, so a promoted member crossed the planet beside it or reached into the next belt: `tests/test_clearance.py::test_a_generated_system_holds_no_overlap` failed on 34 of 95 triple and five-star seeds. A member now has to have its whole orbit between the belt's limits (`belt_planets.in_ring`) and ten mutual Hill radii from every planet that is not of its belt (`belt_planets.clear_of_planets`), else it stays a member and the next heaviest is tried. On 30 sample systems the number promoted is the same (231 of the 24 belts in them); only which members are chosen changes, where a chosen member crossed. The independent verifier (`physics/clearance.py`) reads `planets.belt_order` and does not ask the planets of one belt to clear each other or their own belt, which cross by being a belt; every other pair is asked as before. All 72 sweep systems are clear.

generation-output: changed
generation-affects: fill
