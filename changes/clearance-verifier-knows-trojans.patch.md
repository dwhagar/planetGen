**Headline:** The overlap verifier no longer counts a trojan region among the belts a promoted planet's `belt_order` points into.

### Fixed
- **A trojan region sat in the ring `clearance._came_from` counts (GEN.1885).** A promoted planet's `belt_order` is the belt's place among the belts of its star, and a trojan region (GEN.1150) is made after the moons and is no belt of that count (`irregular_moons.belts_of`); the verifier counted it, so a trojan region nearer in than the first belt moved every belt's place by one and pointed a promoted planet at the wrong belt (or none). The ring now leaves out a region with a `lagrange_point`, as the count does. The rest of GEN.1885, a trojan planet and its region standing on the giant's orbit by design, is GEN.1882, which merged with Train 64. No generated output changes.

generation-affects: none
generation-output: unchanged
