**Headline:** A world whose only liquid is a puddle of concentrated chloride brine now reads "Brine seep" instead of "Surface ocean", and it no longer counts as a sunlit sea for life.

### Changed
- **Brine seep (GEN.908).** A liquid surface with so little water (under 2e-5 of the body's mass) that its chemistry is a concentrated chloride brine, melting at 200 to 270 K and often under no air, is its own water state, `brine seep`, on the page, the map and in the search. It still holds liquid water for a base and its ocean details are shown, but it is not a sunlit sea for life, not a wet world for plate tectonics and not an induction ocean. About two thirds of the worlds that used to read "Surface ocean" in a sun's habitable zone were such seeps.
- **Database.** Migration 0176 lets `planets` and `moons` store the new state and relabels the stored brine oceans.

generation-output: changed
