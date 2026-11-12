**Headline:** The gaps a giant cuts in a belt now come from its mass, orbit and the system's age, and members are never drawn inside them

### Changed
- GEN.1092 (part B of the resonance wiring): the span a giant clears from an asteroid belt is now worked out by the resonance model (`physics/resonance_gaps.py`: the chaotic zone round the giant and the mean-motion resonances 4:1, 3:1, 5:2, 7:3, 2:1 and the rest up to fourth order, each with a width that grows with the giant's mass and eccentricity and a clearing time of 10 Myr). Before, `belt_origin` used its own stand-in numbers (a chaotic zone of 1.3 mu^(2/7), five fixed resonances, a fixed 1.5 percent half-width).
- The cleared share of a belt's mass and the gaps follow the same model, so a belt round a small giant has narrower gaps than before, and the exterior resonances of an inner giant now clear mass from the outer disk. A system younger than the clearing time keeps the gaps its giant has not yet cut.
- A belt's stored members take their semi-major axes outside the cleared spans (`AsteroidBelt.cleared_spans`, `belt_roster.semi_major_axis_au(table, u, cleared)`), so no member sits in a gap of the belt it belongs to. The Yarkovsky size cut is not applied to stored members, which are all far above it.
- Belt borders, belt masses, the number of belts of some systems and the stored members move for systems with a giant.

generation-output: changed
generation-affects: fill
TODO items resolved: GEN.1092
