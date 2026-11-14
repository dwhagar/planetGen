**Headline:** A planet or moon now takes its air, class and spin once, on the orbit and size it keeps: the tides, the Kozai pumping, the forced eccentricity, the giant's radius and the envelope all settle first, so about 100 air remakes in 25 systems are gone.

### Changed
- **The settle pass settles the orbits before any air (GEN.1548, slice 2).** `StarSystem._settle_to_final_age` now: gives every planet its orbit state, runs the secular exchange, gives the belts their planets and the chains their forced eccentricities, settles each planet's tide, giant radius, envelope and oddities (a body built bare, `Planet(defer_climate=True)`, has no air to make again), respaces if the eccentricities closed a gap, and only then makes each planet's air, class and spin (`finish_climate`). Its moons are built bare (`generate_moons(bare=True)`, orbit state only), take the Kozai pumping and the tide, and then make their own air, class and spin once (`Planet.finish_moons`), followed by their day and night shares and the sea tides. A system that bars habitable worlds builds whole bodies and keeps the old remake path. `planets.is_bare` names a body whose air is still owed; `settle_giant`, `settle_envelope`, `remake_after_orbit_change` and the moon history leave it alone.

generation-output: changed
generation-affects: fill
