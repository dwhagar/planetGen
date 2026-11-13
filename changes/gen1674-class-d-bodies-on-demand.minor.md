**Headline:** Every body of a belt at the class D size or larger now has a page and a class when you ask for it, however many the belt holds, made on the spot by the same generator that makes the planets and never stored.

### Changed
- **Bodies past the stored 300 (GEN.1674).** A belt keeps up to 300 bodies, but its size law holds many more over 1e19 kg (about 200 km across): about 230,000 in a median outer disk (Kuiper and scattered disk) and up to over a million. A rank past the stored ones is now a body of the same law, drawn from a stream of the belt's own, the same whenever it is asked for, down to the 1e19 kg floor. `GET /api/objects/<ref>/belt-members/<rank>` and the member page serve it with the same page as a stored body.
- **Its class, on demand.** The member page and the API now give the body's class (`planet`: class, its description, mass, radius, density, gravity, surface temperature and air pressure) from the one planet generator (`belt_planets.planet_of`), for a stored body and a tail body alike; nothing is stored for it. The belt block gives `class_d_count`, the number of bodies the belt holds at the class D size or larger.
- No stored value, no generated body and no generation output changes: the member rows and the planets are as they were. `belt_roster` and `belt_page` give the one function that makes a body's entry to both the stored bodies and the bodies drawn on demand.

generation-output: unchanged
generation-affects: none
