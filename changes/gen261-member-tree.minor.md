**Headline:** A star, black hole or neutron star can now be given a stable tree of companions: stars, brown dwarfs, white dwarfs, neutron stars and black holes in hierarchies of up to twelve members.

### Changed
- A new generator, `draw_member_tree`, draws the members and orbits around a placed object from its own random stream: the number of stars by the object's spectral class (with a tail to twelve stars and a prevalence knob for five or more), each member's state at the system's age (protostar, pre-main-sequence, main sequence, giant, white dwarf, brown dwarf), a compact companion at the table odds, and every orbit top down under the Mardling-Aarseth, Eggleton-Kiseleva and Vynatheya stability criteria and the Kozai-Lidov screen. Nothing the sector fill or the scatter makes changes yet; the next steps call it.
- A new stability module holds the three triple-stability criteria behind one switch, the Kozai-Lidov screen with its relativistic quench, and the Peters merger time. Compact pairs that would already have merged are rejected, and a bound companion beyond 1,000 AU of a neutron star or a black hole that did not collapse directly is rejected.
- The system-type table test now requires a generator to exist once the item that builds it has landed, instead of requiring every generator the moment the module exists.

TODO items resolved: GEN.261 / Added: none
