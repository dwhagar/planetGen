**Headline:** Gas giants, ice giants and brown dwarfs can now be given full layered interiors, from molecular and metallic hydrogen to helium rain, dilute cores and rocky cores, and any body can be asked whether it can have satellites.

### Added
- **`planetgen.physics.hosting_gate` (GEN.579).** `can_host(host)` runs the three parts of the hosting gate in order (a 10 m
  satellite must fit under the mass-ratio limit; a stable orbit must exist between the host's radius or fluid Roche limit and
  0.4895 of its Hill radius, with the hierarchical stability screen for a comparable satellite; the body must reach the 1 km
  size cut and its population and size must have a satellite share), names the first part that fails, and carries the gated
  `satellite-share` table with every row's source state. A moon holds nothing until a sourced row exists, and
  `tidal_survival_radius_m` is the submoon tidal gate of Kollmeier and Raymond 2019. Nothing calls it yet.
- **`planetgen.physics.giant_interior` (GEN.642).** Places the layers of a gas giant (molecular and metallic hydrogen
  envelope, a helium-enriched zone and shell where helium has rained out of a cold giant, a dilute core, a compact core
  of rock and iron), an ice giant (molecular envelope over an ionic or superionic water, ammonia and methane mantle over a
  rock and iron core) and a brown dwarf (no solid core), each with its pressures, mass, state and the mass of every material
  in it, and describes them in words. The layers hold exactly the heavy elements and natal helium the caller fixed. Nothing
  calls it yet: the planet draw wires it in with GEN.428 and the layers are stored by DB.89.
- **Three research notes.** `docs/design/composition/giant-interiors.md` (GEN.641: layers, metallization, helium rain, dilute and compact
  cores, ice giants and brown dwarfs, with each source's state); `docs/design/composition/giant-heavy-element-law.md` (the heavy-element
  mass law of Chachan et al. 2025, now read: GEN.460's input); `docs/design/composition/small-body-size-laws.md` (GEN.566 and GEN.580:
  size-frequency laws of belts, comets, Kuiper belt objects and interstellar objects, and binary and satellite fractions by
  size, with the numbers read in each paper). The size-law note corrects the interstellar object density in
  `unified-star-like-objects.md`: 0.1 to 0.2 per au3 is 7e14 to 1.8e15 per pc3, not 1e12.
- **The unit docs describe the registry (DOC.70).** `docs/design/units-and-number-formatting.md` now says where the units live
  (the unit registry, `Quantity`, the relations, the field declarations, the column registry, the display ladders, the
  browser copy) and that the registry is the only unit table; `docs/database-schema.md` has a "Units of a column" section.
