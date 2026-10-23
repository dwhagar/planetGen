**Headline:** Gas giants, ice giants and brown dwarfs can now be given full layered interiors, from molecular and metallic hydrogen to helium rain, dilute cores and rocky cores.

### Added
- **`planetgen.physics.giant_interior` (GEN.642).** Places the layers of a gas giant (molecular and metallic hydrogen
  envelope, a helium-enriched zone and shell where helium has rained out of a cold giant, a dilute core, a compact core
  of rock and iron), an ice giant (molecular envelope over an ionic or superionic water, ammonia and methane mantle over a
  rock and iron core) and a brown dwarf (no solid core), each with its pressures, mass, state and the mass of every material
  in it, and describes them in words. The layers hold exactly the heavy elements and natal helium the caller fixed. Nothing
  calls it yet: the planet draw wires it in with GEN.428 and the layers are stored by DB.89.
- **Two research notes.** `docs/design/giant-interiors.md` (GEN.641: layers, metallization, helium rain, dilute and compact
  cores, ice giants and brown dwarfs, with each source's state) and `docs/design/small-body-size-laws.md` (GEN.566 and
  GEN.580: size-frequency laws of belts, comets, Kuiper belt objects and interstellar objects, and binary and satellite
  fractions by size, with the numbers read in each paper). The second corrects the interstellar object density in
  `unified-star-like-objects.md` (0.1 to 0.2 per au3 is 7e14 to 1.8e15 per pc3, not 1e12).
