**Headline:** Gas giants, ice giants and brown dwarfs can now be given full layered interiors, from molecular and metallic hydrogen to helium rain, dilute cores and rocky cores.

### Added
- **`planetgen.physics.giant_interior` (GEN.642).** Places the layers of a gas giant (molecular and metallic hydrogen
  envelope, a helium-enriched zone and shell where helium has rained out of a cold giant, a dilute core, a compact core
  of rock and iron), an ice giant (molecular envelope over an ionic or superionic water, ammonia and methane mantle over a
  rock and iron core) and a brown dwarf (no solid core), each with its pressures, mass, state and the mass of every material
  in it, and describes them in words. The layers hold exactly the heavy elements and natal helium the caller fixed. Nothing
  calls it yet: the planet draw wires it in with GEN.428 and the layers are stored by DB.89.
