**Headline:** The nitrogen, sulfur and noble gases of a giant's envelope now follow the ices that the pebbles carried inward, as its carbon and oxygen already do.

### Changed

- `disk_condensation.envelope_enrichments` takes the ice weights of the pebble drift (GEN.1516): the vapor of NH3, N2, H2S, argon, krypton and xenon that evaporated on the way to the giant is its weight times the static vapor, so a giant whose gas the pebbles enriched in a noble gas is no longer fixed at 1 times its star's. Weights of 1 give the factors of GEN.1372 exactly (tested at every temperature).
- A giant that gathered its gas before the vapor could reach it keeps the static factors (GEN.1517), so most giants of Sun-like stars do not change; the giants of dwarfs and those that gather late do.

generation-output: changed
generation-affects: fill
