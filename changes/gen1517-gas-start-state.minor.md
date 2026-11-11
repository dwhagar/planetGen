**Headline:** A giant that gathers its gas before the vapor of a distant ice line could reach it now keeps the disc's own gas, so it no longer comes out with no carbon dioxide or noble gases.

### Changed

- The pebble-drift model (GEN.1473) started the disc's gas ice-free, so a giant that gathered its gas early had none of the vapor of any line far from it (a Sun-like giant: no argon, krypton, xenon or CO2). At the start the gas inside each line holds that line's vapor; it stays so until gas from outside the line has flowed in. `disk_drift.enrichment_ratio` now returns the static weight of 1 while the age is under the time the gas needs from the line to the place, and the lagged pebble ratio after it (GEN.1517, stand-in DE11).
- Giants of Sun-like and heavier stars come out nearly all as the static disc (18 of 24 grid cells for 1 solar mass); the drifted C/O stays for giants of dwarfs and for giants that gather their gas later than the lag. Only the carbon-to-oxygen ratio of the envelope changes in this note (the noble gases follow in GEN.1516).

generation-output: changed
generation-affects: fill
