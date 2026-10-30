**Headline:** The rules for bodies at the Lagrange points L4 and L5 exist as a pure module, so trojan swarms and whole bodies can be wired later.

### Added
- **Lagrange trojans (GEN.1149).** `physics/trojans.py` holds the Gascheau-Routh mass rule (a trojan needs the secondary plus itself under 0.0385 of the three masses), the libration period and tadpole against horseshoe, the split of the scatter's trojan share into the leading (L4) and trailing (L5) regions, and the map label `Jupiter-L4`. Its numbers are in `physics/data/trojans.json` and stand-in rows TJ1 to TJ5 (registry 6.121). `docs/design/systems/lagrange-trojans.md` describes it and what GEN.1150, GEN.1151 and MAP.354 wire.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (GEN.1149 is this module; GEN.1150, GEN.1151 and MAP.354 wire it).
