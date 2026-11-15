**Headline:** A belt body made a planet no longer lies past the orbit the place's passing stars let a planet keep.

### Fixed
- `belt_planets.promote` checked a promoted member's reach against the star's stable limit only, so in a dense place (a core or a cluster, where the survival limit of GEN.325 is a few AU) a promoted belt planet could lie past the limit the ordinary planets keep. It now uses the same ceiling as the planet slots (`StarSystem._orbit_ceiling_au`). In the Sun's neighbourhood the limit is hundreds of AU and nothing changes.

generation-output: changed
generation-affects: fill
