**Headline:** What can accompany or orbit an object is now a set of host classes with gated tables of possibilities, replacing the registry's five fixed system models; nothing the generator draws changes.

### Changed
- `phenomena/registry.py`: the `system_model` column is now `host`, the name of a host class (group node, star-like object, remnant, brown dwarf, planet, moon, small body) in the new `phenomena/hosts.py`. Each class points to a companion table and a satellite table in the format of `possibility_table` (GEN.472), with each row's source and source state and its zone rule.
- Three tables hold what the code already did: `star-count` (the number of stars of a system by spectral class, moved out of `multistar`), `remnant-satellites` (a millisecond pulsar holds planets at 0.7%, the other remnants none; this replaces `tuning.SYSTEM_PASS_ODDS` for stage C) and `planet-satellites` (a rogue planet keeps a captured moon at 20%, a stand-in until GEN.466).
- The system pass, a rogue planet's moon and a multi-star group's star count read those tables. The draws are the same on the same streams: a test compares each against the formula it replaced over tens of thousands of seeds, and the stream is left where it was.
- Every host class carries a `may_host` gate (Boss's cards S5 and S6, 12:10 to 12:13 AM EDT, 2026-10-18, as the coordinator relayed, which superseded "moons are the lowest level"): a small body (comet, asteroid field) has the gate closed and no satellite table; a moon has a `moon-satellites` table with only the empty row until the tidal-survival row of GEN.579. A row states an orbiter's mass as a ratio to its host's where the source gives a ratio. The later recursion ends by the gate, not by a depth of three.
- Left for the item that builds them (U7 of `docs/design/unified-star-like-objects.md`): the star-like object's planets and belts, a brown dwarf's satellites and a planet's bound moons are still made by `StarSystem` and `generate_moons`.

- The scatter summary test (`test_the_scatter_lists_each_shells_phenomena_by_kind_like_the_star_scatter`) compares whole labels: a magnetar line contains the text of the plain neutron-star label and failed some seeds.

TODO items resolved: GEN.552 / Added: none
