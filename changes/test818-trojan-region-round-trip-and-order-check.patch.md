**Headline:** A system read back from its JSON keeps its trojan regions as trojan regions, and the fuzz test's orbit checks know a trojan is no slot of its own.

### Fixed
- TEST.818: `test_every_star_type_generates_a_sane_system` failed on the hypergiants (IA, IB) with "bodies not in orbital order". Since GEN.1150 a giant's trojan regions stand at the end of its star's list at the giant's own distance, by design (`validation.slot_planets` leaves them out of the spacing), but `assert_orbits_sane` still read the whole list as slots. It now checks each trojan (region or planet) against its giant (an L4 or L5 point, a giant of the list, the same distance) and leaves it out of the order and spacing checks.
- Behind that check the same test found a real bug: `StarSystem.from_dict` rebuilt a trojan region as a plain belt, so a system read back from its JSON lost the region's Lagrange point, libration, mean eccentricity and inclination (the store's loader was right). It now builds a `TrojanRegion` and links each trojan to the giant at its distance.

generation-output: unchanged
