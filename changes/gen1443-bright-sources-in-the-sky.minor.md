**Headline:** A planet's or moon's page now names the other suns and the white dwarf in its sky, with their kind, colour, width and share of the light.

### Changed

- The description of a planet or moon gains a sky paragraph (`physics/sky_view.py`, `describe/build.py`): every other source of the system's radiation field (GEN.1441) that gives the body at least the light of a full Moon is named, with its kind (a white dwarf, a K-type star), its colour, its width in degrees and its share of the light. A source that gives at least the twilight of a day (`daylight.DAY_SUN_W_M2`, GEN.1461) is a second (third, ...) sun that makes a day of its own; a fainter one lights the nights, in full Moons. A body that circles a white dwarf is told so.
- A planet of a lone main-sequence star reads the same as before, and nothing is stored: the paragraph is worked out when the page is made, from the ledger of sources that every loaded system already builds. The light and the width of a source follow one distance, so the text agrees with the temperature (the one-companion light is `sky_suns.light_of`, shared with the suns of GEN.1458).
- The neutron stars, black holes with accretion disks and brown dwarfs of the orbit tree give no light in the field yet, so they are not in the sky text: that part of the item stays open (stand-in SKY3, section 6.253).

generation-output: unchanged
