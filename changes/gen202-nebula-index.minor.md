**Headline:** One read now lists every kind of nebula in a box of the galaxy (stored nebulae and remnants, the seed's field clouds, the scattered planetary nebulae and remnants) with no object twice.

### Added
- GEN.202: `db.nebula_index.nebulae_in_box` reads the stored nebulae and supernova remnants of a box, the field clouds no sector has stored and the scattered planetary nebulae and remnants no sector has built, and lists each object once by its Object ID (the stored one first, then the field cloud, then the scattered mark), each entry naming its `source`. The three readers already left out what is stored by their own flags; the index also drops a row that holds the Object ID of a stored one whatever its flag says. A test plants all three kinds with a repeat and shows each object once. Nothing generated changes.

generation-output: unchanged
generation-affects: none
