**Headline:** Rogue planets now come from the systems that threw them out: each is built like any other world, with no star, and filed in the sector it has drifted to.

### Changed
- **A rogue planet row has every planet column.** Surface, air, interior, ocean, life and the rest sit on the rogue row the way they do on a planet row, with the same child tables for its layers, gases, spectrum and life scores. Rogue planets already in the database are not changed; the regenerated galaxy replaces them.
- **A system that threw a body out saves it.** When a system is saved, each body its planet architecture ejected is built as a planet with no star (its own heat warming it, its air kept from the ejection), moved along its kick for the time since, and stored as a rogue planet in the sector where it landed. The rogue keeps the system that threw it out, which of that system's ejections it is, how long ago and how fast, and a sentence that says where it came from. Saving a system again replaces the rogues it made, so none are doubled.
- **A rogue's sector may be one nobody has charted.** The sector gets an empty row, not in the chart ledger, which is filled when the cell is charted. A rogue's ID belongs to the sector of the system that threw it out, so deleting that sector removes its rogues wherever they are, and a sector holding a rogue from elsewhere cannot be deleted.

generation-output: changed
