**Headline:** Two reds on main after Train 24 are fixed in the tests, with no change to what the generator makes.

### Changed
- **The materials reachability scan reads the evaporite sequence.** `physics/evaporite_sequence.py` names the 19 salts, carbonates and perchlorates it dries out of a brine (GEN.1357 and GEN.1366), but it takes its molar masses from its own formula table and does not import the materials table, so the scan skipped it and called the rows unread. It is now in the scan's always-read list.
- **The rogue and bound planet parity test allows the air quantities to differ where only one of the pair holds an air.** A rogue world kept warm by its own heat holds an air; the matched planet at 3 K far from a star has frozen its gases out and stores no density, molar mass or scale height. The test already allowed the same case for its liquid surface and clouds.

generation-affects: none
