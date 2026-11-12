**Headline:** The reds on main after Train 24 in the materials and rogue planet tests are fixed in the tests, with no change to what the generator makes.

### Changed
- **The materials reachability scan reads the evaporite sequence.** `physics/evaporite_sequence.py` names the 19 salts, carbonates and perchlorates it dries out of a brine (GEN.1357 and GEN.1366), but it takes its molar masses from its own formula table and does not import the materials table, so the scan skipped it and called the rows unread. It is now in the scan's always-read list.
- **The rogue and bound planet parity test allows the air quantities to differ where only one of the pair holds an air.** A rogue world kept warm by its own heat holds an air; the matched planet at 3 K far from a star has frozen its gases out and stores no density, molar mass or scale height. The test already allowed the same case for its liquid surface and clouds.
- **The rogue parity test allows `cloud_chemistry` to differ with the cloud.** The sulfur-cloud chemistry (GEN.1425) is read off the acid cloud the air holds, so it is filled on the one of a matched pair whose air holds one.
- **The rogue and bound moon-code test finds a pair that made two moons each instead of pinning seed 0.** `generate_moons` attempts the count and makes fewer where the orbits have no room, so the test searches seeds (TEST.714) for a pair with two each.
- **The stored UV index test skips a planet under a deep envelope.** Such a planet has no ground (GEN.1334), so it stores no dose and no UV index; the test now says so for it and checks the planets that have a ground. The seed's system holds two such class T planets after Train 24.

generation-affects: none
