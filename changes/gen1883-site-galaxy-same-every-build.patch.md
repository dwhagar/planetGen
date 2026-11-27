**Headline:** The sample galaxy the site tests walk is the same galaxy every time it is built

### Fixed

- GEN.1883: two builds of one commit made two sample galaxies (the heaviest planet was 006 in one and 002 in the next; a third build gave 003), so `test_irregular_moon_pages::test_the_system_detail_and_the_system_page_carry_the_irregular_moons_to_the_map`, which reads the heaviest planet, failed by chance. Two causes, both in the fixture `tests/test_web_a11y.py` builds (loaded by about twenty files through the dump cache): `plan` took a new galaxy seed from the operating system on each build (no `--seed`), and `fuzz_support.deterministic_entropy` is not enough to pin the run seed, because the run lines' tickets and ids draw from the same pinned stream a varying number of times, which moved the seed of the next command by a word between two builds. The build now gives `plan --seed` and pins `secrets.randbits` on a stream of its own. The new `tests/test_site_galaxy_same_twice.py` builds the sample twice in two scratch databases and holds the rows of the generated tables to be equal. The product is unchanged: the generator already takes a galaxy seed; only the fixture left it to chance. Tests only.

generation-output: unchanged
