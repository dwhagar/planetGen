**Headline:** The test of the black-widow and redback spiders now expects the post-common-envelope history the generator gives them, which is the physics.

### Changed
- **`test_accreting_compact_binaries` spider case (TEST.539).** It expected `detached` for the pulsar-and-companion pair of a spider, a label from before the pair history rule (GEN.1012, `pair_history.must_be_post_common_envelope`). A pair with a period of hours cannot have held the pulsar's giant progenitor, whose largest radius is far beyond the whole orbit, so the generator is right to call it `post-common-envelope` (the usual route to a low-mass X-ray binary and then a spider); the test is the side that was wrong. The generator is unchanged.

No generator output changes.

generation-output: unchanged

TODO items resolved: TEST.539.
