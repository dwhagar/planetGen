**Headline:** The field-cloud placement tests use a cloud that fills its sphere, so they no longer fail when a drawn nebula's lumpy shape leaves its centre empty.

### Changed

- Three tests in `test_gen_cloud_at_generation.py` put a system at the centre of a drawn nebula and expected the cloud to hold it. A nebula's shape (GEN.75) is a few overlapping blobs and need not hold its own centre; the nebulas those tests draw changed with the generator, and the drawn shape of the test's nebula left the centre empty, so `surrounding_cloud` was `None`. The code was right and the test's assumption was not. The test's nebula now gets a shape that fills its whole bounding sphere, which keeps what the tests check (a placed system finds the cloud, takes the smaller of two nested clouds, turns the sector's local axes into the galaxy frame, and draws the same enrichment for the same placement). A drawn shape's own containment stays tested in `test_nebula_shape.py`. Tests only.

generation-output: unchanged
