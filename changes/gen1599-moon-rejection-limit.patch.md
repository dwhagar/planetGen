**Headline:** Moon generation stops a full disk at the call's own rejection limit, the same limit the too-large moons use.

### Changed

- GEN.1599: `physics/planets.generate_moons` counted the disk-full rejections against the module's `MOON_REJECTION_LIMIT` and the too-large ones against the call's `rejection_limit`. Both now stop at the call's `rejection_limit`. Nothing stored changes: a disk exists only round a giant, whose limit is the same number. A test pins the one limit.

generation-output: unchanged
generation-affects: none
