**Headline:** a star read back from the database keeps its state word, so a sky text names a young or evolved companion the way the generated star did.

### Changed

- TEST.728: `stars.state` and `stars.decay_age_yr` are stored columns but not serialized fields, so a star rebuilt from its row had neither: a planet or moon whose sky holds a stored pre-main-sequence or evolved star said "an M-type star" where the generated system said "a young star". The rebuilt star now carries both. The describe and system-render tests took their seeds from `hash()`, which differs per process, so the mismatch showed in some gate runs only; they now seed from a fixed value (`zlib.crc32`).
