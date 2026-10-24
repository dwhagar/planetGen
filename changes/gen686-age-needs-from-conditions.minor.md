**Headline:** How old a star must be to host a planet now depends on the planet's air, temperature and gravity, not on its class letter.

### Changed
- **Age needs from conditions (GEN.686).** The youngest a star can be and still host a planet or moon is read from the body's stored conditions (`physics/age_needs.py`, `physics/data/age_needs.json`) as one of five tiers (no settling time for giants, airless and small worlds; molten; settled; watered; frozen), each with its minimum age on a fast, normal or slow evolving star. The ages are the old per-class minimums, so about 99% of bodies need the same age as before and a few (some thin-air and cold worlds) need a different one. The `age_ranges` of every class in `tuning.PLANET_CLASSES` are deleted.
generation-output: changed
