**Headline:** A radiation dose now picks its unit by size, so 13,200 mSv/yr reads "13.2 Sv/yr" and not "1.32 × 10⁴ mSv/yr".

### Changed
- **A dose rate has a unit ladder.** The unit registry gained the `dose_rate` kind (nSv/yr, µSv/yr, mSv/yr, Sv/yr, kSv/yr, MSv/yr) and
  `util/ladders.py` the `DOSE_RATE` ladder (also in `unitladders.js` and `unitregistry.js`), so a property sheet's dose rows and a
  Document's dose facts show the largest unit the value is at least 1 of. The sheet had written the number with a fixed mSv/yr and a
  `.4g` format, which turned any dose of 10,000 mSv/yr or more into scientific notation. The stored columns stay in mSv/yr (no migration, no
  epoch); the Document's `dose` kind and its factor of 1000 are gone, in favor of the ladder's `dose rate`.

generation-affects: none
