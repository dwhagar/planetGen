**Headline:** A world's Radiation rating now includes its ground ultraviolet, so an oxygen world with no ozone layer or a very bright star no longer rates Blue.

### Changed
- **GEN.445: ultraviolet counts in the Radiation domain.** The score is the lower of the dose's and the ultraviolet's: Blue under 3 times Earth's DNA-damaging UV, Green to 30, Yellow to 300 (the unfiltered light of an ozone-free world), Red beyond (`uv_lethal_index` and `uv_outer_index` are new engineering rows in the limits table; the levels page shows the new row). A world without a UV figure scores on its dose alone.
- **One cosmic-ray curve.** The second, hand-fitted cosmic-ray curve in `habitability.py` is removed; a world with no stored dose reads the tables of the radiation model (cosmic rays only, no ground), and a test fails if the old curve returns.
- **Generator epoch.** Stored habitability scores, tiers and kit can change for worlds with strong ultraviolet; no draw moves. The golden record is regenerated for it.
generation-output: changed
