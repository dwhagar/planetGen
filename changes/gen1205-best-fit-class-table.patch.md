**Headline:** Every planet class now has its own full list of conditions, ready for a percent-match classifier.

### Added
- **The best-fit class table (GEN.1205).** `physics/data/planet_class_fit.json` gives every class letter its own positive conditions (F and G a warm, 0.5 to 2 g band; M an Earth-like air window; O and H ordinary air), with a shoulder width per quantity, and `physics/class_fit.py` scores a body against them. Nothing stored changes: the first-match rules still decide every letter. Shoulder widths are stand-ins CF1 to CF7 in GEN.725.
