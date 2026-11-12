**Headline:** a body that stores a NaN quantity can no longer be fitted to a class as a perfect match: the class fit scores a NaN as no match and the classifier stops with a plain error that names the quantity.

### Changed

- `class_fit.membership` reads a NaN quantity as 0, as a missing one already was (every comparison with NaN is false, so a NaN used to pass both a minimum and an upper line). A NaN pressure is on neither side of the air gate.
- `class_fit.classify` (the path every build and census takes) raises `planet_classes.ClassError` naming the quantities that are not numbers, instead of choosing a best fit for them. No body stores a NaN today (none in 850 sampled bodies), so no generated output changes.
