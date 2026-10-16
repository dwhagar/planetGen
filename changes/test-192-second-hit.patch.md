**Headline:** The TODO list records a second sighting of the galaxy-grid fuzz failure and a suggested fix.

### Changed
- **TEST.192 gains Core 3's finding.** The same fuzz failure came up once in about 17,900 tests, and the item now carries the cause (the rounding error of the halo subtraction is larger than the tolerance on a bound near 2.7e-20) and a suggested tolerance relative to the model.

TODO items resolved: none / Added: none
