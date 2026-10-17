**Headline:** The scatter layers draw their random numbers in bulk, with the same numbers and the same stream state as before.

### Changed
- **Bulk random draws for the scatter layers (PERF.91).** `Stream.random_array` and `Stream.randbits_list` hand a stream's state to NumPy's `MT19937`, which is the same generator with the same 53-bit float, draw a batch in C and hand the state back, from 1,000 draws on. The phenomenon scatter, the far-reach scatter and the component pass use them for their per-object uniforms and seeds (12 times faster on batches of 60,000), and a test proves same seed, same floats, same bits and the same final state, so no galaxy changes and nothing is reseeded.
