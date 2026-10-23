**Headline:** The generator now has a battery digest: seven tiny canned cases that say whether it still makes what it made, and which case differs when it does not.

### Changed
- **`generator_battery_digest()` (OPS.57).** A function in `planetgen.generation.generator_battery` with no database, no network and no clock: two three-system sectors, a rogue planet, a nebula shape, a bright-star scatter on a tiny outline, a population stream and a known unit seed, each hashed with the GEN.304 encoder into a full, a structural and a float digest, with an overall digest over the lot. It is what OPS.15, OPS.14's mismatch check and TEST.77's first tier compare; nothing calls it yet.
- `tests/test_generator_battery.py`: two runs and two hash seeds give the same digests; changing one canned draw changes exactly one case; a lone float change moves the float digest and leaves the structural one. The design note (`docs/design/generation-determinism.md` section 5.4) records how it was built.

TODO items resolved: OPS.57. Added: none.
