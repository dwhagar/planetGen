**Headline:** Every composition category becomes a gated table of weighted possibilities, filed in the plan.

### Changed

Boss's design requirement for the composition work is recorded in `docs/design/composition-model.md` section 4b: each composition category is a table of possibilities, each row lists the planet or object classes it supports (a class gate), and the draw takes the rows that pass the gate with their weights normalised to 100 percent. The standing columns are possibility, supported classes, weight or weight function, source and licence, and notes. The plan gains the one loader and format (GEN.472) and a test that the gated weights of every registered table add up to 100 percent (TEST.214), and the table and recipe items (GEN.410, GEN.411, GEN.427, GEN.428, GEN.439, GEN.440, GEN.441, GEN.442, GEN.454 to GEN.457, GEN.466) now point at them.

TODO items resolved: none / Added: GEN.472, TEST.214
