**Headline:** A system of several stars is no longer drawn twice to check that its stars clear their orbits, because the tree and the build now read the same stars.

### Changed
- **The check that redrew a whole set of companions when a built star was wider than its orbit (GEN.1015).** `multistar.built_tree_clears`, the clearance branch of `StarSystem._companions_fit`, `StarSystem._built_stars_clear` and the redraw kind `companions.overlap` are gone, with no shim. Each companion is built from the record the tree drew it with (GEN.1009), so a built star is the star its orbit was drawn around and the check never failed; the independent verifier (`physics/clearance.py`) is what the tests now ask. Generated systems are unchanged.
