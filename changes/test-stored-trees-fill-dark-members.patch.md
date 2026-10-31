**Headline:** A test of the stored-tree fill now expects a brown dwarf or remnant in a tree to be a dark member, not a star row.

### Changed
- The stored-trees fill test (`test_gen_stored_trees.py`) counted every member of a stored tree as a star row of the built system, so it failed on any tree holding a brown dwarf, a neutron star or a black hole since GEN.1013 keeps those as dark members. It now expects the star rows to be the lit members and the `dark_members` rows to be the dark ones. Test only: no code, stored value or page changes.
generation-output: unchanged
