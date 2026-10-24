**Headline:** The hosted nebula row tests no longer fail when the built star draws a heavier companion.

### Fixed
- **`test_hosted_nebula_rows` failed with all data sets fetched (TEST.311).** The test built the row's star as a plain `StarSystem`, which draws a member tree of its own (GEN.631), and a heavier companion drawn there became `stars[0]` and lit, or did not light, a nebula the row does not predict. A sector fill never draws a tree for a row with no stored members: it passes `stored_trees.SINGLE`. The test now builds the star the same way.
TODO items resolved: TEST.311. Added: none.
