**Headline:** The Sector Map's arrow-key test waits for the view to settle before it checks that another key does nothing.

### Fixed
- **`test_sector_map_arrow_keys_turn_the_view_and_keep_the_scale` failed in a long run of the map file.** After the last arrow it took its "before" picture while the turn and the nodes it asked for were still landing, so the picture after the `a` key differed for a reason that was not the key. It waits for the view to settle first (TEST.426). The window-resize test named in TEST.426 passed in the same 26-minute full-file run, and the remaining fail was this one.

TODO items resolved: TEST.426 (the resize test passed in a full-file run; the one failure seen was this arrow-key test).
