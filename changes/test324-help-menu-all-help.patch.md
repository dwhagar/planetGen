**Headline:** The top-bar menu test follows the Help menu's new first item, so it opens the Habitability page it means to.

### Changed

- `test_the_top_bar_menus` expected the Help menu to hold Classes, Habitability, Properties, Life and References, and clicked its second item. Since the help registry (the "All help" index first) the second item is Classes, so the click went to the wrong page and the test waited out its 120 second timeout (TEST.324). The test now expects "All help" first and clicks the item named Habitability. The file passes, 11 of 11, in about two minutes.
