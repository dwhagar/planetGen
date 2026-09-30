### Changed

- The Galaxy Map sizes its sector blocks from the screen: each block is the smallest power-of-3 cube of whole sectors that is at least 4 pixels across at the focus, and blocks line up with the sector grid's master wedges wherever that keeps them near one block long. Only the solid's visible surface is built, so views zoom in to single sectors much sooner, and a clicked block shows its exact sector count, leaving out sectors the galaxy's outline doesn't allow.
- The Galaxy Map's wedge lines now follow every master wedge (3 from the core, doubling outward), each zone shown once its lines are far enough apart on screen.
