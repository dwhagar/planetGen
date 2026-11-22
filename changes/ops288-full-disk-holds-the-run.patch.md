**Headline:** A full disk now holds the run with a plain reason, "the disk is full: free space, then Retry held or finish-run", instead of ending it as a process that was gone.

### Changed

- OPS.288: the supervisor watches the free space of the three places a long run writes to (the database's data directory, the run record's Redis folder and the engine's state folder, one drive counted once) against a small reserve (512 MB, or half a percent of a bigger drive). Under two reserves nothing new of the bulk and upkeep lanes is admitted (the memory controller's own hold, reused); under one reserve the live generation run is held, its units stop at their next piece and keep what they committed, and the run says why. Nothing carries it on by itself: free space, then Retry held or `finish-run`. A run that ends with its process gone while a drive is full says the disk is full, not that its process was gone.

generation-output: unchanged
