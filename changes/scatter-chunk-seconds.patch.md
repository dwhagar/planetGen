### Fixed
- The phenomena scatter at default scale no longer sits silent for minutes before its first progress: layers go to the workers in chunks sized by how long they take (about 2 seconds each), so a layer that takes seconds is a task of its own and shows at once. The first progress now comes in about 8 seconds, with two workers or one.
