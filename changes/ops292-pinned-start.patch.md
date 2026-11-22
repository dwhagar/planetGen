**Headline:** A random start now says the arguments that start a run at the same sector again, so two full runs of one seed can be compared.

### Changed

- OPS.292: the line that says a random starting sector as it is saved is followed by "To start another run at this same sector, give its sectors entry: --shell S --band B --wedge W --radius-pc R". A run spec that gives those arguments in place of an empty `sectors` entry starts at the same sector, charts the same neighbors and has the same 19 stages (the existing single-address mode; nothing new in the engine). A new test runs a pinned run whole and a pinned run stopped partway and carried on by `finish-run`, and checks they end with the same galaxy.

generation-output: unchanged
