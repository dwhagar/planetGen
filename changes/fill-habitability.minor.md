### Added

- **`planetgen plan --fill-habitability` refreshes the stored habitability numbers of an existing galaxy (GEN.217).** It re-derives every planet's and moon's flux swing and score (the PHI-6 index, the six domains, the equipment tier and the note) from its stored conditions, so the planet chip, the search page and the facet counts agree with the description after the rules change. Nothing else is touched, and it can be stopped and run again.

### Fixed

- **A lifeless moon stored the habitability score it was first given, before its temperature settled (GEN.217).** The chip and the search page could disagree with the description, which is worked out from the final conditions. A new moon's score is now redone once everything is final.
