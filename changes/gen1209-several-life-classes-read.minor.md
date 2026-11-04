**Headline:** The planet pages, the API, the search filter and the population pass now see every life class on a world, not only the dominant one.

### Added
- **`life_classes` on every planet and moon of the system detail (GEN.1209).** One entry for each class that lives on the body, the dominant one first, with its share of the food production, stage, animal grade and kind scores. A body stored before the per-class rows existed lists its one class. The existing `life_class_id` stays the dominant class.
- **The page names the other classes.** A body with more than one class gets an "Also living here" line in its Life section.

### Changed
- **The life class filter and its counts match a body when any class that lives on it is chosen,** not only the dominant one, and the map's life reasons do the same.
- **The class that built a civilization is the one that reached the technological milestone,** the first to reach it when several did, so a minor class can have the species on a world whose dominant class has not got that far.
