**Headline:** A world's air is now the sum of the life classes that live on it, with the oxygen rise counting the methane its earlier classes made.

### Changed
- **The air of a living world is the sum of its classes (GEN.1213).** `physics/life_air.Living` gives each class its share of the food production; the methane of the classes adds up (the carbon they use counted once), the oxygen source goes with the summed share of the oxygenic classes, and the oxygen rise is where their surplus has filled the store of reduced things (`biogenic_oxygen.rise_age_gyr`). That store now includes the methane the classes before the rise hold, which the single-class air left out. A world with one class gives the same air as before except for the rise age, which moves by a few parts in 1e4 of the delay.
- **The Life and air chart steps at each class's start.** `GET /api/planets/<id>/life-air` and the moon route add `classes` (class, variant, share, producers age).
