**Headline:** The gas dwarf (class T) text now says its air can be thin, which is true of 1.2 percent of bodies; no bound was added, because it would leave them with no class.

### Changed

- GEN.1616: the phrase of class T (both rows, which share it) in `physics/data/planet_classes.json` no longer says "a thick hydrogen and helium envelope, with no surface to stand on"; it says "an air of mostly hydrogen and helium, thick on most" (the giant row of T shares the phrase, so both changed). The 11-16 census (default settings, 5,000 systems) has 8,458 terrestrial-kind T bodies; 851 of them (1.2 percent of all 73,247 bodies; all plain terrestrial bodies, median 24 K, median 1.1 kPa of pure hydrogen) hold a thin, cold hydrogen air with an envelope under 0.1 percent of the mass. Every other terrestrial row of `planet_class_fit.json` needs under 30 percent hydrogen, so an envelope or pressure bound on the T row would leave those 851 with no class at all; the sentence already said "the air can be thin", and the phrase now agrees. The basis line of the row and the row in `docs/design/habitability/condition-based-planet-classes.md` give the reason and the count.

generation-output: unchanged
generation-affects: none
