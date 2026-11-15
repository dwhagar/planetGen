**Headline:** A Class M abundance tier now makes desert worlds beside the temperate ones, and the meter counts all five habitable classes.

### Added
- **The dry part of the water dial (GEN.1721, first slice; Boss, 8:24 PM EST 2026-11-14: classes H, O, P and R present beside M, M in front).** Under a tier other than Science, one in five of the rocky planets the water dial reaches takes the dry water class (under 0.025 percent of its mass, the water of a desert world) and the rest take Earth's as before; the planet's air takes the row of the same class. The dial makes no draw of its own (it splits the water dial's one draw), so the same planets are reached. `tuning.ABUNDANCE_TIER_DIALS` rows are now `(water_share, bare_share, request_share, dry_part)`; registry row 6.289:AT5 (tier T9, tweak L4, a first-pass value for the calibration of GEN.1706).
- **The habitable family in the meter.** `scripts/class_m_meter.py report` prints, per tier, a table of classes M, H, O, P and R: worlds and systems holding one per 1,000 systems, the planets with each class's shape (its own limits of temperature, water cover and gravity) and how many of them pass its limits of pressure, oxygen and carbon dioxide. The limits are read from the class rules.
- **Design:** `docs/design/habitability/class-m-abundance-tiers.md` section 8 (what stops each of the five classes and the dials still to build: groundwater, the early warm phase and hardier oxygen life).

### Changed
- `delivered_water.drawn_for_earth` and the `earth` argument of `sample_fraction` are replaced by `dialed_class` and `dialed` (the class name); `body._water_dialed` holds the class name.

Science is untouched: every dial is 0 there, no draw is made and the golden systems are byte for byte the same (no generator epoch, no migration). A dry planet is Class H only with groundwater and oxygen, which later slices of GEN.1721 address, so this slice alone makes desert-shaped planets, not yet many Class H worlds.

Reuse: the water dial's draw and pin (`bulk_draw._rocky_water`, `mixture_subject`), the classes and rows of `physics/delivered_water.py`, `class_fit.overshoot` and `class_fit.CLASSES` for the meter's limits, and the census thread's `five_classes.py` for the table's shape.

generation-affects: none
