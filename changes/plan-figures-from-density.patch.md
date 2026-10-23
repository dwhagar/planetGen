**Headline:** The plan's crest and planned-share figures now come from the galaxy's own density and follow the Census dial.

### Fixed

- The plan's "Outline" line printed a crest of 11.5 times the one-star density where the model's own thin-disc crest was 8.76, and the same 11.5 with `--richness 1.0` where the model gives 4.37. The crest is now the number the outline walks on (`k_norm` x arm level x (1 + arm contrast) x the systems a sector holds at the real local density, `galaxy/plan_figures.py`), so it follows the galaxy's own system mass and the Census dial.
- The share of the stars in planned sectors (38 percent against the model's 59 for the same galaxy) was a table keyed on the thin crest alone and measured on the Milky Way's structure. It is now one integral of the built density over radius, height and azimuth, so a thick disc, bulge, nuclear disc or stellar halo of any size is in it. The table, the fixed 0.483 solar masses a system and the thin-disc surface-density shortcut are removed. The warning for a crest under 20 and the refusal of an empty outline are unchanged.

TODO items resolved / Added: GEN.626 / Added: none
