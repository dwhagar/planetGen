**Headline:** The life and timeline paragraphs of a world are written from its stored row alone, and an edit keeps the life a world already has.

### Changed
- **The milestone schedule is stored (DB.100, schema v156).** `planets` and `moons` keep the age at which the world's pace reaches each of the five milestones (`life_age_abiogenesis_gyr` to `life_age_technological_civilization_gyr`), set with the pace by `life.apply_life_data` (`evolution.schedule`). A galaxy made earlier gets the ages of its stored pace (migration 0153); a body with no life chemistry has none. The property sheet lists them under Life.
- **The page is built from the stored row.** The words of a pace moved out of the milestone table (`tuning.EVOLUTIONARY_PACE_LABELS`), so `describe/life_facts.py` and `describe/build.py` read the row and key lookups only, and a test renders a saved world with no draw and no pace table.
- **An edit keeps the life a world has.** `validation.reapply_life` (the admin edits that change a class or the climate) keeps the stored chemistry and pace when the new climate and the star still allow them; only one that no longer fits is drawn again. A body the habitability index gives no microbial life no longer keeps a pace an earlier pass drew, so a generated body has a pace only with a chemistry and the generator epoch is bumped.

generation-output: changed
