**Headline:** a requested planet class, the count of Class M worlds and the search now read how well a body matches its class, so a poor fit of a class is shown with its match but does not fill a request for it.

### Changed

- `tuning.CLASS_MATCH_THRESHOLD` 0.5 and `physics/class_outcome.py` (GEN.1315, new): `fits` and `has_class` say whether a stored `class_match` is high enough for the class to count (registry row FT5; the suggestion of the fudge table plan).
- `generation/system.py`: a request for a class (`build_planet`, the forced-class check) is met only by a body with that class at a match that counts, and `count_habitable` counts Class M worlds the same way. The habitable count and the bar on habitable worlds are not changed: GEN.1387 and GEN.1397 already moved them from the class letter to the PHI-6 score and the life of the finished body, so they have no letter left to read.
- `db/search_values.py`: the search takes a range on the match (`planet_min_class_match=0.5`, `moon_max_class_match=0.9`), so a person can ask for worlds that fit their class fully.
- Docs: the condition-based classes design note says what reads the match. A request for Class M is as often unmet as before (the class is rare, GEN.1403, GEN.1554); the threshold changes no census count by itself.

generation-affects: fill
