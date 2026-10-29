**Headline:** The API names the life class of a planet or moon and the score of every class it supports.

### Changed
- **A planet or moon in `GET /api/systems/<id>` names its life by class.** `life_chemical` is gone (the life class switch replaced it); a body has `life_class`, `life_class_variant`, `life_class_best`, `life_class_best_score`, `life_animal_grade` and `life_complexity`, plus two fields made for readers: `life_class_id` (class and variant as one id such as `1a`) and `life_class_scores` (the score, 1 to 9, of every class the body supports). The API version is now 10.

api-version: changed
