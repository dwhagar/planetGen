**Headline:** A database that has no galaxy yet now says so on every page, instead of "no such database".

### Changed
- **A database that holds only `migrate`'s bookkeeping (or nothing) says there is no galaxy yet (MAP.407).** Every page used to answer 404 "Not found: No such database". The API now marks that case (`noGalaxy`, still a 404) and the pages show "No galaxy yet: There is no galaxy in database ... yet. Plan one from the admin Generate page (or run 'planetgen plan') and the pages fill in." A schema that holds something else (the control database, for one) is still "No such database". The web client's new `NoGalaxyError` is a `NotFoundError`, so every page that handles a missing object handles it as before.
- **MAP.406 needed no change.** The sector page and `/sector/<id>/gravity` with no galaxy shape were fixed by MAP.405 (the bug hunt ran on a main from before it); all of the item's red tests (the sector pages, profile, Nearby, route fuzz) pass on main.
