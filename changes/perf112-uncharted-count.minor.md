**Headline:** The Uncharted systems table no longer counts every row of the galaxy to show its total.

### Changed
- **`count_uncharted_systems` answers from the stored counts (PERF.112).** The total behind the Uncharted systems table and `/api/systems/uncharted` went through a full `COUNT(*)` over the stars-without-a-system join on every request. It now uses the same stored count as the other tables (`db/countcache.py`): the last answer is served at once while a background recount refreshes it, and a galaxy with no stored answer yet gets a capped estimate (`total_capped` is true in the API and the page shows "more than N"). The count is stamped with the galaxy's content stamp, so it is recounted when the galaxy changes.

TODO items resolved: PERF.112. Added: none.
