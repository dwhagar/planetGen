**Headline:** the API's routes answer their 404s and 400s through `ApiError` and share one paging helper, with the same response bodies.

### Changed

- API.48: thirty routes built their own `jsonify({"error": ...}), 404` (or 400); they raise `ApiError` now, whose handler writes the same `{"error": message}` body and status. The two jump-by-value lists (uncharted systems, phenomena) shared the past-the-end answer and the `after` cursor parsing; both are one helper each (`_past_end_page`, `_after_number`), and `facility_list` reads its integer filters through `_int_arg` (same wording, blank means none) instead of a local function. Response bodies and statuses are unchanged.

generation-output: unchanged
