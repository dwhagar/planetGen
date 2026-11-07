**Headline:** The API route snapshot lists the Queue page's held-jobs route, so the API version test passes again.

### Changed

- **API snapshot gains `/api/admin/work/held GET` (TEST.675).** The Queue page listing of the jobs held for a job they require (PERF.248) added the route without rewriting `tests/fixtures/api_routes.json`, so `test_api_version` failed on main. The route is new and nothing existing changed, so it is not a breaking change: `API_VERSION` stays 13 and `docs/reference/api.md` gets no line. `/api/uncharted-systems/<id>/generate POST` is already in the snapshot. Test fixture only.

generation-output: unchanged
