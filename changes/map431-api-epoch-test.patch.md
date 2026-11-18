generation-output: unchanged

MAP.431 (bug): the `test_api` test that a scatter drawn again makes every map node stale now names an epoch that changed the scatter; the epoch it made up changed nothing, which the map now reads as no change.
