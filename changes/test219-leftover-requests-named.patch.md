**Headline:** The Generate page's line under the dials is worked out once for the same dials, and the test server's warning about requests left running names each route.

### Changed
- **TEST.219: what the browser pages left running after the server stopped (kind: perf and test).** The TEST.193 warning now names every request still running by method and path. In `test_web_generate_alignment.py` and `test_web_browser_controls.py` they were `GET /admin/generate/shape-report` (10 at once in one page, each working the whole plan outline out: 3 to 10 s warm, 12 to 116 s cold on a busy machine) and a few Shoelace files under `/static/` that a closing page had stopped reading. The shape report is now kept (`generate_page._shape_preview`, 32 answers): the same dials are answered at once and a repeat of dials still being worked out waits for the one answer instead of making its own. Requests under `/static/` hold no database connection, so the server stop no longer waits for them or counts them. Nothing the generator makes changes.

generation-output: unchanged

TODO items resolved: TEST.219.
