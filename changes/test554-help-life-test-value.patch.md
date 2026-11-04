**Headline:** The old-URL error tests know the new life help route, so the HEAD and OPTIONS checks run on it again.

### Changed
- **Test fix, no code change (TEST.554).** `test_web_old_urls_errors.py` builds a path for every page route from a table of argument values and stops on a route whose argument it does not know. `/help/life/<topic>` and `/help/life/classes/<int:number>` (DOC.84) had no `topic` or `number` value, so the HEAD and OPTIONS tests failed for the anonymous and the admin user. The table now has `topic` set to `abundance` and `number` to 1, both real pages of the life help, so the tests check those pages too.
generation-output: unchanged
