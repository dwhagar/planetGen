**Headline:** The warm-up no longer reports an answer built on every run for a sector page whose Profile the API refuses.

### Fixed

- **The warm-up's "answers built" counts only answers that were built (MAP.289).** The sector page of a sector with a grid address asks the API for its Profile; where the API refuses (a database with no planned galaxy answers 404) the question was counted as an answer built, and asked again on every run, so a warm run over an unchanged galaxy always read "1 built" in its objects phase instead of "unchanged, nothing built". An answer the API refuses is now counted as neither built nor cached. A galaxy with a plan already cached its Profile and is unaffected; a test now holds that for several sectors. No stored value or generated output changes.
