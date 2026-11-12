**Headline:** The HEAD and OPTIONS check of every page route has a test value for the irregular-moons page, which turns main green again.

### Changed

- `test_web_old_urls_errors.py` gives `planet_ref` (`planet:1`) in the values it fills page routes with, so `/object/<planet_ref>/irregular-moons` and `/object/<planet_ref>/irregular/<rank>` (GEN.1086) are checked like the belt and field pages.

generation-output: unchanged
