**Headline:** The retired `/<name>.py` page addresses and the old NAV parameter style now answer 404 instead of redirecting.

### Removed
- The `/<name>.py` redirects (`index.py`, `search.py`, `sector.py` and the rest) and `web/old_urls.py`; those addresses now give the normal 404 page.
- The old NAV parameters (`from_id`, `to_id`, `from_kind`, `to_kind`, `from_type`, `to_type`), which used to redirect to the current `from`/`to` form.
- The redirect tests for both, about 190 cases; the 404, header, malformed-body, HEAD and OPTIONS sweeps stay.
