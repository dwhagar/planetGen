**Headline:** Dead code left by earlier moves is gone.

### Fixed
- **Dead code left by earlier moves is gone.** `web/lib/xfetch.py` still held copies of `is_cacheable` and `shared_key` that named a constant and a module it never imported (the live ones are in `pagecache.py`), and `cli/orbits.py` defined `_stage_label` twice.

TODO items resolved: none. Added: none.
