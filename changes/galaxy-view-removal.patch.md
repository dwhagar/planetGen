### Fixed
- Removed the old `/galaxy/view` Galaxy Map endpoint (`html/galaxy_view.py`, `GET /api/galaxy/view`). The map has used cube tiles since 5.47.0, but a browser holding a cached copy of the old map script kept calling it, and its unbounded query timed out and tied up the API. Such a request now gets a quick 404 instead.
- The Galaxy Map page now loads its script as `static/galaxymap3d.js?v=<version>`, so each release reaches browsers on their next page view instead of a stale cached copy.
