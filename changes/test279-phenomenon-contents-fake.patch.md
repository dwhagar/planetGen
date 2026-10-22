**Headline:** A phenomenon page opens even when the contents list cannot be loaded.

### Fixed
- **A nebula or supernova remnant page was a 404 when the web server answered its contents request with a 404** (TEST.279). The "What is inside it" list is extra, but only an unreachable API was tolerated; a not-found answer now leaves the block out too. The phenomenon page tests also fake the contents call, so they no longer depend on whether a web server is running on the machine.

### Changed
- `test_web_system_phen.py` and `test_class_pages.py` fake `get_phenomenon_contents`.
