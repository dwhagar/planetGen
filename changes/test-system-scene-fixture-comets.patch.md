**Headline:** The system-scene tests find a system with comets in minutes instead of hours.

### Fixed
- **`test_system_scene` no longer searches for hours for a system with comets.** Under the default science comet tier (GEN.1765, GEN.1775) few stars hold comets, so its fixture's search over 400 fixed streams, run again for each of its 9 tests, took longer than the 6,000 s test limit. The fixture now makes its systems under the abundant comet tier and makes them once for the whole run (12m 31s alone for all 9 tests). Test-only; `generation-output: unchanged`.

generation-output: unchanged
