**Headline:** Four test files find a system with comets in seconds instead of timing out in their setup.

### Fixed
- **`test_api_object_ids`, `test_parent_link`, `test_orbit_velocity` and `test_orbit_from_vector` ask for comets.** Since GEN.1775 a star holds comets only where its outer disk was scattered, so their fixtures' searches over 300 to 400 fixed streams for a system with comets ran past the 1,500 s test limit (`test_api_object_ids` failed with four setup errors in the Train 53 gate). Each fixture now sets `COMETS = True`, as `test_api_link_ids`, `test_body_positions` and `test_serialization` already do. Alone: 9 passed in 1m 51s, 6 in 9m 49s, 8 in 3m 56s and 7 in 3m 24s. Test-only.

generation-output: unchanged
