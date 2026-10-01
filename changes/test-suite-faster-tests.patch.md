### Changed
- **The slowest tests run in seconds, and two tests that skipped at
  random now always run.** The radius-neighborhood generation test builds
  a sparser neighborhood (70 s to 3 s) and also checks every generated
  sector is inside the radius; the sector-enumeration fuzz test checks
  tolerances with `math.isclose` (30 s to 2 s). Tests hash passwords with
  1,000 PBKDF2 rounds instead of 600,000, except the tests of the hashing
  setting itself (marked `real_password_hashing`). The moon life-data
  round trip used to skip almost every run (a moon with its own life data
  is too rare to wait for) and the black-hole test skipped on one draw in
  ten; both are now seeded and always run. The full suite on 4 workers
  went from 5 min 13 s to 4 min 27 s.
