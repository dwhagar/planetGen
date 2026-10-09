### Added
- The path each star system, rogue planet and interstellar comet takes through its sector is saved (`sector_paths` and `sector_path_knots`, schema v62): worked out when a sector is generated and recomputed at each orbit update for the sectors something moved into or out of. Run `sudo ./update.sh` to migrate; existing sectors get their paths at their next orbit update that touches them.
