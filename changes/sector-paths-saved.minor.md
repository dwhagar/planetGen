### Added
- The path each star system, rogue planet and interstellar comet takes through its sector is saved (`sector_paths` and `sector_path_knots`, schema v62): worked out at each orbit update for the sectors something moved into or out of (not at generation, where they would depend on the order workers fill neighbouring sectors). Run `sudo ./update.sh` to migrate.
