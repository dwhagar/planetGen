### Fixed
- **install.ps1 and update.ps1 no longer fail when no Redis answers**:
  the Redis check stays a warning instead of leaving its exit code as
  the script's.
