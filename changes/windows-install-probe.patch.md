### Fixed
- The Windows installer no longer stops when planetGen isn't installed yet: the check for an existing editable install asks Python where the package is without importing it, so it prints no traceback.
