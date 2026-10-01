### Added
- **One module to validate a planet, a lunar system and a star system
  (ADM.5).** `stellarObjects/validation.py` holds checks that report what
  is wrong without changing anything, the orbit-spacing pass generation
  already ran (moved there unchanged), and a stabilize pass that re-spaces
  an edited system from its moons outward, ready for the admin overrides
  (ADM.6, ADM.7).
