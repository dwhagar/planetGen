**Headline:** The galaxy tile test expects the habitable-system `worlds` section, and the accessibility fixture names the response when a system read fails.

### Fixed
- **`test_galaxy_tiles_carry_no_comets_rogue_planets_or_asteroid_fields` failed every time on main (TEST.342).** MAP.220 added the `worlds` section (habitable-system markers) to every galaxy tile and the test's expected key set did not list it. It does now.
- **The browser fixture of `test_web_a11y.py` died with a bare `KeyError: 'planets'` when a system read did not return the system (TEST.329).** It now asserts the status of `/api/systems/<id>` and shows the response text, so a failed read names itself.
TODO items resolved: TEST.342. Added: none.
