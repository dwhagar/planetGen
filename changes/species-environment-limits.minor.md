### Added
- Every species now has environment limits: the human habitability table with some rows replaced, made from its homeworld so a species always passes its own homeworld. The six habitability domains can be scored against a species' limits, and a species' limits show on its page (POP.12).
- Lists of what could live where: the species that could live on a planet or moon (unaided, or with equipment their tech level can make), the worlds of a system or sector one species could live on, a "Species that could live here" panel on the system page, and API routes for each. The limits are generated from the homeworld until DB.30 adds the columns to store them.

### Fixed
- The system page no longer answers 404 when the API says there is no owner, facility types or species list to show; its tests now fake those calls.
