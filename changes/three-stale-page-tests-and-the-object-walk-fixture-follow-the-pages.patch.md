**Headline:** Three tests that still described pages as they were before UX.173 and the Radiation table now match the pages, and the object page walk makes far fewer API calls.

### Fixed
- **`test_api.py` description test expected one facts table per planet (TEST.354).** A planet's description now holds its Planet table, a Radiation table and, for a scored rocky body, a Survival equipment table. The test now asks for exactly one Planet table among those three captions.
- **`test_class_pages.py` nebula D page test crashed on a nebula with no type (TEST.355).** Every stored nebula has a `nebula_type` (its class's family); the test's thin detail had none, so the builder called `.capitalize()` on `None`. The detail now carries the family of class D.
- **`test_web_search.py` star result test looked for `/system/7#star-3` (TEST.356).** A body result now opens the body's own page (`/object/star:3`, UX.173 and UX.176), so the test expects that and renamed to say so.
- **`test_object_page_walk.py` called the API twice for each of the sample galaxy's 45 systems (TEST.349 follow-up).** The walk only opens six objects of a kind, so its fixture now stops once every kind has six (the facilities call is skipped once six are found). That is far below the API's 50-an-hour limit that caused the TEST.349 429, and the module runs in 46 s instead of 125 s.
TODO items resolved: TEST.354, TEST.355, TEST.356. Added: none.
