**Headline:** The admin edit tests post to the body's own page, and the phenomenon class-link tests stub the property sheet the page now reads.

### Fixed
- **`test_system_page_regenerates_a_planet` and `test_system_page_changes_a_class` failed on main (TEST.347).** UX.173 moved a planet's Admin menu from the system page to the planet's own page, so the system page now refuses a planet target ("That isn't something on this page") and has no Change class button. The tests now post to `/object/planet:<id>` and look for Change class there; the system page still has Change star.
- **`test_phenomenon_page_links_its_class` (five cases) and `test_phenomenon_unknown_class_stays_text` failed on main with a 502 (TEST.348).** The phenomenon page reads the object's property sheet (`get_object_sheet`), which the fake API of `test_class_pages.py` did not stub, so the page's own call went to a missing server. The fake now answers it with no groups.
TODO items resolved: TEST.347, TEST.348. Added: none.
