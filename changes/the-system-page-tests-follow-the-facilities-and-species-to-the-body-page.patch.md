**Headline:** The facility and population system page tests follow what UX.173 and UX.176 moved to the body's own page.

### Fixed
- **Eighteen tests of `test_web_facilities.py` and two of `test_web_population.py` failed on main (TEST.350).** Their fake API did not answer the calls the system page's tree now makes (`get_object_children`, `get_object_sheet`, and for the facilities page `get_system_habitable_species` and `get_system_text`), so the page answered 502. The fakes now answer them.
- **Four of those tests looked for rows UX.176 removed from the system page.** The per-body facility list, the Inhabited chip and the dominant species are on the body's own page. The facility and colony tests now read the planet's page for them; the system page keeps the Facilities table and the territory. A new test in `test_item_pages.py` covers the dominant species on a life world's page (a lifeless body is never looked up).
TODO items resolved: none (TEST.350 stays open for its browser case). Added: none.
