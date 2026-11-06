**Headline:** The page of an uncharted sector now lists what is placed in it, each object linked to its own page, and says plainly that the sector is uncharted.

### Added
- The uncharted-sector page (`/sector/at/<shell>/<band>/<wedge>`, and the page of a placeholder sector row from the Sectors list) has a "What is placed here" panel: a table of the rogue planets and other objects stored in the cell, the bright stars the scatters left waiting in it, the unbuilt scattered phenomena (black holes, neutron stars, planetary nebulae, supernova remnants, hypervelocity stars) and the nebulae its waiting stars light. Each row links to the object's page (`/phenomenon/...`, `/systems/uncharted/...`, `/uncharted-phenomenon/...`). A cell with nothing in it says so. Nothing is generated; the lists are read from what is stored.
- New tests: `test_web_uncharted_contents.py` (the page with a faked API) and a browser check in `test_web_uncharted_sectors_browser.py` (the panel on a placeholder's page fits the phone, tablet and desktop and passes the accessibility check).

generation-output: unchanged
