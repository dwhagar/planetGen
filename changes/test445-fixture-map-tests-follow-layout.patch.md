**Headline:** Three browser map tests follow the layout changes of the zoom redesign and the contents menu.

### Changed
- **The point phenomena test (TEST.445).** A step in is as far as the nearest thing along the ray allows (MAP.334), so the fixture's black hole, pulsar and quasar still stood within a few pixels of one another after one; the test steps in until each is 40 px from the next before it clicks them.
- **The nebula click test (TEST.445).** The contents list under the title (UX.215) takes height from the map in a 1000x800 window, and no spot on the shorter map named the nebula; the test runs in 1000x900.
- **The map card and title test (TEST.445).** A page with a contents menu is wider and puts its title beside the menu (UX.215), so the Galaxy page's title (which has no menu since TEST.442) is compared with the Home page's, not the Sector page's.
