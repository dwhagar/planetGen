**Headline:** The "/9" beside a life score is now easy enough to read on white, and the accessibility test waits for the Generate system menu to finish building.

### Fixed
- **The "/9" after a life score (the habitability tree and system page) had a contrast of 3.18 to 1 on white**, under the 4.5 to 1 the accessibility test requires (UX.229). It was the score's colour faded to 70% opacity; it now uses the page's muted text colour, which passes in both themes. Eight `test_web_a11y` cases on the system and object pages passed again.
- **`test_web_a11y` could check the Generate system page before its menu had built** (UX.230). Axe then saw the options without their list around them ("aria-required-parent", once, at phone/light). The test now waits until every Shoelace element on the page is defined and has rendered before it runs axe. The page itself is unchanged: the menu is built by Shoelace and passed in the other three cases.
TODO items resolved: UX.229, UX.230 (TODO closes them after reading this). Added: none.
