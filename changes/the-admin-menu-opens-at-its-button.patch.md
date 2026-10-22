**Headline:** The Admin menu on the sector, system and phenomenon pages now opens under its button instead of far above it, off the screen.

### Fixed
- **The Admin menu (and the Navigate menu) opened off the visible window (UX.145).** The action row is a size container (`container-type: inline-size`, UX.27), which also makes it the containing block of the menus inside it. The menus' placement (Shoelace's Floating UI reads `offsetParent`) did not know that and measured against the page, so a menu opened as far above its button as the row is from the top of the page: off the top of a short window, and detached from its button on any window. The row is now `position: relative`, so the offset parent and the containing block are the same element. A browser test opens the Admin menu on every page that has it, at a laptop, a small laptop, a tablet, a phone and a short window, in both themes, and checks that it hangs at its button, lies inside the window and that every entry can be reached.

TODO items resolved: UX.145. Added: none.
