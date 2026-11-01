**Headline:** The update's cache warm-up no longer fails on the Galaxy Map page while the site is held in maintenance mode.

### Changed
- **Warm-up stage 1 failed with "/galaxy answered 302".** The update holds the site closed while it warms the caches, and `/galaxy` sent the warm-up's own visit to the maintenance page, which counted as a failed page. The warm-up's visits are now let through, so stage 1 gets through on an empty galaxy and with maintenance mode on. Visitors are still sent to the maintenance page.
