**Headline:** One wheel scroll now runs from the whole galaxy down to a sector without stopping, and the map no longer swoops to a preset view or throws away your turn.

### Changed
- **The wheel carries on across stages (MAP.263).** Scrolling in past the end of a window's zoom takes the block under the pointer (the window one size down, or the sector) and scrolling out takes the parent, with no pause between one stage and the next, no gesture cooldown, and no flight to a preset camera: the camera stays where it is, with the turn you gave it and its distance, and the next window's blocks fade in over it. The next window's data is asked for before the zoom reaches the end of its range, so the crossing tick finds it. The map's middle can be moved anywhere: the limit of 1.5 fits from the stage's own middle is gone. Opening a sector or a system with the wheel keeps the camera too. Clicking a block, Up, Home, Back and the breadcrumb still fly to their view.

TODO items resolved: MAP.263. Added: none.
