**Headline:** The Galaxy Map's cards fit the map's window and a system's card lists its nearest systems.

### Changed
- **The card follows the map's window (MAP.295).** The popup a picked object opens inside the map was a fixed size, so on a smaller map window part of it could run off the edge. It now takes a share of the map's width (more on a wide map, less on a narrow one), is never taller than the window, and its text is a little smaller than before and shrinks a step at a time until everything it says fits. What cannot fit even at the smallest text (a very short window) scrolls inside the card: the wheel over it scrolls the card while it has more to show, then goes on to zoom the map. The card also sits under the map's Menu panel so it never covers the Menu's buttons.
- **A system's card names its nearest systems (MAP.296).** Under its position a star system's card (the habitable-world markers' card) now lists the systems nearest to it, nearest first with their distances, each name a link to that system's page. They come from the same stored neighbours the system page's "Nearest:" line names.

TODO items resolved: MAP.295, MAP.296. Added: none.
