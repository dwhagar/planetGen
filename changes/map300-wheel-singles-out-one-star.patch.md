**Headline:** The wheel now zooms in far enough on a dense cluster to single out one star from its neighbours.

### Fixed
- **MAP.300: the mouse wheel stalled short of a star in a dense cluster.** With no sector open, the wheel only knew the stars of an open sector, so a pointer on a drawn star zoomed toward the focus plane or the column behind it and the eye stopped well short of the star. The wheel now takes any star the map draws from its nodes' lists as the thing under the pointer, so the camera comes in on it to the least distance the map allows (a thousandth of a parsec) and two stars under a tenth of a parsec apart are each the star under the pointer on their own spots. No stored value changes.
