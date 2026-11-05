**Headline:** Opening a saved course now says when a jump through unknown space has since been filled in with known sectors, or the reverse.

### Changed
- **A saved course keeps which of its hops were jumps through unknown space (NAV.39).** The record the Save button keeps now carries
  an `unknown_space` flag on each hop of the route. When the course is opened again and the route runs through the same systems,
  the list of what has changed adds a line for each hop whose flag differs ("Hop 1 of 2 is now through known space; it was a jump
  through unknown space."). A course saved before the flag was kept says nothing about it, so no old record is flagged wrongly.
