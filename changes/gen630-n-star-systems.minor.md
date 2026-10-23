**Headline:** A star system can now hold more than two stars, stored and loaded with each star's place in its tree of orbits.

### Added
- **The stars of a system beyond the pair (GEN.630).** `StarSystem.member_stars` holds each star after the first two and `star_slots` says which child of which orbit node every star is. A system of three stars is saved with a `member` row of `stars` for the third (schema v133: `stars.role` takes `member`, and the new `stars.slot` is 0 or 1; a stored pair gets its slots from its roles) and is loaded back with it, and it goes through the JSON form (serialization version 6). The fill still builds single stars and pairs, so no seed changes; the draw of stars with a third member is GEN.631.

### Changed
- **The readers that took a system to have one or two stars handle more.** The sector list orders a system's stars primary, secondary, then members; a member is placed in the galaxy and in the system scene at its offset in the tree, the scene's tree reads each star's stored node and slot, the system map draws the pair and the stars table lists every star, a facility around a close pair counts only the pair, and the stars of a system number their body slots by role and place among the members.
