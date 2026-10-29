**Headline:** "Near this system" now measures from where the system really is, and four more tests that were red on main are green.

### Fixed

- **A system's page and the near command put the system in the wrong place (TEST.453).** The galaxy position an object page and `near` read for a system added the system's offset to its sector centre without turning it by the sector's orientation, so it sat up to a couple of parsecs from the point stored with it (the position every map and search uses). `near` and `nearest` from a system now measure from the stored point, and the system's own sector-local position on its page is unchanged.
- **A moon whose tide works out to nothing stored 0 and not none (TEST.453).** An icy moon on a circular orbit saved a tidal flux of exactly 0.0 while a rocky one saved none; both are none now, so a saved moon's tide is empty or positive.
- **A world around a black hole could stop a system from generating (TEST.453).** A black hole's temperature gives no visible light, and the atmosphere colour of a giant lit by it raised "the spectrum has no light between 380 and 780 nm"; a body lit by a star with no visible light now has no atmosphere colour, as one with no star temperature already did.
- **Two planet page tests expected a class name the page no longer has (TEST.453).** The summary paragraph is `summary inline-links` since UX.220; the tests now look for that.
