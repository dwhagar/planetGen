**Headline:** A planet's ring is now as icy as its age and the dust that fell on it leave it: Saturn-like rings stay bright, thin old rings go dark.

### Changed
- **The stored ring ice share is the share the ring has now (GEN.1055, R6).** `rings.apply` runs `rings.darken` after the ring is built: the dark dust that falls on the ring over its age (`physics/ring_darkening.py`, scaled from Saturn's infall by distance from the star and the planet's pull) lowers the ice share the ring began with, much further for a thin ring than for a massive one. A ring with no age (a resonance ring of a small body) and a ring with no ice are unchanged. No draw is added or moved; the values change.
- The ring shadow (R8, `ring_darkening.shadow_transmission`) stays unread until the latitude-band climate (GEN.1108) is stored.

TODO items resolved: GEN.1055 part R6 (R8 stays open until the climate reads the shadow).
