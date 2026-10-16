**Headline:** Gamma-ray bursts and the galactic nucleus now count toward a planet's ozone-loss flag.

### Changed
- **The ozone-loss flag only knew about nearby supernovae.** The lethal-event rate on a star now adds gamma-ray bursts (2 per Gyr at the Sun's distance from the centre, falling off outward; Melott 2004, Thomas 2005) and nuclear-activity episodes (10 per Gyr inside 2.9 kpc of the centre). The nuclear rate is a documented design estimate, not a measurement; the sources and the arithmetic are in the design note.
- `activity.supernova_rate_per_gyr`, `grb_rate_per_gyr` and `agn_rate_per_gyr` give the three parts; `lethal_event_rate_per_gyr` is their sum.

TODO items resolved: GEN.189 / Added: none
