**Headline:** A planet's surface is now seen through its air: the ground takes on the sky's color, hazy worlds look softer, and cloud cover shows.

### Added
- **Surface tint, softening, and clouds from the atmosphere (MAP.350).** The look record of a solid world with air carries `veil` (from the air's optical depth at 550 nm), `haze` (from the stored haze and dust depths), and `cloudCover` and `cloudTau` (from the stored cloud). The solid-world painter tints the ground toward the sky color at the same brightness, lowers its contrast, blurs it a second time under a thick haze, and lays a seeded whitish cloud layer over it. A world with thin air, or with none of those values stored, is drawn exactly as before.
