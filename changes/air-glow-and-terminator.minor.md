**Headline:** A world's atmosphere now shows in its picture: no air has no glow and a hard day and night line, thin air almost none, thick air a wide glow and a soft dusk.

### Changed
- **The glow of an atmosphere follows its air (MAP.277).** The look record carries `glowStrength`, `glowWidth` and `twilight` worked from the stored pressure, gravity and gases (`bodylook.air_glow`). The body view draws the halo at that strength and width and softens the terminator by the twilight; the System Map and the 3D system view draw the glow shell at that strength and width and light the night side by the twilight (a world with no air has a dark night side). Before, every body with air had the same glow and every body the same terminator. A rogue planet has none of these (no star lights its air).
