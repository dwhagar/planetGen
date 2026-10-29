**Headline:** Superbubbles on the Galaxy Map now have a crisp rim and a fainter fill instead of three crossing rings.

### Changed
- **Clearer bounds and more transparency for superbubbles (MAP.236).** A bubble was three faint great circles (opacity 0.14), so its bound read as two crossing ellipses and the inside had no edge. It is now its sphere: a faint even fill (0.06 on the dark map, 0.07 on the light) and a crisp rim about two pixels wide where the sphere meets the line of sight, so the silhouette is the bound, in a lighter blue on the dark map and a darker one on the light. A camera inside a bubble, as the Sun is, sees the shell round it with no doubled layer. The look is in the new `static/superbubble.js`; no stored value changes.
