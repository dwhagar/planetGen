**Headline:** Stars, glows and charted regions seen through a nebula on the Galaxy Map take its hue.

### Changed
- **MAP.215.** The nebulae's raw sum (the volume layer's target, before the cap) is read by the star shader, the unresolved star glow and the charted-region layer: what is behind a nebula leans toward the sum's hue (`nebulavolume.tinted`, the same formula in the shaders as `nebulaTinted`), more the denser the core and levelling off at 0.6 so the tint stays subtle, and loses a fifth of its light at full tint so a dense cloud dims what is behind it. Nothing is tinted while no nebula is drawn, in the Sector Map, or with the nebula layers off. The density views are unchanged.

generation-output: unchanged

TODO items resolved: MAP.215.
