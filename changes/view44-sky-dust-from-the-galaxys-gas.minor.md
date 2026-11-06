**Headline:** The sky of a system is seen through the dust of its own galaxy's gas, so a small galaxy's sky is no longer smothered by the Milky Way's dust written in fixed parsecs.

### Changed
- **The sky route, the glow band and the sprites use the dust of the galaxy's own gas.** `api/system_sky.py` builds the stored galaxy's `DustProfile` (`dust_model.dust_for_shape`, VIEW.43) and gives it to the stars (`sky_from_database`, which takes it by default from the stored shape: `sky.default_dust`), the nebula sprites, the glow band and the "wanted" reach. Before, every galaxy had the Milky Way's law in fixed parsecs (9.6 magnitudes per kiloparsec in the plane of the 481 pc test galaxy at its site, 0.05 now). `av_per_kpc_at_sun = 0` still means a dust-free sky.
- **The glow band's cache key follows the dust tables** (the digest of the profile, `DustProfile.key()`), and `SKY_VERSION` is 4, so every stored sky and band is built again.
