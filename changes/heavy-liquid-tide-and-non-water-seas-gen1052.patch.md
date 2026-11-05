**Headline:** Seas of methane, ethane and other liquids now have an equilibrium tide, and a new module holds the exact tide numbers of a deep global liquid layer.

### Changed
- **GEN.1052.** `ocean_tides` took every sea for water. A body whose stored `surface_liquid` is another non-lava liquid (a methane, nitrogen, ammonia or carbon dioxide sea) and covers part of the surface now gets its spring and neap range and its `ocean_liquid` (the liquid named, its density at the body's temperature in the loading of the floor). In 60 systems 15 bodies gain the columns, all small (a locked moon feels the star's tide only). A liquid too heavy for a thin layer still gives no tide. New `physics/heavy_liquid_tide.py`: the Love numbers `k2` and `h2` of an incompressible body with a core under a global liquid layer, solved exactly with the core moving or fixed, checked against a homogeneous fluid (3/2, 5/2), a fluid surface (`1 + k2 - h2 = 0`), the fluid Earth (0.93 to 1.04) and the thin-sea loading factor, which it shows is exact for a floor that does not move at any depth. Mercury and the other liquid metals were already rows of the liquid table; no generated body has a molten core under a global liquid, so the module is registered `by design`.

generation-output: changed

TODO items advanced: GEN.1052.
