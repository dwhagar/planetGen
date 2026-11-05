# The sky behind a system: what the System View draws as its background

How the System View's background is built: the stars, the glow, the nebulae and remnants that the system itself would see, each star at its true direction and distance and placed where it was when its light left it. The note holds the numbers the build uses, where each part lives in the code, and what is a stand-in.

Informs: VIEW.20 to VIEW.38 (goal #4907); VIEW.39 and VIEW.40 are optional and open

Status: built for stored stars (VIEW.20 to VIEW.37, 2026-11-05). The sky tiers that charting builds (VIEW.25, VIEW.26) and the calibration against a real catalogue (VIEW.30) are open; until they land the sky holds the stars that are stored and says how far it is built.

Plan and research: `reports/system-sky-from-the-system-2026-11-04.md`. The wire format is [system-sky-wire-format.md](../interface/system-sky-wire-format.md); the longer, full-chart design for the view from a planet (9.1) is [sky-view.md](sky-view.md), and its items VIEW.7, VIEW.8, VIEW.11, VIEW.12 and VIEW.13 stay in 9.1. Sources are cited as `[R0123]`; the list is on the References pages.

## 1. What Boss decided

Boss, 5:17 PM EST on 2026-11-04: the System View's background is the sky that system sees, each star at its true direction and distance and placed where it was when the light left it, not where it is at the galaxy's present time. Boss, 5:29 PM EST the same day (answers 1.2, 2.1, 3.3):

1. A lighter 8.2 backdrop from stored stars and the glow. The items VIEW.7, VIEW.8, VIEW.11, VIEW.12 and VIEW.13 are built in 8.2 as foundations of the view from a planet (9.0), and this part is built so they can reuse it.
2. The limiting magnitude is a setting. The default is 7.5 and a phone gets a lower one (6.5).
3. Missing stars are built when a system is charted, not by a visit. Reading the sky never starts a job.

Defaults taken without asking: aberration of the observer's motion is off (VIEW.39 would add it); stars fainter than the limit are the glow with a seeded grain, never invented one by one; light travel time is applied to every star; no database migration.

## 2. Where each star appears

Light leaves a star at time `t - tau` and reaches the observer at `t`, so the star is seen where it was `tau` earlier, with `tau = |apparent offset| / c`. Stars move in straight lines in the stored model (position, velocity, epoch). Per star, with `d0` the offset from the observer at the star's stored epoch, `v_rel` its velocity minus the observer's and `t` the years since that epoch:

    apparent offset = d0 + v_rel t - v_star tau        tau = |apparent offset| / c

It is solved by iteration: two passes, three for a star faster than 1,000 km/s (`skysphere.FAST_KMS`). Each pass shrinks the error by `v/c`, about 1e-3 at 300 km/s and 3e-3 at 1,000 km/s. This is the standard light-time correction of astrometry [R5112]. The same arithmetic runs three times: in numpy for the catalogue (`physics/light_travel.py`, `galaxy/sky.py`), in the vertex shader that moves the stars every frame (`skysphere.js`) and in the page's picking (`skypick.js`), and tests hold the three to one another.

How much it moves a star, whatever its distance: about `v/c` radians. That is 20 arcseconds at 30 km/s, 2.5 arcminutes at 220 km/s and 11.5 arcminutes at 1,000 km/s, which on a 60 degree view across 1,000 pixels is 0.1, 0.7 and 3 pixels. The path of a star bends away from a straight line by `(v T)^2 / (2 R)` over a light time `T`; for a star at 10 kpc that is 0.003 pc, about 0.01 arcsecond, so no orbit is integrated for the sky.

Three more effects use the same `tau`:

* **Age.** A star seen at distance `d` is seen at its age minus `tau`. If that is negative it had not formed, and it is left out (counted in the report). A star whose phase ended less than `tau` ago is seen in the phase before; the sky keeps its stored luminosity and temperature and counts it (`other_phase`). The stellar track at the earlier age is VIEW.12's, in 9.1.
* **Remnants.** A supernova remnant younger than the light time has not been seen. `physics/light_travel.py` carries the rule; the sprite builder marks such a remnant unseen (`remnant_seen`).
* **The galaxy's present.** Every stored position carries its own epoch. The orbit clock (up to 100 years a second) moves the stars, so nearby ones drift across the view and one at 100 pc does not move, which is the true picture. The list stays right across the monthly orbital update, so the cache is not cleared when objects move.

Not done: aberration (172 arcseconds at 250 km/s, under a pixel) and gravitational lensing. The galactic nucleus is in the glow; a dark black hole is not drawn.

## 3. Which stars

The rule of sky-view.md stays: resolved stars come from the real generator and everything fainter is light, and no resolved star is invented for space a later fill could not reproduce.

1. **Stored stars (exact).** Every star of a charted system, the pre-placed `bright_stars`, the scattered phenomena (hypervelocity stars with their stored velocity, neutron stars and black holes at their map luminosity, quasars) and cluster members. A close pair fuses to one point. The observer's own system is left out; its stars are bodies.
2. **Sky tiers.** The stars between the scatter floor and the limit, drawn by the existing deterministic backfill and stored when a system is charted (VIEW.24 is the table, VIEW.25 the build, VIEW.26 the command for sectors charted earlier). Until a system has them, the sky shows what is stored and a status line says "built to X pc, filling in"; the page asks again every 30 seconds, at most 20 times.
3. **The glow band** of everything fainter than the limit (section 5).
4. **Extended objects** as sprites (section 6).

How many: the limit is a magnitude, not a distance. At 7.5 about 25,000 stars are resolved in a real sky, at 6.5 about 9,000 and at 8.5 about 75,000; each magnitude of limit multiplies the count by 2.5 to 3. At 7.5 the density in the galactic plane is a few stars per square degree and a pixel of the view is 0.06 degree, so no two stars blend. The page draws at most 60,000 stars and a phone 15,000 (`skysphere.MAX_STARS`, `MAX_STARS_PHONE`), keeping the brightest. A star fades in over the last 1.5 magnitudes before the limit so none pops in when a tier finishes.

## 4. Brightness and colour

* **Magnitude.** `m = M + 5 log10(d / 10 pc) + A_V`, with `d` the retarded distance (`physics/sky_photometry.py`). The bolometric correction is the one of `physics/photometry.py`; stars under 4,100 K take the main-sequence table of Pecaut and Mamajek [R3469].
* **Dust.** The smooth galactic law: `A_V` of 1.0 mag per kpc in the plane at the solar circle, scale length 2.26 kpc and height 134 pc [R1158]. The server adds the dust of a star's line of sight into the list's magnitude (the packed column), so the browser does no integral. Reddening by dust is not applied (VIEW.11, 9.1); the law for it would be [R0651].
* **Look of a point.** The Galaxy Map's laws (`starmagnitude.js`, `starlight.js`) with `skylight.js` on top: the middle of a star is at full strength at magnitude 4.5 and brighter, its Gaussian width is 0.7 pixel at 4.5 plus 0.45 pixel for every magnitude brighter, a point is drawn 4 sigmas wide and never wider than 12 pixels. The colour is the blackbody of the temperature, blended toward white as the star dims (all of the colour is kept at magnitude 1.5 and brighter, 35% of it at 5.5 and fainter).

## 5. The glow band

The light of the stars too faint to draw: for a direction `u` from the observer `P`, `J(u) = integral of j0 density(P + s u) f_faint(s) 10^(-0.4 A_V(s)) ds`, in the galaxy's own stellar density (so it follows the model's disc, bulge and arms), 400 steps spaced as the square of a uniform fraction out to 30 kpc. `f_faint(s)` is the share of a population's V light that comes from stars too faint to reach the limit at distance `s`, so the band adds to the drawn stars without counting them twice. A lower limit (a phone) draws fewer stars and leaves more light in the glow.

Surface brightness is `SB = -0.19 - 2.5 log10(Omega J)` in mag per square arcsecond; a Milky-Way-like galaxy at the Sun runs 21.6 to 23.0, so the band is values of about 0.4 to 1.5 in the image's unit (1 is 22 mag per square arcsecond). The server keeps one image per 50 pc cell of the observer's position (moving 20 ly changes it by at most 0.004 mag), as a cube map of 6 faces of 32 pixels, and sends a 256 by 128 equirectangular image of three half floats a texel (`rgb16f`). The page samples it by direction in a shader, draws it as `0.4 (1 - exp(-0.5 v))`, multiplies by a seeded grain of plus or minus 20% over 0.12 degree cells so it does not look like smooth paint, and dithers the 8-bit output by one level. Nebula dust and the full ray-march with the nebulae are VIEW.13 (9.1).

## 6. Nebulae, remnants and compact objects

Each stored nebula and supernova remnant becomes a sprite record: direction, angular radius, surface brightness, extinction (`galaxy/sky_sprites.py`, at most 600, the widest first, none narrower than 0.002 rad or 0.115 degree). The page draws a glow sprite above the limiting surface brightness of 23 mag per square arcsecond (full strength 5 magnitudes brighter) and a dark sprite that dims what is behind it for a cloud with no light of its own. A luminous neutron star, black hole (disc) or quasar is a point at its map luminosity and temperature. Which stars a sprite hides is decided by the depth test, because the sky is far past the bodies: a planet in front of a star hides it.

## 7. How it fits the page

* Only the standalone System View (`systempage3d.js`). In the Galaxy Map the real stars are the sky; VIEW.40 (optional, after the continuous zoom) would cross-fade the map's stars into the sky as the camera nears a system.
* One points object, one band sphere and the sprites, all centered on the camera and drawn first (a sky at infinity). Camera moves of tens of AU are ignored: 40 AU shifts a star at 1 pc by 40 arcseconds.
* A system's orbital plane is not tied to the galactic frame in the stored model, so the galactic plane is laid in the scene's reference plane (the galaxy's north pole is the scene's up). A stand-in (SKD5), replaced when a system gets an orientation.
* The preference (on unless the browser remembers off) is on the page; a failed request, or a system with no place in the galaxy, leaves the plain backdrop and says why.
* **Picking (VIEW.37).** Hovering a sky star shows what it is (the charted star's name or the kind), its class, distance and magnitude and "seen as it was N years ago"; a click opens a card, and a charted star links to its page. The pick uses the same light-time arithmetic as the card, within 10 pixels of the pointer. Scatter and tier stars have no ID until GEN.1254 and a name until GEN.1255, so they are told but not linked.
* Cost: a sky of 25,000 stars is about 0.65 MB packed (about 0.4 MB gzipped), one draw call and no per-frame upload. The first visit costs a query and a pack; the rest is a cache hit (`Cache-Control: no-store` on the route, the cache is the server's, cleared from the same admin page as the map tiles). Gaia Sky streams a catalogue by magnitude in the same way the tiers do [R5113].

## 8. What is exact and what is a stand-in

Exact: the direction and distance of every stored star (to the 1 mpc rounding of the stored position); the light-time position of a star given its stored velocity; the magnitude given luminosity, temperature, bolometric correction and distance. Everything else below is a stand-in in the GEN.725 registry (`docs/design/research/stand-in-values.md`), with its tier and tweak level:

| Piece | Rows | Section of the plan |
|---|---|---|
| Photometry: dust law, cool-star table | SKP1 to SKP4 | 6.171 |
| The sky catalogue: velocity of an unbuilt star, the star with no epoch | SKC1, SKC2 | 6.175 |
| The tier table | SKT1 to SKT6 | 6.176 |
| The glow band: luminosity density, faint-light shares, tint | SKB1 to SKB3, GLW1 | 6.174 |
| Sprites | SPR1 to SPR3, SKV2 | 6.179 (route) |
| The route's limits | SKV1 to SKV3 | 6.179 (route) |
| Drawing of band, sprites and plane | SKD1 to SKD5 | 6.180 |
| The look of a star, the star caps, the page's reach | SKL1 to SKL8 | 6.181 |

The velocity of an unbuilt scatter or tier star is the circular speed at its place plus a dispersion of (35, 25, 20) km/s drawn from the star's own key, so it is a pure function; a wrong 30 km/s is 20 arcseconds, invisible. GEN.339's seeded dispersion should key on the scatter seed so the fill agrees later.

## 9. Left for later

* VIEW.25 and VIEW.26 (charting builds the tiers, and the one-off command), VIEW.30 (calibration against a real catalogue on a developer machine: about 9,000 stars at 6.5).
* VIEW.39 (aberration, off by default) and VIEW.40 (cross-fade).
* GEN.1254 and GEN.1255 (an ID and a name for scatter stars), then the links.
* A system orientation to the galaxy (SKD5).
* In 8.2 (foundations of the 9.0 view): reddening and desaturation (VIEW.11), the nebula dust (VIEW.8), the stellar track at the earlier age (VIEW.12), the full ray-march (VIEW.13), the view from a planet's surface.

## Sources

[R5112] aberration and light-time corrections; [R5113] streaming a star catalogue by magnitude; [R1158] the dust disc; [R0651] the extinction law that VIEW.11 will use; [R3469] cool main-sequence stars.
