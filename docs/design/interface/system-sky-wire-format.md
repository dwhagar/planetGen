# System sky wire format: what the System View downloads for its background (VIEW.31)

Status: built (VIEW.31 and the glow image of VIEW.13, 2026-11-05). Plan: `reports/system-sky-from-the-system-2026-11-04.md` (section 7, label S13). Code: `src/planetgen/api/system_sky.py` (the build), `src/planetgen/web/lib/skypack.py` (the packed list), `src/planetgen/web/lib/skycache.py` (the disk cache), `src/planetgen/galaxy/sky_sprites.py` and `src/planetgen/db/sky_sprites.py` (the sprites). Readers: `src/html/static/skypack.js` (the stars, VIEW.33) and `src/html/static/skysprites.js` (the sprites, VIEW.36). Tests: `src/tests/test_system_sky_route.py`, `src/tests/test_skypack.py`, `src/tests/js/skypack.test.mjs`.

The System View's background is the sky its system sees: each star at its true direction and distance, placed where it was when its light left it. One request gives a page everything it draws behind the system, as one packed file.

## 1. The routes

| Route | Answers |
|---|---|
| `GET /system/<id>/sky?limit=` (the web layer; the page's script calls this one) | the packed list, `Content-Type: application/vnd.planetgen.sky`, `Cache-Control: no-store`, with `X-Sky-Source` (`cache` or `built`) and `X-Sky-Version` (the version key, section 5) |
| `GET /api/systems/<id>/sky?limit=` | JSON: `system`, `limit`, `epoch_unix`, `version_key`, `count` (stars), `sprites` (records), `report` (section 4), `seconds` (what the build took), `format` (`"PGK1"`) and `data` (the base64 of the packed list) |

`<id>` is the system's printed ID. `limit` is the limiting magnitude, 2 to 9, default 7.5 (the page lowers it on a phone; it is a user setting there, `skypack.DEFAULT_LIMIT` is only the route's default). A limit that is not a number or is out of range is a 400, an unknown system or one with no place in the galaxy a 404, an API failure a 502 (web route), all as `{"error": "..."}` JSON.

The sky is a snapshot at the system's own epoch (`star_systems.epoch_unix`; the orbit simulation's when that is empty): the page's orbit clock moves it on, so a sky stays right across orbital updates.

## 2. The packed list

One file, little-endian. The layout is the docstring of `skypack.py` and is held by `test_skypack.py` and `skypack.test.mjs` to the same bytes. In short: the magic `PGK1`, the star count, the sky's epoch, the limit, the velocity scale, the observer's velocity, the length of `meta`, `meta` (UTF-8 JSON), then the columns of the stars: offset x, y, z (float32 parsecs from the observer at the epoch), velocity x, y, z (int16, km/s relative to the observer), the absolute magnitude plus the dust of the line of sight (int16, steps of 1/512), the temperature (uint8, log over 100 K to 1,000,000 K, 0 for none), the kind (uint8: 0 charted, 1 bright star, 2 phenomenon, 3 cluster) and the index (uint32). A star is 26 bytes; 25,000 stars are 0.65 MB, about half that gzipped. A black hole, a neutron star or a quasar with a map luminosity is a star of kind 2, so it is drawn like one.

## 3. `meta`

`meta` is a JSON object, always present:

| Field | Meaning |
|---|---|
| `system` | the printed ID of the observer's system |
| `version_key` | the key the file is good for (section 5) |
| `limit` | the limiting magnitude it was built to (also in the header, as a float32) |
| `galaxy` | names the galaxy the stand-in velocities were drawn for |
| `nebulae` | the sprite records, below |
| `band` | the glow image, section 3.2: `{key, width, height, format, href}` |
| `charted` | the charted stars of the list by their index (the `index` column), section 3.3: `{"<index>": [printed star ID, name or null]}` (VIEW.37). Absent in a list built before it; the page then links nothing |
| `report` | the tier state and counts, section 4 |

### 3.1 The sprite records (`meta.nebulae`)

The nebulae, the supernova remnants and the dark and diffuse clouds that the system sees, widest first (at most 600, `system_sky.MAX_SPRITES`). Each record is the output of `physics/nebula_sprites.sprite_data` for one object, with its direction worked out from where the object was when its light left it (`light_travel.retarded_offsets`, the same law as the stars). This is the shape `skysprites.js` reads.

| Field | Type | Meaning |
|---|---|---|
| `id` | text | the object's printed ID (a nebula's or a remnant's) |
| `kind` | text | `"nebula"` or `"supernova_remnant"` |
| `family` | text | a key of `tuning.NEBULA_FAMILIES` (`emission`, `reflection`, `planetary`, `dark`, `diffuse`, ...), `"supernova-remnant"` or `"plerion"`; sets the colour |
| `direction` | `[x, y, z]` | the unit vector to its centre in the galaxy frame (z up), light-time retarded; the page turns it into its own axes as it does for the stars |
| `distanceLy` | number | the length of the retarded offset |
| `angularRadiusRad` | number | the radius of its disc, `asin(R / d)`; `pi` when the observer is inside |
| `inside` | boolean | the observer is within the ball (the page draws nothing; the haze is the dust of the stars) |
| `extinctionMag` | number | how much it dims what lies behind it along its middle, magnitudes; the page scales it by the chord |
| `surfaceBrightness` | number or `null` | its mean V surface brightness, magnitudes per square arcsecond, dimmed by the dust between; `null` when it has no light of its own (a dark cloud, a diffuse one, a plerion) |
| `visible` | boolean | the surface brightness is under 23, the naked-eye limit |
| `seen` | boolean | for a remnant, its explosion has been seen (it exploded longer ago than its light takes to arrive); always true in the list, because an unseen remnant is not sent |

A record is sent when it is wider than 0.002 rad (SPR1), its explosion has been seen, and it either shines under the limit (`visible`) or dims by at least 0.02 magnitude (SPR2). The counts of what was left out are in `report.sprites`.

### 3.3 The charted stars (`meta.charted`, VIEW.37)

A star of kind 0 (charted) has a page, and the page links it when the star is picked. `meta.charted` maps the star's `index` (as text) to `[printed star ID, name]`; the browser builds the link from the page's `data-sky-star-url` (`/object/star:{id}`). Up to 5,000 stars are named (the brightest first). A scatter, cluster or compact-object star is not here: until GEN.1254 and GEN.1255 give it an Object ID and a name it is told ("Star not yet charted", "Star of a globular cluster") but not linked. The sky's version (`SKY_VERSION` 3) changed with this, so a list stored before it is built again.

### 3.2 The glow image (`meta.band`, VIEW.13)

The light of the stars too faint to be listed, as one small image the band layer (`skyband.js`, VIEW.34) lays behind the stars. It is not in the packed list: `meta.band.href` (`/system/<id>/sky/band?limit=<limit>`) is its route, and the page asks for it with `&key=<meta.band.key>` so a stored copy is served by name.

| Route | Answers |
|---|---|
| `GET /system/<id>/sky/band?limit=&key=` (the web layer) | the bytes, `Content-Type: application/octet-stream`, `Cache-Control: no-store`, with `X-Sky-Source` (`cache` or `built`), `X-Sky-Band-Key`, `X-Sky-Band-Width`, `X-Sky-Band-Height` and `X-Sky-Band-Format` |
| `GET /api/systems/<id>/sky/band?limit=` | JSON: `key`, `limit`, `cell` (the 50 parsec cell of the observer), `width`, `height`, `format` (`"rgb16f"`), `seconds` and `data` (the base64 of the bytes) |

The image is 256 by 128 texels (`sky_band.EQUIRECT_WIDTH`, `EQUIRECT_HEIGHT`), equirectangular in the sky's galactic axes, as `skyband.js` `equirectUV` reads it: row 0 at the north pole, column 0 at longitude -180 (longitude is `atan2(y, x)`, latitude `asin(z)`), a texel's value taken at its centre. Each texel is red, green and blue as little-endian half floats (`rgb16f`, 196,608 bytes), linear light in the page's unit: `10^(-0.4 (SB - 22))` for a surface brightness `SB` in V magnitudes per square arcsecond (a texel of 1.0 is 22 mag/arcsec^2, 0 is no light), times the tint of a 4500 K blackbody scaled to a luminance of 1 (red 1.30, green 0.95, blue 0.65; stand-in GLW1). The surface brightness is `sky_band.band_for`'s: the integrated light of the stars fainter than the limit, marched through the galaxy's density and dust from the centre of the observer's 50 parsec cell, so every system in a cell, at one limit, shares one image. It is worked out in the API process the first time (tens of seconds on a galaxy with many clusters) and then kept in memory; the web layer keeps `b<key>.bin` beside the skies and `skycache.clear()` removes it.

`key` is 32 hex characters: a hash of the cell, the galaxy (its seed and the numbers of its shape), the limit, the band's own constants and the image's (size, reference, tint, `SKY_VERSION`). It does not follow the cluster field, which the band reads but the key leaves out; after clusters are regenerated, clear the cache.

## 4. The report

`meta.report` (and `report` in the API answer) says how much of the sky is built, for the page's "filling in" note:

| Field | Meaning |
|---|---|
| `limit` | the limit it was built to |
| `visible` | stars in the list |
| `stars_in` | stars read from the stored rows before the limit cut |
| `tiers` | per tier (`charted`, `bright_star`, `phenomenon`, `cluster`): `in`, `visible`, `faintest_magnitude`, `farthest_pc` |
| `excluded_own_system` | stars of the observer's own system, left out |
| `no_photometry`, `unknown_epoch`, `not_formed_yet`, `other_phase`, `unbuilt_velocity` | the counts `galaxy/sky.py` keeps: stars with no luminosity or temperature, no stored epoch (taken to hold at the sky's time), light that left before the star formed, stars seen in an earlier phase, and stars whose velocity is the stand-in |
| `built_to_pc` | how far out the stored stars are complete to the limit, parsecs (`sky.built_to_pc`) |
| `complete_above_sol` | the luminosity above which every cell of the sphere is complete, solar, or `null` when none was drawn |
| `wanted_pc` | how far the stored stars must be complete for this observer and limit: the reach of a star of the galaxy scatter's own floor (`generation/sky_tiers.py` `reach_pc`), past which the pre-placed scatter supplies every star that shows; `null` for a galaxy with no scatter level |
| `filling` | `built_to_pc < wanted_pc`: the tiers are not all stored yet, so the page shows its note and asks again later |
| `sprites` | `in`, `sent`, `too_small`, `too_faint`, `not_seen`, `inside` |

The cell levels of a sphere of more than 200,000 sectors are not listed (`CellLevels.truncated`, SKV3): such a shell counts as complete only where its luminosity floor is above the galaxy scatter's level. `filling` is only as exact as the two ledgers it compares: VIEW.25 (charting builds the tiers) records the levels the tier table asks for, and sharpens `built_to_pc` to the table's own floors.

## 5. The cache

The web layer keeps each packed list on disk beside the Galaxy Map's nodes (`web/lib/tilecache.py`, PERF.129): `<cache dir>/<db>/<generation>/k<system>_<limit in hundredths>_<version key>.bin`. The version key is 16 hex characters from the galaxy's tile stamp (the sectors, systems and bright stars stored, and the generator epoch) and `skypack.SKY_VERSION`, so a sky is rebuilt when stars are stored or the generator changes, and not on an orbital update. A new generation of the tile cache drops every sky with the nodes; the size budget (`cache.max_mb`) counts them and drops the least recently asked for; a sky built for a newer key removes the older copies of the same system and limit. `skycache.clear()` removes every stored sky: it is the clear the admin page of PERF.134 calls beside the nodes'. One request at a time builds a given sky (`singleflight`); the rest wait for its file. A cache that is off or unwritable only means each request builds it. The API response is never put in the shared page cache (`pagecache.UNCACHED_SUFFIXES`).

## 6. Stand-ins and what is not here

The sizes and cuts above are registered in GEN.725 (`docs/design/research/gap-filling-plan.md` section 6.179: SPR1, SPR2, SPR3, SKV1, SKV2, SKV3). The physics of a sprite is VIEW.28's (section 6.173), the stars' is VIEW.23's (6.175). Not here, by decision (Boss, 5:29 PM EST on 2026-11-04): reddening by dust (VIEW.11), the full dust of nebulae (VIEW.8), the HYG calibration (VIEW.12) and the full glow ray-march (VIEW.13's full version), which stay in 9.1. The glow image the band layer samples (VIEW.34) is built here (section 3.2): the light version of VIEW.13, `galaxy/sky_band.py`. The 9.1 twins (VIEW.7, VIEW.8, VIEW.11 to VIEW.13) read the same list and the same records.
