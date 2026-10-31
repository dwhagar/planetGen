# Rendering bodies from composition (design draft, 2026-10-24)

Asked by Boss, 10:15 PM EDT 2026-10-23: the 3D rendered objects get texture, colour and atmospheric glow from composition and the other computed factors. Planets with oceans get a very rough vector-style picture of ocean, cloud and land colours; gas giants show colour and layers like Jupiter and Saturn, roughly; all rendered planets and moons are replaced with newly rendered surfaces; stars get a full paint job from luminosity and spectral class; phenomena (comets, rogue planets, nebulae, remnants) use as much of the calculated data as possible.

Status: approved by Boss (10:28 PM EDT 2026-10-23, decisions in section 10). Built so far: R1 `bodylook.py` (MAP.271), R3 `bodytexture.js` (MAP.272), R4a `solidpaint.js` (MAP.273), R5 land colour from the stored ground types (MAP.274, section 4.1) and R7 `starpaint.js` (MAP.276) on the System Map, and R9a (MAP.278: rogue planet, comet, black hole, neutron star and asteroid field looks from stored data; MAP.286: remnant shells by age on the Sector Map, shell hollow, density, opacity and colour interpolated on a log age scale across the generator's 50 yr to 100 kyr span, a stand-in look at tweak scale L0, the shell texture and emission-class colour still open); the rest follows the plan in section 11. Items: MAP.270 to MAP.285, DB.96, DB.99, GEN.670, UX.158.

## 1. What renders today (read from the repo)

| Surface | File | What it does now | Data it reads |
|---|---|---|---|
| System map spheres (stars, planets, moons) | `src/html/static/systemmap.js` | One scissored shared canvas. Solid worlds are a flat `MeshStandardMaterial` tinted by `data-color`. Gas giants get `makeBandTexture` (a 4 x 256 sine-sum stripe, latitude only) and a ring texture. Glow is a back-face shell from `bodyRendering.js`, coloured by `glowColorForTemp` (surface temperature text parsed from a label). | `data-color`, `data-bodytype`, `data-hasatmosphere`, `data-surfacetemp` |
| Planet colour on the page | `src/planetgen/web/maps/systemmap.py` `_CLASS_COLORS` | A hand-picked colour per class letter. This is a fixed table, not computed. | the class letter only |
| Star spheres | `bodyRendering.js` `makeStarSurfaceTexture` | 128 px sine granulation, one base colour. | one colour |
| Sector scene | `src/html/static/sectorscene.js` | Stars are points of light; nebulae, asteroid fields, black holes, neutron stars are radial sprites or small spheres with fixed recipes. | kind, a few colours |
| Phenomenon page | `src/html/static/phenomenonrender.js` | Neutron star, black hole, quasar, rogue planet, interstellar comet. Rogue planet texture is a hard-coded blue-grey (`hsl(210, ...)`) with random craters or a grey-blue band; comet is a lumpy icosahedron. | `view` JSON (`gas`, `heat`, a few numbers) |

The key finding: the colour work is built in Python but not stored or sent.

- `physics/atmosphere_colour.py` (GEN.416 to GEN.418) returns an `AtmosphereColour` (sRGB, a colour word, ozone, haze, dust, cloud) from the air column, the star's spectrum and the ground albedo.
- `physics/giant_colour.py` (GEN.419) does the same for a giant at the cloud-top reference level.
- `models.py` has no colour, albedo, haze or cloud columns on `planets` or `moons`; the marker fill still uses `_CLASS_COLORS`. The limb glow is derived (MAP.244, built): `bodylook.air_colour` works the atmosphere colour of the stored mixture and ozone column under the body's own star (`_star_temperature_k`, put on each row by `systemmap.put_star_temperatures`) at render time, about 3 ms a body, and the look record's `glow` is its hex; a giant's is the colour of the air above 1 bar (no stored cloud deck yet). A row with no star temperature or no stored mixture has no `glow`. Storing the colour as a column (R2) is still open and would only cache this. MAP.277 (R8, built): the look record also carries `glowStrength`, `glowWidth` and `twilight` from `bodylook.air_glow`. The strength is `0.85 (1 - exp(-column / (0.6 Earth columns)))` of the air's column `P / g` (the Earth's air 0.69, Mars's 0.02, Venus's and Titan's the most; under 0.01 the air has no glow fields and the body is drawn bare), a giant's is a fixed 0.7, the width is the scale height `k T / (m g)` over the radius times 50 held to 0.03 to 0.2 radii (left out when the gases are not stored, and each view then takes its own default), and the twilight (0 to 0.2 of a cosine) moves the day and night line that far onto the night side and softens it by the same amount in the body view (`lit` in `bodyview.js`) and sets the night-side fill of the System Map and the 3D system view (0.15 and 0.04 for no air, the old 0.35 and 0.12 for the Earth's). The numbers are display choices (a picture, not physics), picked so thin air shows almost nothing and thick air most.
- Stored today and usable: `ocean_fraction`, `land_fraction`, `ocean_depth_km`, `ocean_class`, `ice_shell_km`, `hydrosphere` (dry, vapour, ice, ice-covered ocean, surface ocean, hycean), `water_mass_fraction`, the partial pressures `p_*_kpa`, `uv_surface_index`, surface temperature, radius, mass, rotation and obliquity where stored, magnetic data (`dipole_class`), and the planet class.
- Not yet stored: cloud cover, albedo (GEN.412), haze/dust/cloud depths (GEN.414), ground (surface) colour, ice fraction (the polar cap size), mineral surface type.

## 2. Principles

1. Everything is a pure function of stored values plus the body's Object ID (96-bit) as the seed. Same body, same picture, on every visit and every machine.
2. No fixed Solar System content: no Earth continent shapes, no Jupiter Great Red Spot image. Shapes come from seeded noise; colours come from the computed composition.
3. No raw datasets in the repo. The only tables the renderer needs are small derived ones (blackbody colour by temperature, a mineral colour table, a few dozen numbers), each labelled stand-in in the gap plan registry where not sourced.
4. Rough on purpose. Boss asked for "very rough vector-style": flat or two-tone fills, soft edges, no photoreal shading. That also keeps it cheap.
5. Where the data does not exist yet, fall back to today's look and say so in the info panel (a "colour estimate" note) rather than inventing.
6. A planet shows the colour of the thing it is, not of its class letter. `_CLASS_COLORS` is retired when the new path lands (no shim, per the standing rule).

## 3. One architecture for all bodies

**Server side (Python, new module `planetgen/web/maps/bodylook.py`)** turns a stored row into a small "look" record, a few dozen numbers, and puts it in a single `data-look` JSON attribute (replacing `data-color`, `data-bodytype`, `data-hasatmosphere`, `data-surfacetemp`, which were parsed from labels). Doing the physics server-side keeps the JS thin and testable, and the same record serves the map, the page, and later the view from a planet (9.1).

**Client side (`bodyRendering.js`, split into `bodytexture.js` for painters)** turns a look record plus a seed into one equirectangular texture (and a few uniforms). One painter per body family, all deterministic from the seed.

**Texture generation:** painted once per body when it first becomes visible, on an offscreen canvas, 128 x 64 for the far tier, 256 x 128 for the near tier (see section 9), then cached by Object ID in an LRU. Painting is 2D canvas, not WebGL noise, so it works in the same scissored-canvas setup used today, with a later option to move to a fragment shader (section 9).

**Shared shader pieces:** the existing `GLOW_*` shell stays. A new small planet material adds (a) a day/night terminator from the star direction, (b) a cloud layer as a second slightly larger sphere with its own alpha texture, (c) the atmosphere limb colour from the computed atmosphere colour instead of `glowColorForTemp`.

**Seeding:** `look.seed` is the low 32 bits of the Object ID hash; the existing `makeRandom` LCG in `phenomenonrender.js` is moved to the shared module. (Server-side the draw.py wrapper is the one random source for anything generation-time; this is display-only noise, so a JS LCG is fine, but it must not change between releases without a note, since a changed look is a visible change.)

## 4. Solid worlds (planets and moons, same classes)

### 4.1 Inputs

| Look element | Computed from |
|---|---|
| Ocean colour | `ocean_class` (ice-sealed, chloride brine, acid sulfate, soda, neutral), `ocean_depth_km`, atmosphere colour (the sky tints a shallow sea), `ocean_ph` as a minor hue shift |
| Land colour | built (MAP.274): the stored areal mix of the land ground types (ultramafic, anorthositic, basaltic, felsic, clay, oxidised, carbon, evaporite, sulfur) weighted by their lab reflectance in the 450, 550 and 650 nm bands (`surface_composition.band_albedo`), read as blue, green and red and converted to sRGB; metal and the intermediate mixes have no measured spectrum and are left out (`landKnown`); with no stored ground, no land or under half the land known, the temperature and class estimate stays. The materials table's `hematite_fine` is nearly flat across the three bands (0.092, 0.086, 0.085), so an oxidised world comes out dark grey, not red: the broad bands average over the iron absorption edge, and a red tone needs finer bands or a sourced hematite spectrum (open) |
| Life tint | only if `phi_bio` / life scores are Yellow or better. Green for chlorophyll-type life; where the life chemistry differs, the tint is the colour of that life's pigment chemistry (Boss 10-24 answer to card 5), mixed in on land at a faint strength. Pigment colour follows the host star's spectrum (stand-in, flagged) and the life type, so the life type must be read from the stored life data (check what the four life scores and the class lane store; a gap item if the pigment is not derivable) |
| Ice and snow | surface temperature by latitude (a simple energy-balance band from obliquity and mean temperature), `hydrosphere`; polar caps, or full ice for `ice` |
| Cloud cover and colour | cloud depth and cover (GEN.414), cloud material (water, sulfuric acid, CO2, dust) from the dominant condensable in the mixture, hence white, pale yellow, or grey-tan |
| Atmosphere limb glow | the computed atmosphere colour (GEN.418) scaled by column thickness: thin air (Mars, under ~10 kPa) gets almost no glow; no atmosphere gets none (a real change from today, where `hasAtmosphere` is a flag) |
| Night side | emitted heat for lava worlds and hot bodies (the `heat` cracks idea in `phenomenonrender.js`), a dim city layer is NOT proposed (no population data on this surface until 9.2) |
| Brightness | geometric albedo (GEN.412) sets overall lightness so a dark carbonaceous moon is dark |

### 4.2 Painter (rough vector look)

1. Seeded low-frequency noise (3 to 4 octaves of value noise, wrapped in longitude) gives an elevation field. The sea level is chosen so the area under it equals `ocean_fraction` (a histogram quantile, so the numbers on the page and the picture agree exactly). Land is the rest, matching `land_fraction`.
2. Posterise to 3 or 4 flat tones, with a thin lighter coast outline: ocean deep, ocean shallow (`ocean_depth_km` sets contrast), lowland, highland. Soft edges via a 1 to 2 px blur. This is the "very rough vector" look; no hillshade.
3. Land tones come from the land colour; polar caps and high-altitude snow overlay by latitude where cold enough.
4. Cloud layer: a second texture, seeded noise with latitude banding (Hadley-style bands: more cloud near the equator and 50 to 60 degrees, clear near 25 to 30 degrees), thresholded to cover fraction, white or tinted, alpha 0.7. It rotates slightly faster than the surface for a cheap sense of depth. A world with an opaque cloud deck (Venus-like) skips the surface entirely and paints the cloud deck as the visible surface.
5. A dry world has no sea level step: elevation tones only, with crater speckle scaled by surface age (when the age is known) and gravity.
6. A world in vapour or `hycean` state shows the cloud or haze deck as the surface, tinted by the atmosphere colour.

### 4.3 Moons

Same painter, same classes ("moons same classes as planets", Boss 10-23). Tidally locked bodies put the large-scale land/ocean split by the sub-stellar point (an "eyeball" gradient: ice-free patch facing the star on cold locked worlds) as a later refinement, not first pass. Small irregular moons and asteroids use a lumpy mesh (as `interstellarComet` does) instead of a sphere, scaled by the tri-axial shape if one is stored.

## 5. Gas and ice giants

Inputs: the giant colour (GEN.419) as the base tint at the reference level; the mixture (methane gives the blue-green of the ice giants, ammonia clouds the cream of the Jupiter class, sulfur and phosphorus the browns and yellows when present, silicate and iron cloud at the hot end); equilibrium temperature and the Sudarsky label (text only, never read for colour, as the code states); rotation period and radius (fast rotators have more, thinner bands; the band count scales with the rotation and the Rossby deformation scale, a one-line law with a seeded range, stand-in); internal heat; obliquity; ring data if present.

Painter: latitude bands only in the base texture (as now) but driven by data, not a fixed sine pair.

1. Band count from rotation: `n ~ clamp(round(a * (R / R_J)^b / P_hours), 4, 24)` with a seeded jitter (stand-in coefficients to be fitted in the research item, not taken from Jupiter's image).
2. Alternating belt/zone tones: two lightness steps around the base colour, with chemistry shifting hue (more saturation of the minor chromophore where the mixture has it).
3. Zonal wind shear is drawn as a sinusoidal longitude displacement of each band edge (small), so edges are wavy; a few seeded oval storms (1 to 4, more with fast rotation and high internal heat), tinted by the minor chromophore, are the only 2D features.
4. A cold giant with a clear upper atmosphere (Neptune-like, methane absorbing red) is nearly featureless: low band contrast, as the optics predict, a deliberately different look from a cream, ammonia-clouded one.
5. Hot giants: very low contrast, the glow replaced by a night-side thermal emission shell; the very hottest get a visible thermal glow tinted from the blackbody colour.
6. Rings: today's gradient texture is replaced by a radial density profile from the stored ring data (inner and outer radius, optical depth, ring colour by composition: water ice pale, silicate dark red-brown). The Cassini-style gaps come from resonances with the moons that are already computed (a gap at a moon's 2:1 resonance radius), so gaps are real consequences of the system, not a cosmetic sine.

## 6. Stars

Inputs: effective temperature and spectral class (colour), luminosity and radius (apparent size and corona), age and evolutionary stage (giant, white dwarf, main sequence), [Fe/H] (line blanketing nudges colour redder, a tiny effect, optional), activity/flare data (`flare_irradiation_index` is stored for the host), rotation, magnetic star flag.

Colour: a single blackbody-to-sRGB function (Planck spectrum integrated against the CIE matching functions; `atmosphere_colour.star_spectrum` already builds a stellar spectrum and `_xyz` an XYZ conversion, so reuse it so the star's rendered colour and the colour the sky function assumed are the same number). This replaces `SPECTRAL_CLASS_COLORS` words and the per-star hex.

Texture by class (rough, all two or three octaves):

- O/B/A: nearly featureless, very bright, tiny granulation contrast, limb darkening mild.
- F/G: fine granulation cells, a few dark spots scaled by activity; limb darkening coefficient from temperature.
- K/M: larger cooler cells, more and bigger spots (M dwarfs especially), contrast higher because cool; flare stars get an occasional bright patch.
- Giants and supergiants: few, huge convection cells (the cell size scales with radius over gravity, a ratio the code has), very soft limb.
- White dwarfs: a smooth hot disc, with a faint photospheric tint set by the atmosphere type (the stored pollution flag lets a polluted dwarf get a faint line-colour tint); no granulation.
- Brown dwarfs: banded like a gas giant plus a dull red or magenta glow from the brown-dwarf type (T dwarfs methane blue-grey), reusing the giant painter with an emissive material.
- Neutron stars and black holes keep their current special builders in `phenomenonrender.js`, with the sizes and magnetic geometry tied to stored values.

Glow: the corona shell is sized from luminosity (bright stars get a wider, brighter glow, capped), tinted by the star colour, with the existing back-face shader. Fix MAP.269 first (stars beyond ~3 sectors are 1.5 px with no halo and a white core): the far-star point sprite should take the blackbody colour and a luminosity-scaled halo sprite, which is the same data and the same function.

## 7. Phenomena

Use the `view` JSON the server builds (`phenomenonrender.py`) enriched from the stored phenomenon row:

| Phenomenon | Painted from |
|---|---|
| Rogue planet | the same planet/giant painters, lit only by a very dim ambient light plus internal heat; surface temperature is low, so no ocean (frozen or none) unless a heat source is stored; ice, rock or giant look from its composition, not the fixed blue-grey; the `heat` crack layer stays and scales by internal heat flux |
| Interstellar comet | nucleus colour from the comet composition (GEN.488 coma classes, `comet-composition-tuning.md`): dark organics and red for volatile-poor, grey-white for ice-rich; coma colour from the dominant emitting species (CN and C2 green, ion tail blue from CO+, dust tail yellow-white), brightness and tail length from the activity and distance to the nearest star (zero tail far from any star, as is true); lumpy shape from the stored axis ratio |
| Nebula | the sprite colours from the nebula class (`nebula-and-asteroid-field-classes.md`): emission red/pink for H II, blue for reflection, dark for absorption, green-teal for planetary nebula O III, with the star inside setting the bright core; turbulence texture layers from seed; size and density from the stored values |
| Supernova remnant | an expanding shell texture with the age setting shell radius and filament contrast (`remnant-age.md`), colour by the dominant emission class |
| Pulsar and neutron star | spin rate and field strength set the beam cone geometry and the glow tint |
| Black hole | accretion disc colour and brightness from the accretion rate if stored (else off), jets only if the stored flag says so |
| Asteroid field | grain colour from the belt composition (carbonaceous dark, silicate tan, metal grey), density from the stored count |
| Quasar | unchanged shape; colour and jet scale from stored luminosity |

## 8. Dependencies on unfinished work

| Part | Needs | State |
|---|---|---|
| Any colour from composition on the page | colour stored or derived and sent to the page (new `bodylook` module); GEN.418 and GEN.419 produce it but nothing wires it | GEN.418/419 functions exist; wiring is new work |
| Planet land colour from bulk composition | the composition spine (GEN.411, GEN.428, GEN.429) | open; the spine is the biggest blocker in the report `top-5-blockers-direct-2026-10-24.md` |
| Cloud cover and cloud material | GEN.414 stored values | open |
| Albedo to set brightness | GEN.412 | open |
| Ocean and land picture | `ocean_fraction`, `land_fraction`, `ocean_class` | stored (GEN.88), usable now |
| Ocean fraction accuracy | GEN.441, GEN.461, GEN.462 | open (Blocker 4 lane) |
| Gas giant bands by chemistry | GEN.419 done; band strengths for CH4/NH3/H2O are still coarse (GEN.422), so colour carries `colour_quality = partial` | partial |
| Heavy-element/chromophore colours (sulfur, phosphorus, tholins) | atmosphere remainder (lane L) and sky colour (lane M, 9.0) | open |
| Atmosphere glow thickness | atmospheric pressure and mixture | stored |
| Star spectral colour | `star_spectrum` and temperature | exists |
| Star surface activity | host activity (`host-activity-tuning.md`) | mostly stored |
| Comet composition | GEN.488 (sourcing) and comet admission | open |
| Ring profile and resonance gaps | ring data and moon periods | check which are stored; ring generation may be a gap |
| MAP.269 far-star glow | independent bug | open |
| Colours of the 9.0 headline (atmosphere colour as sky) | 9.0 | the planet-from-space look here is the "from outside" half of the same function |

What can start now with stored data alone: the ocean/land/ice painter, the new star colour, the glow tied to pressure, the comet and nebula tints from class, the data attribute plumbing. What waits: land composition colour (spine), cloud and albedo (GEN.412/414), chemistry bands for giants (GEN.422), comet species (GEN.488).

## 9. Performance on the zoom map

Budget facts: the system map uses one shared canvas with scissored spheres; the sector scene draws many bodies at once; the infinite-zoom camera (MAP.146, done) cuts detail by projected size.

- **Tiers by projected pixel diameter:** under ~3 px a point sprite coloured by the look (no texture); 3 to ~24 px a flat tinted disc with a baked limb glow (the current thin sprite path, but with the right colour and brightness); 24 to ~120 px the 128 x 64 texture; over ~120 px the 256 x 128 texture plus the cloud layer plus the terminator; over ~500 px (a planet filling the screen, relevant to the view from a planet in 9.1) a procedural fragment shader for detail instead of a bigger canvas texture.
- **Texture memory:** 256 x 128 x 4 bytes = 128 KB per body. An LRU of 128 textures is about 16 MB; far less than the per-frame budget concerns. The painter takes an estimated 2 to 6 ms for a 256 x 128 canvas, so it is queued: at most one paint per frame, nearest to the camera first, with the tinted-disc tier shown until the texture is ready.
- **No texture work per frame:** rotation is a uniform, the cloud layer rotation is a uniform, everything else is baked.
- **Draw calls:** a system with a dozen planets and moons is a dozen meshes, trivial. A sector scene with thousands of stars stays on the single `THREE.Points` plus a halo points layer; only the nearest few stars get a sphere.
- **Server cost:** the `bodylook` record is computed from stored columns in the page build; no new queries beyond columns the page already reads, once the colour columns exist. Storing the computed colour (rather than re-running the two-stream model per request) is recommended: three small columns per body (RGB as an integer, haze and cloud cover are already in GEN.414), with a unit entry in `db/column_units.py` in the same PR.
- **Reduced motion:** one still frame, as `phenomenonrender.js` does today.
- **Determinism test:** a Python-side test that the look record for a fixed row is stable, and a JS test (the repo's `fakedom.mjs` setup) that the painter returns the same pixel hash for a fixed seed and look. This guards against accidental changes to how a world looks.

## 10. Decisions (Boss, 10:26 PM EDT 2026-10-23)

1. Version: not 9.0 only. A full plan spread across versions, the first groundwork in 8.1 and the rest in later versions (section 11).
2. Painter: 2D canvas textures cached per body.
3. Colour: computed on the server, stored next to the body, sent as one look record.
4. Spine: the finished land-colour feature waits for the composition spine, but at least one task is done before it (R1, R3, R7, R9a, R0 and R10 below need no spine).
5. Life tint: yes on Yellow-or-better worlds, but not always green: the tint is the pigment colour of that life's chemistry.
6. Small bodies: lumpy meshes under about 400 km.
7. Night side: thermal glow on hot worlds only, no city lights until 9.2.

## 11. Full plan, items by version (TODO assigns numbers)

The version rule is the roadmap's: a blocker is done earlier in the same phase or an earlier version.

### 8.1 groundwork (nothing here needs the spine or a new stored value)
| Tag | Item | Needs |
|---|---|---|
| R0 | (bug) Far stars on the maps take the blackbody colour and a luminosity halo (extends MAP.269) | none |
| R1 | `bodylook` record: server builds one `data-look` JSON per body from stored columns, replacing `data-color`, `data-bodytype`, `data-hasatmosphere` and `data-surfacetemp`. In 8.1 it carries what is stored today (water fractions, ocean class, hydrosphere, pressure, temperature, class); fields for later data stay empty | none |
| R3 | Shared painter module: seeded noise, LRU texture cache, one-paint-per-frame queue, size tiers, determinism hash tests (R11 for these) | none |
| R7 | Star painter: blackbody colour shared with `star_spectrum`, class texture, spots by activity, luminosity glow, white and brown dwarf variants | R1, R3 |
| R9a | Phenomena painted from classes: nebula class tints, remnant shells by age, asteroid field grains, black hole and pulsar geometry from stored values; rogue planet and comet get the shared painters with today's data | R1, R3 |
| R4a | Ocean, land, ice painter on stored data (sea level from `ocean_fraction`, caps from temperature), land hues from temperature and class family, "estimate" line in the info panel. This is the pre-spine task for card 4; it ships dark behind the look record if Boss wants the planet change held until the spine | R1, R3 |
| R10 | Research: giant band, storm and Rossby laws, life pigment colour by life chemistry and host spectrum, ring colour by composition (stand-ins registered in the gap plan) | none |

### 8.4 (new worlds and composition)
| Tag | Item | Needs |
|---|---|---|
| R5 | Land colour from bulk composition (oxidation, rock type, mineral colour table) replaces the R4a estimate, finishing the planet look | spine GEN.411, GEN.428, GEN.429; R4a |
| R2a | Store the ground colour and ice fraction on planets and moons (migration, epoch, `db/column_units.py`) | R5 |
| R6a | Giant painter: bands, storms and thermal glow from rotation, temperature and the GEN.419 colour; rings with resonance gaps from the moons | R1, R3, R10, GEN.419 |
| R13 | Small bodies: lumpy meshes for moons and asteroids under about 400 km, from stored shape | R3 |

### 9.0 (composition and atmosphere colour)
| Tag | Item | Needs |
|---|---|---|
| R2 | Store the atmosphere colour (RGB integer) and cloud cover from GEN.418 and GEN.419; MAP.244 follows | GEN.414, GEN.418, GEN.419 |
| R4b | Cloud layer, terminator and cover on solid worlds; cloud-deck worlds paint the deck as the surface | R2, GEN.414 |
| R8 | Atmosphere limb glow from pressure and composition (none for no air, thin for thin air) | R2 |
| R6b | Giant chemistry bands once band strengths exist | GEN.422 |
| R14 | Life tint: the pigment colour of the life chemistry, faint, on Yellow-or-better worlds | R10, R4b, life data |
| R9b | Comet colour from coma species | GEN.488 |
| R12 | Info panel line: how the colour was worked out (inputs, partial or full) | R1 |
| R15 | Night-side thermal glow on hot worlds | R4b |

### 9.1 (view from a planet)
| Tag | Item | Needs |
|---|---|---|
| R16 | Close-up tier: fragment-shader detail over about 500 px, shared with the sky view | R3, 9.1 sky data |

### Not built (Boss decision)
City lights wait for 9.2.

## 12. First 8.1 groundwork to start

R0 and R3 start first and in parallel (no blockers). R1 starts with them (it only needs the data attribute change). R7 and R9a follow as soon as R1 and R3 merge; R4a and R10 can start in parallel from the first day. R11 tests ride with each item.


## 13. First on-page use: the object page (UX.171)

The object page of a star, planet or moon (see [object-pages.md](object-pages.md)) paints a body with the System Map's own
painters, so what the page shows and what the map shows agree: `web/lib/bodyview.py` makes the look record from the body's
row, `static/bodyview.js` paints it (`solidpaint.js`, `starpaint.js`) and shades it as a sphere lit from the upper left on a
2D canvas (no WebGL needed), with the atmosphere's limb glow or a star's corona around it. Giants, small bodies and phenomena
arrive on the page as their painters merge. The caption is the look note; the colour account (UX.158) replaces it.
