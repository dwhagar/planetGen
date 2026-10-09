# Other galaxies: neighbours in the sky and a second galaxy

Plan for GEN.9. Stage 1 puts other galaxies into the generated galaxy's universe as objects seen from it, so a planet's sky (see [sky-view.md](sky-view.md)) or a map can show Andromeda, the Magellanic Clouds and the rest of the Local Group. Stage 2 adds a second, fully generated galaxy with its own skeleton, sectors and systems. The document gives the Local Group table, a `neighbor_galaxies` schema, the frames between galaxies, the URL and page changes, and how an existing single-galaxy database migrates. It is a plan only; no code exists.

Informs: GEN.9, VIEW.2, VIEW.3 (neighbour galaxies in the sky)

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result (mostly Wikipedia pages, URLs in Sources), [C] computed in this research, [R] recalled and not confirmed. The research environment could read search-result text but not the papers or VizieR, so no catalogue value here has been checked against McConnachie 2012 or NED; the evidence notes list what to verify.

## Summary

- **Defaults for GEN.9:** (1) neighbouring galaxies are the real Local Group, a hand-built table of about 25 rows shipped as a data file and loaded into a `neighbor_galaxies` table by `planetgen plan`, with an option to draw fictional ones from the seed; (2) no `galaxy_id` on any row, because (3) a second galaxy is a second database, which is already how the project works (`?db=`, `list_databases`, `galaxy_naming` keyed per database, Alembic per database, `galaxy_shape` a singleton). The only schema change for stage 1 or stage 2 is one new table.
- **Do not hand-key magnitudes from memory.** Recalled apparent magnitudes were wrong by 0.7 mag (IC 1613: 9.2 recalled, 9.9 found [S]) and 1.9 mag (Fornax: 7.4 recalled, 9.3 found [S]). Store absolute magnitude and distance, compute the apparent value, and check every row against McConnachie 2012.
- M31, M33 and farther objects look the same from anywhere in the Milky Way to within 5% in size and 0.1 mag. The Magellanic Clouds and the Milky Way's other satellites change a great deal and must be 3D objects whose direction, size and brightness are recomputed per observer [C].
- **The generated galaxy rotates the opposite way to the real one** (counterclockwise from +Z against clockwise from the real north pole, [galactic-potential.md](galactic-potential.md) section 5). Placing real neighbours by their galactic (l, b) therefore needs a defined embedding; section 2.2 gives the proper one (b flipped, longitude 90 along the rotation) and shows the plain one (b kept) is a mirror image.

## Decisions already taken

- Boss, 2026-10-01 (GEN.9): "Lay the groundwork for different galaxies within the same DB, this is a planning item. Develop a plan to have multiple galaxies first as objects in the neighborhood (i.e. we can see Andromeda from Earth kind of thing) but also I might need to have another galaxy at some point. So just plan and save as a planning document."
- One galaxy per database is the current state ([galaxy-coordinate-system.md](galaxy-coordinate-system.md), "Open questions still standing"); no decision to change it has been taken.
- Navigation frames are Boss's design of 2026-09-30 ([navigation-frames.md](navigation-frames.md)): "North" points toward the local dominant centre of mass. Section 6 extends that rule one level up.

## 1. Two stages

| Stage | What | Needs |
|---|---|---|
| 1 | Other galaxies as objects: direction, distance, size, brightness and type, so a sky or a map can show them | the `neighbor_galaxies` table, a data file, a step in `planetgen plan`, a drawing rule in the sky and on the Galaxy Map |
| 2 | A second fully generated galaxy with its own skeleton, sectors and systems | a second database, a `linked_database` value on the neighbour's row, the Intergalactic Frame, a galaxy list page |

Today the database holds one galaxy: one skeleton (`galaxy_shape`, `galaxy_layer`, `galaxy_column`), one sector grid, and coordinates in that galaxy's own frame (origin at its core, +Z north, +X the zero meridian).

## 2. Local Group facts

### 2.1 Table

Confirmed [S] and derived [C] values:

| Galaxy | Distance | V mag | Apparent size | M_V (derived) | Mean SB (derived) |
|---|---|---:|---|---:|---:|
| LMC | 49.6 kpc (49.59 +/- 0.09 stat +/- 0.54 sys) [S] | 0.13 (older sources 0.6) [S] | 10.75 x 9.17 deg [S] | -18.35 | 22.6 mag/arcsec^2 |
| SMC | 62.4 kpc (203.7 kly) [S] | 2.7 [S] | 5.33 x 3.08 deg [S] | -16.28 | 23.3 |
| M31 | 765 kpc (2.50 Mly) to 778 kpc [S] | 3.44 (older editions 4.4) [S] | 3.167 x 1 deg [S] | -20.98 | 22.2 |
| M33 | 2.73 to 2.88 Mly, about 837 to 883 kpc [S] | 5.72 (effective 6.6) [S] | 70.8 x 41.7 arcmin [S] | -18.89 | 23.0 |
| NGC 6822 | 500 +/- 10 kpc [S] | 9.3 [S] | not found | -14.2 | n/a |
| IC 1613 | 730 +/- 20 kpc [S] | 9.9 [S] | not found | -14.4 | n/a |
| Fornax dSph | 143 +/- 3 kpc [S] | 9.3 (one source 8.1) [S] | not found | -11.5 | n/a |
| M81 | 3.675 +/- 0.049 Mpc [S] | 6.94 [S] | not found | -20.9 | n/a |

Mean surface brightness is `V + 2.5 log10(area in arcsec^2)` over the quoted ellipse. About 23 mag/arcsec^2 is the naked-eye limit for large objects [R], consistent with M33 being a famous difficult target [S, the "effective 6.6" remark]. Sagittarius dSph sits at 22 to 24.8 kpc in the sources found [S]. McConnachie (2012, AJ 144, 4) is the standard compilation: more than 100 galaxies within 3 Mpc with positions, distances, magnitudes and structure; a FITS table is on the author's page [S, arXiv 1204.1562]. VizieR was blocked from the research environment; it is the source to check against.

### 2.2 Placing real neighbours in the generated frame

The generated galaxy is a Milky Way analogue with a fixed Sun: R0 = 8.2 kpc at the inter-arm azimuth `model_terms(shape)['solar_angle_rad']` (5.8575 rad for the default shape [C]). Define `ex` as the unit vector from that Sun to the core and `ez = +Z`. A neighbour at galactic longitude `l`, latitude `b` and distance `d` sits at `sun + d (cos b cos l ex + cos b sin l ey + sin b ez)`, with `ey` and the sign of `b` set as follows.

- In the real galaxy longitude 90 points along the Sun's motion, which is the direction of rotation. In the generated frame the rotation at the Sun is `ez x sun_hat`, which equals `-(ez x ex)`. The first research run used `ey = ez x ex` and kept `b`; that has longitude 90 pointing against the rotation (dot product -1 [C]), a mirror image of the real sky.
- The proper embedding is a 180-degree rotation about the Sun-core axis: `ey = ex x ez` (dot product +1 with the rotation [C]) and `b -> -b`, so the project's +Z plays the real galaxy's south. M31 (l = 121, b = -21.6) then lies north of the generated plane. Distances from observers on the Sun-core line are unchanged by the choice; distances from an observer off that line (above the plane, for instance) change with it.
- Recommendation: use the proper embedding, so arms, rotation sense and the neighbours' places agree with the real galaxy.

### 2.3 How they look from other positions

Computed in the first-run frame (b kept, `ey = ez x ex`). Observers on the Sun-core line and the LMC-direction row are the same in the proper embedding; only the off-plane row depends on it. Angular sizes use the catalogue major axes; V is the catalogue value adjusted by the distance modulus [C].

| Observer | LMC | SMC | M31 | M33 |
|---|---|---|---|---|
| Sun | 49.6 kpc, 10.8 deg, V 0.13 | 62.4 kpc, 5.3 deg, V 2.70 | 765 kpc, 3.17 deg, V 3.44 | 837 kpc, 1.18 deg, V 5.72 |
| core region (1 kpc) | 49.0 kpc, 10.9 deg | 60.0 kpc, 5.6 deg | 3.15 deg | 1.17 deg |
| far side, 20 kpc | 53.2 kpc, 10.0 deg, V 0.28 | 57.7 kpc, 5.8 deg, V 2.53 | 3.11 deg, V 3.48 | 1.16 deg |
| 5 kpc off the plane, away from the Clouds' side | 52.5 kpc, 10.2 deg | 66.0 kpc, 5.0 deg | 3.16 deg | 1.18 deg |
| 40 kpc out on the far side (outside the generated galaxy) | 63.7 kpc, 8.4 deg, V 0.67 | 62.3 kpc, 5.3 deg | 3.07 deg, V 3.51 | 1.14 deg |
| a planet in the LMC's neighbourhood (60 kpc out toward it) | 10.4 kpc, 48 deg, V -3.3 | 22 kpc, 14.9 deg, V 0.45 | 3.03 deg | 1.15 deg |

On the Clouds' side of the plane (5 kpc off) the distances are LMC 47.1 kpc, SMC 59.0 kpc and M31 763 kpc [C]. In the proper embedding that side is +Z, the first-run frame's other side.

General rules, exact: angular size scales as 1/d, total flux as 1/d^2, surface brightness not at all at these distances. M31 and farther objects can therefore be one 3D position treated as nearly fixed in direction, though M31's direction itself moves up to a few degrees across the galaxy (50 kpc across 765 kpc is 3.7 degrees). The Magellanic Clouds and Sagittarius dSph (22 to 25 kpc) are 3D objects recomputed per observer exactly like stars. M31's disc inclination changes by 3 degrees at most from anywhere in the galaxy (77.0 to 74.2 degrees, from PA 37 and inclination 77 [R] inputs) [C]. From M31 the Milky Way is a 30 kpc disc (if D25 is 30 kpc [R]) 2.2 degrees across at about V = 3.5 with M_V = -20.9 [R]; a second generated galaxy at a similar distance would look the same.

## 3. The `neighbor_galaxies` table

One row per galaxy other than the host, all positions in the host galaxy's frame so nothing else needs converting.

```sql
CREATE TABLE IF NOT EXISTS neighbor_galaxies (
    id                     BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    uid                    BIGINT UNSIGNED,                  -- GEN.69 style, unique
    name                   VARCHAR(255) NOT NULL,            -- "Andromeda Galaxy", or a codec name for a generated one
    catalog_ids            VARCHAR(255),                     -- "M 31; NGC 224"; NULL for a generated galaxy
    source                 VARCHAR(16) NOT NULL CHECK (source IN ('catalog', 'generated', 'linked')),
    linked_database        VARCHAR(64),                      -- the galaxy database that holds this galaxy in full (stage 2)
    -- position in the host galaxy's frame
    position_x_pc          DOUBLE NOT NULL,
    position_y_pc          DOUBLE NOT NULL,
    position_z_pc          DOUBLE NOT NULL,
    -- shape and light
    morphology             VARCHAR(24) NOT NULL,             -- Hubble/de Vaucouleurs text: 'SA(s)b', 'dE', 'IBm', 'dSph'
    family                 VARCHAR(16) NOT NULL,             -- spiral, barred_spiral, lenticular, elliptical, irregular, dwarf_spheroidal, dwarf_irregular
    t_type                 DOUBLE,                           -- de Vaucouleurs T, -6 (E) to 10 (Im); drives colour
    absolute_mag_v         DOUBLE NOT NULL,                  -- the stored fact; apparent mag is computed per observer
    diameter_major_pc      DOUBLE NOT NULL,                  -- D25 or the V isophote
    axis_ratio_thin        DOUBLE NOT NULL,                  -- b/a seen edge on (about 0.15 spirals, 1.0 ellipticals)
    normal_x DOUBLE, normal_y DOUBLE, normal_z DOUBLE,       -- disc normal in the host frame (from PA, inclination)
    major_axis_x DOUBLE, major_axis_y DOUBLE, major_axis_z DOUBLE,   -- a reference direction in the disc plane
    color_bv               DOUBLE,                           -- or derive from t_type
    -- motion (optional, for a long-time sky)
    velocity_x_kms DOUBLE NOT NULL DEFAULT 0, velocity_y_kms DOUBLE NOT NULL DEFAULT 0, velocity_z_kms DOUBLE NOT NULL DEFAULT 0,
    seed                   BIGINT UNSIGNED,                  -- for generated ones
    UNIQUE KEY uq_neighbor_galaxies_uid (uid)
);
```

Geometry. The apparent ellipse from any place comes from the disc normal `n` and the position. With `w` the unit direction to the observer and `c = |n . w|`, the apparent axis ratio is `q = sqrt(q0^2 + (1 - q0^2) c^2)` (`q0` the edge-on thickness ratio) and the major axis is perpendicular to the projection of `n` on the sky. To fill `n` from a catalogue position angle `PA` (north through east) and inclination `i`, with `u` the direction from us to the galaxy and `N`, `E` the sky north and east at `u`: major axis `m = cos(PA) N + sin(PA) E`, `k = u x m`, `n = -cos(i) u + s sin(i) k`. The sign `s = +/-1` ("which side is near") is unknown in most catalogues; pick it from the seed for generated galaxies and record it. The construction returns exactly 77.0 degrees for M31 seen from the Sun [C]. M31's catalogue axis ratio is 0.316 against 0.268 with a 0.15 thickness [C], so the spiral thickness constant should be nearer 0.2 [R]. Mean surface brightness for visibility is computed, not stored.

`family`, `t_type` and `color_bv` set the sprite colour (ellipticals redder, about B-V 0.9; Sc and irregulars bluer, about 0.4 to 0.5 [R]); `diameter_major_pc` and `axis_ratio_thin` set size and shape; the profile (exponential disc, de Vaucouleurs bulge) is a rendering choice. Inclination and position angle come from McConnachie for dwarfs and RC3 or NED for large galaxies [R]. The `uid` follows the GEN.69 rule ([object-ids.md](object-ids.md)); the 76-bit position ID of the galaxy's centre, in Mpc units, is the alternative (the layout already has Mpc and Gpc units), with a new type code 13 for `galaxy`.

## 4. Data source: a hand-built table

A hand-built, version-controlled data file (JSON or CSV) of about 25 rows, for example `src/planetgen/data/neighbor_galaxies.json`, loaded by `planetgen plan` into `neighbor_galaxies`. The version key then covers the file by hash and Boss can edit a row. No runtime lookup of an external catalogue (no network).

Suggested contents; every value is [R] unless confirmed in section 2.1:

- Milky Way satellites: LMC, SMC, Sagittarius dSph, Fornax, Sculptor, Carina, Draco, Ursa Minor, Leo I.
- M31 group: M31, M32, M110 (NGC 205), NGC 147, NGC 185, M33.
- Other Local Group members: NGC 6822, IC 1613, IC 10, WLM.
- Beyond: M81, M82, NGC 253, Centaurus A, M83, IC 342, Maffei 1 (hidden behind the plane), the Virgo cluster as one extended object at about 16.5 Mpc.

Enter `M_V` and distance, never apparent magnitude; check every row against McConnachie 2012's VizieR table (J/AJ/144/4) and NED. Real coordinates are galactic (l, b) from the Sun, converted with section 2.2. Beyond the Local Group the (l, b) are as seen from the Sun and the distances are large enough that the Sun's place in the galaxy does not matter.

A generated neighbour (a fictional second galaxy, or filler beyond the catalogue) draws positions within 3 Mpc from the seed, `M_V` from a Schechter-like luminosity function and morphology from the morphology-density relation. Build that only when a second galaxy needs it. Catalogue galaxies keep their catalogue names; generated ones take a codec name.

## 5. One database per galaxy

### 5.1 Option 1: a `galaxy_id` column in one database

The tables keyed by the address or the singleton: `sectors` (unique `(ring_index, layer_index, ring_slot_index)`), `galaxy_shape` (`CHECK (id = 1)`), `galaxy_layer`, `galaxy_column`, `bright_stars`, `sector_stats`, `sector_paths` and `sector_path_knots`, `phenomenon_scatter`, `id_blocks`, `orbit_simulation_state`, `sector_name_registry`, `system_name_registry`, `generation_runs` and its arguments, `nearest_systems`, plus the sector column of every phenomenon table. Child tables (systems, stars, planets, moons) inherit through foreign keys, so "only the top-level rows" is the right answer for this option. But the top-level tables are the big ones (`bright_stars` 26.9 million rows at the default 1,000 L_sun floor, 60 to 220 million at 500 to 100); adding a column and rebuilding the address indexes on 10^8 to 10^9 rows is a long ALTER. Sector uids are the bare address (`geometry.provisional_sector_designation`: `ring << 33 | (layer + 4096) << 20 | slot`) and would collide between galaxies, so they would need galaxy bits. The per-database seed, version key, naming key and settings backup (OPS.13, OPS.18) would all become per-row, and the map caches and `query.galaxy_content_state` would need the id in every key.

### 5.2 Option 2: one database per galaxy (the current model)

Checked in the code:

- `db/schema.sql` describes one database per galaxy ([architecture.md](architecture.md)); `store.list_databases(prefix=...)` lists every schema on the server whose name starts with the prefix (`planetgen`, `planetgen_alpha`, ...); the API and pages accept `?db=<schema>` (`api/routes.py::_resolve_requested_db_config` validates it against that list); jobs carry `PLANETGEN_MYSQL_DATABASE` (`web/jobs.py`). No table has a `galaxy_id` column.
- The control schema (`control_schema.sql`, version 10) is deployment-global and exists so admins are not duplicated into every content schema. `galaxy_naming` (v9) has one row per galaxy database, as does `generation_size`; the version key and the galaxy seed are per database ([reproducible-galaxies.md](reproducible-galaxies.md)).
- Alembic runs per database against the same revision tree (`alembic_runner.upgrade`; head 0065 at this writing); a second database changes no revision. A second galaxy is created by `planetgen plan` on a new schema name.
- Cost: no cross-galaxy SQL joins (none are needed; neighbours are rows in `neighbor_galaxies`), the UI needs a galaxy switch, and all schemas live on one MySQL server today.

Recommendation: Option 2 for stage 2. Stage 1 adds one table and no `galaxy_id` anywhere; stage 2 sets a neighbour's `linked_database` and lets `plan` create the second schema. Migration of an existing single-galaxy database: one Alembic revision creates `neighbor_galaxies` (no change to any existing table), then `planetgen plan` or `planetgen migrate` fills it from the data file. No rewrite of any billion-row table.

## 6. The Intergalactic Frame

Keep each galaxy's own frame (origin at its core, +Z its north, +X its zero meridian) and add one level above Boss's nested frames in [navigation-frames.md](navigation-frames.md): the **Intergalactic Frame**, origin at the host galaxy's core by default (or a Local Group barycentre), +Z the host's north. Its "North" follows Boss's rule, toward the local dominant centre of mass: outside any galaxy that is the nearest galaxy centre, and in the Local Group the host.

`course_between` already takes `frame_center_gal` and `plane_up_gal`; an intergalactic course needs a new frame type and a rule for choosing it (endpoints in different galaxies, or farther apart than a few galactic radii). Between two galaxy databases the relative pose lives in one place: the host's `neighbor_galaxies` row (position, disc normal and axis in the host frame) and the mirrored row in the other database, computed from the same pose. A control-database table `galaxy_links(db_a, db_b, offset, rotation)` can hold the pair once and generate both rows. A double holds 1e6 pc (1 Mpc) without loss, and the GEN.64 position ID has Mpc and Gpc units, so a galaxy ID of type 13 fits the same 76-bit layout.

## 7. URLs, pages and generation

- Keep `?db=`. Add `/galaxies` (a list and map of neighbours and databases) and `GET /api/galaxies` (neighbour rows, with a `database` field when linked). Object references get a `galaxy:<id>` kind: `galaxy/objectref.py` has `GALAXY = "galaxy"` as the root with no id, and a neighbour is a different object, so it needs a new kind name, not the root.
- Selecting a linked galaxy changes the `db` parameter; selecting an unlinked one opens an info page only. A path prefix (`/g/<slug>/...`) is a larger change that gains nothing the query parameter does not give.
- Neighbours appear on the Galaxy Map as an extra layer when zoomed far out (billboards with the computed size and magnitude). The Generate page's plan form gets "neighbours: real / generated / none". `planetgen plan` is the only step that touches them.
- In a planet's sky, neighbours are extended sprites: an ellipse with total magnitude and angular size, shown when mean surface brightness is below about 23 mag/arcsec^2 ([sky-view.md](sky-view.md) section 2.10).

## 8. GEN.9's three questions

1. **Real or generated neighbours?** The real Local Group (about 25 rows): Boss asked for "Andromeda from Earth", the generated galaxy is a Milky Way analogue with a fixed Sun position, and the real positions are free. Offer generated neighbours as a per-plan option; a second generated galaxy is one neighbour promoted to `linked`.
2. **A galaxy id on every row, or top-level only?** Neither: with a database per galaxy there is no `galaxy_id`. If one database is insisted on, top-level only (the list in 5.1), with galaxy bits added to the sector uid.
3. **Same database, or a second one chosen at login?** A second database, chosen by `?db=` (the existing mechanism), not at login: admin sessions are deployment-wide and the pages already carry the schema in the URL. A per-session default galaxy can be a control-database setting later.

## Evidence notes

[R] items to verify when VizieR, NED, arXiv or Wikipedia are reachable:

- Every catalogue value not marked [S]: l and b of each neighbour, M31's PA and inclination (37 and 77 degrees), the M31 and Milky Way disc diameters, Hubble types, and the whole dwarf and beyond-Local-Group list (Sagittarius, Sculptor, Carina, Draco, Ursa Minor, Leo I, M32, M110, NGC 147, NGC 185, IC 10, WLM, M82, NGC 253, Centaurus A, M83, IC 342, Maffei 1, Virgo). Two recalled magnitudes proved wrong, so none of the unmarked ones should be trusted.
- The 23 mag/arcsec^2 naked-eye surface-brightness limit; the colour and thickness constants of section 3.
- The embedding of section 2.2 follows from the rotation senses stated in [galactic-potential.md](galactic-potential.md) (itself [C], arithmetic); it was checked by computing the dot product of the longitude-90 direction with the rotation direction at the model Sun, not against an external source.
- The ICRS-to-galactic matrix used for the demo positions was recalled and matches astropy to 1.2e-7.

## Sources

- Local Group and neighbours: https://en.wikipedia.org/wiki/Andromeda_Galaxy ; https://en.wikipedia.org/wiki/Large_Magellanic_Cloud ; https://en.wikipedia.org/wiki/Small_Magellanic_Cloud ; https://en.wikipedia.org/wiki/Triangulum_Galaxy ; https://en.wikipedia.org/wiki/NGC_6822 ; https://en.wikipedia.org/wiki/IC_1613 ; https://en.wikipedia.org/wiki/Fornax_Dwarf ; https://en.wikipedia.org/wiki/Messier_81
- McConnachie 2012: https://www.arxiv.org/abs/1204.1562 ; https://www.cadc.hia.nrc-cnrc.gc.ca/en/community/nearby
- Sagittarius dSph distance: https://ar5iv.arxiv.org/html/0903.3040 ; https://arxiv.org/abs/astro-ph/9605148
- Repo files read: `docs/TODO.md` (GEN.9), `docs/design/galaxy-coordinate-system.md`, `navigation-frames.md`, `object-ids.md`, `reproducible-galaxies.md`, `galactic-potential.md`, `docs/database-schema.md`, `db/schema.sql`, `db/control_schema.sql`, `db/store.py`, `db/alembic_runner.py`, `api/routes.py`, `galaxy/density.py`, `galaxy/objectref.py`, `galaxy/uid.py`, `names/object_id.py`
