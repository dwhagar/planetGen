# Galaxy Coordinate System — Design Proposal

**Status:** proposal, not implemented. This document is the design pass
`TODO.md` Phase 4 calls for ("needs its own design pass") before
`galaxyGen.py` is written. No generation code changes here — the schema
migration below is written out concretely so it can be reviewed, but it has
not been applied to `stellarObjects/schema.sql`.

**Scope:** the galaxy-scale coordinate system and the `sectors` schema
changes it needs. Explicitly **out of scope**: `galaxyGen.py` itself,
batch-generation performance, and any model of *where in the galaxy stars
are actually dense* (spiral arms, bulge, metallicity gradients). Where
those touch this design, they're flagged as open questions rather than
decided here — see §7.

## 0. Decisions already made (recap)

Per the brief, these are fixed and not relitigated below:

- Sectors are arranged **radially** from the galactic center, not on a flat
  grid/lattice.
- A sector's galactic position is defined by **spherical coordinates
  measured from the galactic center**, but **stored/exposed as xyz**
  (the sector's center point in a galaxy-frame Cartesian system).
- Sector-*internal* geometry stays exactly as it is today: milliparsecs,
  relative to the sector's own cubic center (`star_systems.position_x/y/z_mpc`,
  `sectors.edge_mpc`). This design only adds a *second*, larger-scale
  coordinate — it does not touch the first.

## 1. Coordinate system definition

### Origin and axes

- **Origin**: the galactic center, `(0, 0, 0)`.
- **Axes**: right-handed Cartesian. `X`/`Y` span the galactic plane; `+Z`
  points toward galactic north (out of the disk). This is a generation-space
  convention for this project only — it does not need to match any real
  astronomical reference frame (galactic longitude/latitude, J2000, etc.),
  since nothing here is meant to correspond to real sky coordinates.

### Spherical convention

Reuses the same `(theta, phi)` naming `spaceSector.py`'s own
`_random_unit_direction` already uses — **not** the ISO physics convention,
which swaps the two names. Kept for internal consistency with existing
code, since anyone extending this into `galaxyGen.py` will already be
reading that function:

- `r` — radial distance from the galactic center, in **parsecs** (see §2).
- `theta` — azimuthal angle in the `X`/`Y` plane from `+X`, `[0, 2*pi)`.
- `phi` — polar angle from `+Z`, `[0, pi]` (0 = galactic north pole, `pi/2`
  = in-plane, `pi` = galactic south pole).

Conversion to Cartesian (identical form to `_random_unit_direction`, scaled
by `r`):

```
x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)
```

The inverse (needed anywhere code wants "how far out, and how far off the
plane, is this sector" from a stored xyz row):

```
r     = sqrt(x^2 + y^2 + z^2)
theta = atan2(y, x)
phi   = acos(z / r)                 # undefined at r == 0; treat as phi = 0
```

### What's stored per sector

Per the fixed decision, **xyz is the persisted form** — `r`/`theta`/`phi`
are a derivation used by the placement algorithm (§3) at generation time,
not separate columns, with one exception: `r` (as `galactic_radius_pc`) is
worth persisting anyway, for the same reason `star_systems.quadrant` is
persisted even though it's derivable from `position_x/y/z_mpc` — "sectors
within radius R of the core" is an obviously common query, and a stored,
indexed column serves it directly instead of needing `sqrt(x*x+y*y+z*z)` in
every `WHERE` clause (which SQLite can't use an index for). `theta`/`phi`
have no equivalent query use, so they're not stored. See §4 for the exact
new columns.

## 2. Unit choice: parsecs (pc)

**Recommendation: store galaxy-scale distances in parsecs (`_pc` suffix
columns), not milliparsecs and not light-years.**

### Why not milliparsecs

`schema.sql`'s own header comment gives the rule that put `mpc` on
sector-internal columns in the first place: km gets unwieldy at that scale.
The same argument now cuts the other way. The Milky Way's disk radius is
commonly cited around 15 kpc (15,000 pc). In milliparsecs that's:

```
15,000 pc * 1000 mpc/pc = 15,000,000 mpc      -- and that's the radius alone
```

Two sectors on opposite sides of the galaxy would be up to **30,000,000
mpc** apart — an 8-digit number for what should be this schema's *biggest*
distance, while `sectors.edge_mpc` (~3,526) is one of its smallest. That's
the exact "ugly number" problem the brief calls out, just moved up one
level.

### Why not light-years

Light-years are the generator's own *native* unit for sector-internal
geometry (`SpaceSector.edge_ly`, `hill_radius_ly`, everything in
`spaceSector.py` before it gets converted to `mpc` at the database
boundary). It's tempting to reuse it here too, and the numbers are at least
plausible — 15,000 pc is about 48,923 ly — but `ly` and `mpc` are
**different unit families**: converting between them needs the existing
`ly_to_milliparsecs`/`milliparsecs_to_ly` helpers, which round-trip through
AU (`LY_TO_AU`, `AU_PER_MILLIPARSEC`). Every calculation that relates a
sector's galaxy-scale position to its own internal geometry — "how many
sector-widths from the core is this shell", "does this sector's cube reach
past the next shell in" — would need that AU round-trip for no reason.

### Why parsecs

A milliparsec *is* a parsec, scaled by `1/1000` — the existing convention
already anchors sector-internal geometry to the parsec family, just at
`10^-3` scale. Parsecs continue that same family at `10^0` scale:

```
edge_pc = edge_mpc / 1000          # exact, no AU round-trip
```

Worked example, `DEFAULT_SECTOR_EDGE_LY = 11.5`:

| Quantity | Value |
|---|---|
| `edge_ly` (generator-native) | 11.5 ly |
| `edge_mpc` (stored today) | ~3,526 mpc |
| `edge_pc` (proposed) | ~3.526 pc |
| Galaxy radius, same unit family | 15,000 pc |
| Galaxy radius / sector edge | ~4,254 — a plain, dimensionless ratio, no unit-family crossing |

Compare the same ratio computed via `ly`: `48,923 ly / 11.5 ly` — same
answer (~4,254), but only because you happened to convert both sides back
to the same family first; nothing about the stored columns (`mpc` vs. `ly`)
would have let you divide them directly.

### Conversion helpers needed (`stellarObjects/utils.py`)

Following the existing `ly_to_milliparsecs`/`milliparsecs_to_ly` pattern
exactly:

```
mpc_to_pc(mpc)   -> mpc / 1000                          # exact, trivial
pc_to_mpc(pc)    -> pc * 1000                           # exact, trivial
pc_to_ly(pc)     -> pc * (AU_PER_PARSEC / LY_TO_AU)     # ~3.2616 ly/pc, for display/prose
ly_to_pc(ly)     -> ly * (LY_TO_AU / AU_PER_PARSEC)     # inverse
```

`AU_PER_PARSEC` and `LY_TO_AU` already exist in `physical_constants.py`
(lines 73 and 84) — no new physical constant is needed, only the four thin
wrapper functions above, mirroring `ly_to_milliparsecs`'s own
docstring-and-one-line-body shape. `pc_to_mpc`/`mpc_to_pc` are exact
(power-of-ten scaling, no floating-point drift beyond a single multiply);
`pc_to_ly`/`ly_to_pc` exist purely for human-readable prose (e.g. rendering
a sector's summary as "~48,923 ly from the galactic core") the same way
`table_*` columns keep `ly`/`AU` display strings alongside `km`-stored raw
values elsewhere in this schema.

## 3. Radial sector-tiling scheme

### The geometric problem, stated plainly

A cube cannot tile the surface of a sphere without gaps or overlaps — this
isn't a shortcut being avoided, it's not possible in general (the reason
there's no "cubical tiling of S²," unlike the flat plane, which cubes-as-
squares tile perfectly). Any radial, shell-based layout of cubic sectors
has to either accept small seams at shell boundaries, or abandon cubes.
This proposal keeps cubes (§0 fixes `edge_mpc` as the sector shape, and
changing it is out of scope) and **explicitly accepts small gaps/overlaps
between adjacent sector cubes** as the practical resolution — justified in
"Why this is fine" below, not swept under the rug.

### Shells

Divide radial distance into shells of thickness equal to the sector edge
length, so cubes stack outward along the radial direction without gaps
*in that one direction* — the only direction cubes can tile exactly, since
adjacent shells share a flat spherical-ish boundary the way stacked plates
do:

```
shell k spans radius [k * edge_pc, (k+1) * edge_pc)     for k = 0, 1, 2, ...
nominal shell radius:  r_k = (k + 0.5) * edge_pc
```

### How many sectors per shell

Approximate each sector's footprint on the shell's spherical surface as
`edge_pc^2` (its cross-section perpendicular to the radial direction), and
divide it into the shell's surface area:

```
N_k = round(4 * pi * r_k^2 / edge_pc^2)
    = round(4 * pi * (k + 0.5)^2)          # edge_pc cancels — count depends only on k
```

This is a clean, dimensionless formula: the count per shell doesn't depend
on the absolute sector size at all, only on how many shells out you are.
Worked values:

| Shell `k` | `r_k` (sector-edges out) | `N_k` (sectors in this shell) |
|---|---|---|
| 0 | 0.5 | 3 |
| 1 | 1.5 | 28 |
| 2 | 2.5 | 79 |
| 3 | 3.5 | 154 |
| 5 | 5.5 | 380 |
| 10 | 10.5 | 1,385 |
| 100 | 100.5 | 126,923 |
| 4,253 (outermost shell reaching 15,000 pc at `edge_pc ~= 3.526`) | 4,253.5 | ~227,000,000 |

Growth is quadratic in `k` — expected, since a shell's surface area grows
with `r^2` while its thickness (and hence sector footprint) stays fixed.
Summed across all shells out to the galaxy's edge, total addressable
sectors come out close to `galaxy volume / sector volume`:

```
(4/3) * pi * 15,000^3 pc^3  /  3.526^3 pc^3   ~=   3.2 x 10^11 sectors
```

**~320 billion addressable sectors** if the whole galactic volume were
treated as a uniform sphere at this sector size. This number matters for
§7 (eager vs. lazy generation) — it's not a count anyone should try to
generate up front.

### Placing sectors within a shell

Within shell `k`, distribute its `N_k` sector centers evenly across the
sphere of radius `r_k` using a **Fibonacci (golden-angle) sphere
sequence** — a standard, deterministic way to place N points roughly
evenly over a sphere with no pole-clustering, using the same "uniform
area, not uniform angle" concern `_random_unit_direction`'s own docstring
already flags for a single random direction, just made deterministic and
extended to N points instead of 1 random one:

```
for i in 0 .. N_k - 1:
    phi_i   = acos(1 - 2 * (i + 0.5) / N_k)      # polar angle, pole-to-pole
    theta_i = (2 * pi * i / golden_ratio) mod (2 * pi)   # golden-angle azimuth step (~137.5 deg)
    center_i = r_k * (sin(phi_i)*cos(theta_i), sin(phi_i)*sin(theta_i), cos(phi_i))
```

This is deterministic and reproducible: `(shell_index, shell_slot_index)`
is a stable address that always maps to the same galaxy position, which is
useful for `galaxyGen.py` later (a shell can be regenerated or extended
without perturbing sectors already placed in other shells) — but that
determinism is itself a design choice worth flagging; see §7.

### Why this is fine despite the gaps/overlaps

- Sectors are **sparse content cells**, not solid masonry. This project's
  own stellar density (`LOCAL_STELLAR_DENSITY_LY3`, ~1 system per 352
  cubic light-years) means a default-size sector (~1,521 ly³) holds only a
  handful of systems on average (`expected_system_count`) — most of a
  sector's cube is already empty space. A seam of a few percent of an
  edge-length between two adjacent cubes has no content or gameplay
  consequence; it isn't a wall or a floor that needs to line up.
- "Nearest sector" / adjacency should be defined by **center-to-center
  distance**, exactly the pattern this codebase already uses one level
  down (`SpaceSector.nearest_neighbors`/`distance_between` for systems
  within a sector). It never needed shared-face cube topology at the
  system level, and doesn't need it at the sector level either.
- The alternative — a sector shape that *does* tile a sphere exactly — has
  no good answer for "cubic," since no polyhedron tiles a sphere's surface
  the way a cube tiles a plane. Any exact solution would mean giving up
  the fixed cubic sector shape entirely, which §0 already puts out of
  scope.

### The disk-vs-sphere gap (real galaxies aren't spherical)

Shells as defined above populate a full sphere uniformly in `phi.` A real
disk/spiral galaxy is flattened — the Milky Way's stellar disk has a scale
height on the order of ~300 pc (thin disk) to ~900-1,000 pc (thick disk),
tiny next to its ~15,000 pc radius. Treating the whole sphere as equally
"real" would put as many addressable sectors far off the galactic plane
(large `|z|`, i.e. `phi` near 0 or `pi`) as near it, which is not how the
galaxy this package is modeling actually looks.

This design deliberately **separates addressing from population**:

- The shell/Fibonacci-sphere scheme above defines the full set of
  addresses a sector *could* occupy — this is the coordinate system, and
  it's what this document is scoped to.
- Whether `galaxyGen.py` actually generates content for a given
  `(shell_index, shell_slot_index)` address is a separate, later decision
  — e.g. gating generation on `|z| / r` against a disk-density envelope,
  or a radial density falloff for the bulge/disk/halo. That's a galactic
  population/astrophysics model, not a coordinate system, and per this
  doc's scope (§0) it's flagged as an open question (§7) rather than
  decided here.

For scale, restricting to a cylinder of the galaxy's full 15,000 pc radius
but only ~1,000 pc total thickness (a generous stand-in for "thin + thick
disk") cuts the addressable volume by about 20x versus the full sphere
(`pi * 15,000^2 * 1,000` vs. `(4/3) * pi * 15,000^3`) — from ~320 billion
down to ~16 billion sectors. Still enormous, but it's the first concrete
argument that a disk envelope isn't just astronomically motivated, it's
also a large, free reduction in how much of the address space ever needs
sectors actually generated.

### Cube orientation

Each sector's cube needs an orientation, not just a center point, since
`edge_mpc` alone doesn't say which way the cube's local axes point in the
galaxy frame. The natural default needs no extra stored state at all: align
the cube's local "up" (its own `+Z`, i.e. the direction its own octant
scheme's `phi`-sign axis points) with the **radial direction** from the
galactic center to the sector's own center (`normalize(x, y, z)`) — like a
brick pointing outward, consistent for every sector, and fully derivable
from the already-stored center point. The remaining one degree of freedom
(rotation *around* that radial axis, i.e. "roll") needs a convention too;
the simplest is to derive it from a fixed reference (e.g. project galactic
`+Z` onto the plane perpendicular to the radial direction, and align the
cube's local `+X` to that projection — breaks down only exactly on the
galactic axis where a sector's radial direction *is* `+Z`, an edge case
with no sectors in practice at galaxy-core-adjacent shells).

Because this whole scheme is derivable from the center point with a fixed
convention, **no orientation column is proposed** in §4. If per-sector
random roll (rather than a fixed convention) turns out to be wanted later,
that would need a stored angle — flagged as an explicit open question in
§7 rather than added preemptively.

## 4. Schema migration proposal (v3 -> v4)

New nullable columns on `sectors`, mirroring exactly how
`star_systems.position_x/y/z_mpc`/`quadrant` are nullable ("never placed"
is a real, valid state — every sector generated by today's
single-sector-only tooling has no galaxy position at all):

```sql
ALTER TABLE sectors ADD COLUMN center_x_pc REAL;
ALTER TABLE sectors ADD COLUMN center_y_pc REAL;
ALTER TABLE sectors ADD COLUMN center_z_pc REAL;

-- Derived from center_x/y/z_pc (sqrt(x^2+y^2+z^2)) but persisted anyway,
-- same "quadrant" precedent -- lets "sectors within radius R of the core"
-- use a plain indexed range scan instead of a non-indexable expression.
ALTER TABLE sectors ADD COLUMN galactic_radius_pc REAL;

-- Which radial shell (see Sec 3) this sector's center falls in, and its
-- placement index within that shell's Fibonacci-sphere ordering. Together
-- these are a stable, reproducible "address" a shell can be regenerated
-- or extended from without perturbing sectors already placed elsewhere.
ALTER TABLE sectors ADD COLUMN shell_index INTEGER;
ALTER TABLE sectors ADD COLUMN shell_slot_index INTEGER;

CREATE INDEX IF NOT EXISTS idx_sectors_galactic_radius_pc ON sectors(galactic_radius_pc);
CREATE INDEX IF NOT EXISTS idx_sectors_shell_index ON sectors(shell_index);
```

Because SQLite's `ALTER TABLE ADD COLUMN` can't add a multi-column `CHECK`
after the fact, the "all four galaxy-position columns are NULL together or
all set together" invariant (matching the existing
`quadrant`/`position_x_mpc` null-together convention on `star_systems`)
would need to move into `schema.sql`'s `CREATE TABLE sectors` body directly
the next time this table's DDL is rewritten wholesale, e.g.:

```sql
CHECK (
    (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
    (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
    (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
)
```

`shell_index`/`shell_slot_index` are left out of that CHECK deliberately —
they're bookkeeping for the placement *algorithm*, not physically implied
by a center point the way `galactic_radius_pc` is, so a sector could in
principle have a galaxy position without having been placed by this
particular shell scheme (e.g. a hand-authored one-off position).

`PRAGMA user_version` moves from `3` to `4`. Per `db/README.md`'s existing
history format:

> **v3 -> v4**: added `sectors.center_x/y/z_pc`, `galactic_radius_pc`,
> `shell_index`, `shell_slot_index` — a sector's position in a galaxy-scale,
> spherical-coordinates-from-center layout (stored as Cartesian; see
> `docs/design/galaxy-coordinate-system.md`). NULL for a sector never
> placed in a galaxy (every sector generated by today's `sectorGen.py`
> single-sector tooling, and any migrated from an older database).
> `migrate_database`'s `_migrate_v3_to_v4` would need to be added, following
> `_migrate_v2_to_v3`'s exact shape (attach the old file read-only, copy
> every table straight across except `sectors`, which needs the explicit
> pre-v4 column list instead of `INSERT ... SELECT *` since the current
> schema has more columns than a v3 database's `sectors` table does), and
> per `migrate_database`'s own docstring, **`_migrate_v1_to_v2` and
> `_migrate_v2_to_v3` would both need updating too** to leave the new
> columns NULL, since each of those functions maps straight to the
> *current* schema in one hop rather than chaining through intermediate
> versions.

No backfill step is possible for `_migrate_v3_to_v4` beyond leaving the new
columns NULL — there is no way to recover a legacy sector's intended galaxy
position after the fact, exactly the same situation `location`/`quadrant`
are in for a `star_systems` row that was never placed in a sector at all.

## 5. Interaction with existing physics — flagged, not resolved

`stellarObjects/physical_constants.py` already has:

```python
GALACTIC_CENTER_DISTANCE_LY = 25800  # Distance from Sol to the Galactic Center in light-years
```

`Star.calculate_system_perimeter` (`starData.py:341-361`) and
`BinaryStarProxy`'s equivalent (`doubleStar.py:270-282`) both use this
fixed constant for **every** generated system, regardless of which sector
it's in — i.e., today, every star in the game is implicitly assumed to sit
at Sol's own galactocentric distance (~25,800 ly, ~7.9 kpc) when computing
its Hill sphere relative to the galaxy (`system_perimeter`). That value
then drives `hill_radius_ly`/`required_separation_ly` in `spaceSector.py`
— the minimum-separation physics used to place systems within a sector.

Once sectors have a real galactic radius (`galactic_radius_pc`), this
constant becomes **wrong in a specific, physically meaningful way** for
any sector not near 7.9 kpc out: a sector near the galactic core should
have a *smaller* Hill sphere ceiling per star (stronger galactic tidal
shear, tighter Oort-cloud limits, systems packed closer together), and a
sector out in a sparse halo shell should have a *larger* one. This is a
real interaction this design surfaces but does not resolve — changing
`calculate_system_perimeter` to take a distance parameter (falling back to
`GALACTIC_CENTER_DISTANCE_LY` when no sector position is known, for
standalone/non-galaxy-placed generation) is generation-code, out of this
document's scope, and is called out explicitly in §7 as a judgment call:
is this worth doing now, alongside the coordinate system, or deferred as
its own follow-up once `galaxyGen.py` exists and it's clear how often
sectors actually land far from 7.9 kpc?

## 6. Worked end-to-end example

A sector placed at shell `k = 50`:

```
edge_pc         = 3.526                          (DEFAULT_SECTOR_EDGE_LY)
r_50            = (50 + 0.5) * 3.526             = 178.06 pc
N_50            = round(4 * pi * 50.5^2)         = 32,047 sectors in this shell
```

Picking slot `i = 12,000` of those 32,047 via the golden-angle sequence
(illustrative, not solved to full precision here):

```
phi_i   = acos(1 - 2 * 12000.5 / 32015)   ~= 1.318 rad  (~75.5 deg from galactic north)
theta_i = (2*pi*12000 / 1.618...) mod 2*pi ~= 2.564 rad (~146.9 deg)

center_x_pc ~= 178.06 * sin(1.318) * cos(2.564) ~= -144.4
center_y_pc ~= 178.06 * sin(1.318) * sin(2.564) ~=   94.2
center_z_pc ~=  178.06 * cos(1.318)             ~=   44.6

galactic_radius_pc = sqrt((-144.4)^2 + 94.2^2 + 44.6^2) ~= 178.1   (matches r_50, as it must)
```

In light-years for a display string (`pc_to_ly`): `178.06 pc * 3.2616 ~=
580.8 ly` from the galactic core — small enough that
`GALACTIC_CENTER_DISTANCE_LY`'s fixed 25,800 ly assumption (§5) would be
off by a factor of ~44x for a system generated in this sector, which is
exactly the kind of case that flag matters for.

## 7. Open questions — for the user's judgment, not decided here

1. **Eager vs. lazy shell/sector generation.** §3 shows shell counts grow
   quadratically (~227 million sectors in the single outermost shell of a
   15,000 pc galaxy) and total addressable volume runs to ~320 billion
   sectors (or ~16 billion restricted to a disk-shaped envelope). `galaxyGen.py`
   almost certainly cannot generate "the whole galaxy" eagerly — but
   whether it should generate whole shells on demand, individual sectors
   on demand (e.g. as a player/query "visits" a region), or some other
   unit of batching is a `galaxyGen.py` design question this document
   deliberately leaves open (per TODO.md, that script's own design is
   gated on this one, not the reverse).
2. **Disk-density envelope.** §3 flags that a pure uniform sphere doesn't
   match a real disk galaxy's shape, and sketches the *concept* of gating
   generation by `|z|`/scale height, but does not specify a concrete
   profile (thin/thick disk split, bulge, spiral-arm overdensity,
   metallicity-vs-radius gradients feeding into star-type selection).
   That's a galactic astrophysical population model, arguably its own
   design pass, layered on top of the coordinate system this document
   defines.
3. **`GALACTIC_CENTER_DISTANCE_LY` becoming position-dependent** (§5) —
   whether `calculate_system_perimeter` should start taking the owning
   sector's actual `galactic_radius_pc` (converted to ly) instead of the
   fixed Sol-distance constant, and if so, whether that lands alongside
   this migration or as a separate follow-up once real galaxy-scale data
   exists to see how much it matters in practice.
4. **Cube roll/orientation.** §3 proposes a fixed, zero-extra-storage
   convention (radial-outward `+Z`, projected-galactic-north `+X`). If a
   per-sector random roll is wanted instead (e.g. for visual/narrative
   variety), that needs a stored orientation value (a single roll angle
   suffices, given the radial axis itself is already fixed by the center
   point) — not proposed here since it's pure additional complexity with
   no identified requirement yet.
5. **Deterministic Fibonacci-sphere placement vs. randomized placement.**
   §3's scheme is fully deterministic — the same `(shell_index,
   shell_slot_index)` always yields the same position, which makes shells
   independently regenerable/extensible, but also means sector placement
   has none of the organic irregularity `spaceSector.py`'s own
   system-within-a-sector placement deliberately has (Poisson-disk growth,
   `secrets`-backed randomness). Whether galaxy-scale sector placement
   should also carry some jitter (e.g. a small random offset within each
   Fibonacci-sphere cell, sized so it can't push a sector's cube into a
   neighboring shell) is an open aesthetic/design call, not a geometric
   necessity.
6. **Non-uniform sector size.** This entire scheme assumes one fixed
   `edge_pc` for every sector (today's actual default, `edge_mpc` is
   already a per-row column so it technically *can* vary — but §3's shell
   math assumes a single edge length shared by a whole shell's sector
   count formula). Larger, cheaper sectors farther from the core, or
   smaller/denser ones near it, would need a materially different
   version of §3's `N_k` formula (shell thickness and sector footprint
   would no longer be shell-wide constants). Not addressed here; assumed
   out of scope unless the user wants it folded in now.
7. **One galaxy per database.** This design assumes a single galactic
   center/origin per database file (there is no `galaxy_id` concept
   anywhere in the schema or in `TODO.md`'s framing, which only ever says
   "the galaxy," singular). If multiple independent generated galaxies
   ever need to coexist in one database, every new column in §4 would need
   a `galaxy_id` companion — flagged in case that's a real future
   requirement rather than a hypothetical one.
8. **Sector-to-sector adjacency as a queryable relationship.** §3's
   "nearest sector by center-to-center distance" answer works fine as a
   one-off computed query (mirroring `SpaceSector.nearest_neighbors`), but
   if `galaxyGen.py` or a future web UI wants persistent "neighboring
   sectors" (e.g. for travel/lore, similar to how `star_systems.location`
   persists each system's nearest neighbors today), that would need its
   own join table and its own write-time computation — not proposed here,
   since nothing in the current brief calls for it yet.
