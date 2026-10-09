# Generation determinism: drift, fingerprints, epochs and net diffs

What makes a seed rebuild the same galaxy on another Python version, operating system or CPU, measured on the repository's own code and lock file, and the design that follows: a deterministic-math helper rule with a lint, rounding before hashing, the fingerprint encoder and digest tree (GEN.58), a version key extended with a generator epoch and a battery digest (DB.6, DB.7, OPS.13 to OPS.15), the net-diff delta format (GEN.59, GEN.61) and the structure of the golden-seed test (TEST.77). It extends [reproducible-galaxies.md](reproducible-galaxies.md), which holds the overall design and the as-built notes.

Informs: GEN.55, GEN.57, GEN.58, TEST.77, OPS.12, OPS.13, OPS.14, OPS.15, DB.7, ADM.17, ADM.18, GEN.59, GEN.61, API.16, API.17 (GEN.56 is done)

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result, [C] computed by the research run against `main` at c30d74a, [R] recalled and unconfirmed. The research environment could only read search-result text, not the papers or documentation themselves, and its web-search budget ran out before the first search, so there are no [S] items here and every prior-art claim is [R]. Code paths are under `src/planetgen/`.

## 1. Summary

1. **A seed does not give the same galaxy on the project's own CI legs today.** The repo's hash-pinned `requirements.lock` was installed into Python 3.9 to 3.13 environments and `src/tests/reproducible_sector_probe.py`'s logic run over 30 sectors of 40 systems. Python 3.10 to 3.13 agree on all 30; Python 3.9 differs on all 30 [C]. The cause is `math.hypot` and `math.dist`, which Python changed between 3.9 and 3.10. Replacing them with `sqrt(x*x + y*y)` made all five agree, although the lock gives them different numpy and astropy versions [C]. CI has a 3.9 leg, so TEST.77 as written would be red there on day one.
2. **`util/draw.py` is sound** (GEN.56): identical output on Python 3.9 to 3.14 over a 20,000-round battery [C]. Only `Stream.gauss` depends on the C maths library (`log`, `cos`).
3. **Floating point across OS and CPU is a small, nonzero risk** (`exp` differs between glibc and musl in 0.07 percent of inputs, `atan2` 17.7 percent, by 1 to 11 ulp [C]); rounding to 9 significant digits before hashing hides it. Windows and macOS were not tested.
4. **The version key tracks releases, not output.** Add an integer `generator_epoch` and a battery digest (section 5). The OPS.13 seed is not changed on update (section 5.3).
5. **Fingerprint:** a 40-line canonical encoder, SHA-256, a leaf digest per sector stored at save time, a two-level tree for regions (section 4). **Net diff:** a custom address-path delta with `was` and `now` pairs (section 6).

## 2. What breaks bit-for-bit reproducibility

### 2.1 Interpreter versions

Run on Linux x86-64 (glibc 2.39) with CPython 3.9.25, 3.10.20, 3.11.17, 3.12.3, 3.13.16 and 3.14.6, plus a musl-libc 3.12.13 (python-build-standalone) for a second C maths library [C].

| Item | Result across Python 3.9 to 3.14 |
|---|---|
| `Random(seed).random()` for int, str, bytes and float seeds; `getstate()`; and `getrandbits`, `randrange`, `randint`, `choice`, `choices`, `shuffle`, `sample`, `uniform`, `gauss`, `normalvariate`, `expovariate`, `triangular` and the other variates tested | Identical |
| `paretovariate` | 3.9 differs from 3.10 to 3.14 |
| Seeding with a tuple or frozenset; `sample(set)`; `randrange(10.0)` | Work in 3.9 and 3.10 (the last to 3.11), then `TypeError` |
| Seed sign | `Random(-42)` equals `Random(42)`; `Random(1.0)` equals `Random(1)` |
| `draw.Stream` battery, 20,000 rounds | Identical hash on all six |
| `math.hypot(x, y)`, `math.hypot(x, y, z)`, `math.dist` | **3.9 differs from 3.10 to 3.14**, which agree |
| `sum()` of three or more floats | **3.11 and earlier differ from 3.12 and later** (3.12 sums floats with compensation) |
| `sum` of two floats, `math.fsum`, `statistics.fmean`, `sqrt(x*x+y*y)`, `**`, `exp`, `log`, `sin`, `atan2`, `round`, float `repr` | Identical on glibc |
| `unicodedata.unidata_version` | 13.0 (3.9, 3.10), 14.0 (3.11), 15.0 (3.12), 15.1 (3.13), 16.0 (3.14) |
| Set order of strings | Changes with `PYTHONHASHSEED`; int and tuple hashes do not |

On 200,000 random floats in 0.1 to 3 [C]: `sum` of 3 floats differs 18.4 percent between 3.11 and 3.12, `sum` of 4 floats 26.2 percent, and `math.hypot(a, b)` differs 34.5 percent between 3.9 and 3.12. A spot check on 2026-10-09 against `sqrt(a*a + b*b)` over 100,000 inputs gave: 3.9 hypot differs 37.3 percent, 3.10 and 3.12 hypot 16.8 percent (and each other 0 percent); `sum([a,b,c])` against `(a+b)+c` 0 percent on 3.9 and 3.10, 18.6 percent on 3.12 [C]. Two consequences. Swapping a version-dependent builtin for a plain-IEEE expression changes values on every interpreter, so adopting the helpers is a generation-output change and bumps the epoch (section 5). The 30-sector probe showed 3.11 and 3.12 agreeing, so no path in that sample sums three or more floats into a stored value, but `generation/system.py:1514` (`draw.uniform(0, sum(weights))`), `generation/binary.py` (`sum(s.mass ...)` for three-star systems) and `star_population.py` do sum floats, and a larger run will hit them.

About 43 `math.hypot` and `math.dist` call sites exist under `src/planetgen/` today; roughly 25 are in `galaxy/`, `population/`, `physics/` and `db/` and affect stored values, the rest are in `web/` drawing code where bits do not matter. `galaxy/geometry.py` (`sector_address_at`) already uses `sqrt(x*x + y*y)` with a comment saying why, in that one place.

### 2.2 The C maths library (OS and CPU)

CPython's `math.exp`, `log`, `sin`, `cos`, `pow` and `**` call the platform C library. IEEE 754 requires correct rounding only for `+ - * /`, `sqrt` and conversions. Measured on 100,000 random inputs per function, glibc 2.39 against musl 1.2 [C]:

| Function | glibc != musl | Max ulp apart | Function | glibc != musl | Max ulp apart |
|---|---|---|---|---|---|
| exp | 0.074% | 1 | atan | 4.25% | 1 |
| log | 0.003% | 2 | asin / acos | 6.0% / 7.6% | 1 |
| log10 | 4.37% | 11 | tanh | 7.7% | 3 |
| log2 | 0.020% | 6 | atan2 | 17.7% | 1 |
| sin / cos | 3.0% / 3.2% | 1 | pow (`**`) | 0.13% | 6 |
| tan | 3.8% | 1 | hypot (libm) | 0.000% | 0 |

Against a correctly rounded reference (mpmath, 200 bits), glibc misrounds `exp` 0.06 percent, `sin` 0.13, `pow` 0.055, but `sinh` 25 percent, `cosh` 23, `expm1` 8.5 and `log1p` 6.8 [C]. Treat `sinh`, `cosh`, `expm1`, `log1p` and `tanh` as least stable and `exp`, `log` and `pow` as most stable. The musl build's `math.sqrt` misrounded 0.09 percent of inputs (glibc 0), cause unknown: "sqrt is exact" is a property of a platform build. Windows (UCRT) and macOS libm are expected to differ from glibc at similar rates [R, medium].

Rounding hides these differences. Of 58,526 glibc/musl results that differed, the share still different after rounding to N significant digits [C]:

| N | 17 | 16 | 15 | 14 | 13 | 12 | 11 | 10 | 9 | 8 | 6 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| still differ | 100% | 58.1% | 6.4% | 0.68% | 0.075% | 0.007% | 0 | 0 | 0 | 0 | 0 |

Each dropped digit divides the chance by about 10. A quantised fingerprint misfires only when a differing value sits next to a rounding boundary; at 9 digits that is about 10^-7 per differing value in the one-ulp case, and errors that grow through a calculation need the margin. A discrete decision (`x >= threshold`) flips with probability about 10^-16 per comparison per ulp, so a whole galaxy is safe in practice. A loop that iterates to a tolerance (the Kepler solver in `physics/kepler.py`) can change its iteration count and move a result by up to its tolerance, so keep tolerances below the 9-digit quantum or accept that the digest hides it.

### 2.3 numpy, BLAS and astropy

The statement in [reproducible-galaxies.md](reproducible-galaxies.md) section 4 that no numpy or other numeric library is used was wrong. `galaxy/nebula_shape.py`, `physics/sector_path.py` and `physics/kepler.py` import numpy (the last also scipy `optimize`), `nebula_shape.py` needs scikit-image for marching cubes, and `physics/constants.py` imports astropy for 19 constants (GEN.66). The lock pins a different version per Python leg:

| Library | 3.9 | 3.10 | 3.11 | 3.12 and later |
|---|---|---|---|---|
| numpy | 1.26.4 | 2.2.6 | 2.4.6 | 2.5.3 |
| astropy | 6.0.1 | 6.1.7 | 8.0.1 | 8.0.1 |
| scikit-image | 0.24.0 | 0.25.2 | 0.26.0 | 0.26.0 |
| scipy | 1.13.1 | 1.15.3 | 1.17.1 | 1.18.1 |

- **Astropy constants drift.** Comparing every numeric constant in `planetgen.physics.constants` under astropy 6.1.7 and 8.0.1, exactly one differs: `HYDROGEN_ATOM_MASS_KG` (1.67353286206015e-27 against 1.67353286432139e-27, relative 1.35e-9) [C]. It did not affect the probe sector, but a 9-digit fingerprint sits at the edge of that difference. Freeze the constants as literals in `constants.py` and assert the astropy values equal them within tolerance in a test.
- **BLAS kernels.** numpy's matrix product goes through OpenBLAS, which picks a kernel by CPU. With `OPENBLAS_CORETYPE` set to SANDYBRIDGE or NEHALEM instead of HASWELL, `(100000,3) @ (3,3)` and `(50,3,3) @ (50,3,3)` gave different bits, and a 300x300 product also changed with thread count [C]. `np.sin`, `cos`, `exp`, `log`, `power`, `tanh`, `arctan2`, `sqrt`, `np.sum(axis=)`, `einsum` and `norm` did not change when AVX512, AVX2 and FMA3 were masked [C]; that proves little about other CPUs (Apple Accelerate, Windows builds). The repo uses `@` in `nebula_shape.py` (lines 107 and 163) and `sector_path.py` (line 232). Over 150 nebula shapes under three kernel settings, 6.5 percent of field values differed (max relative 1.6e-11) but all 150 `scale` values and all 450,000 `contains_many` answers were identical [C]: the effect is real and currently absorbed by thresholds. Replacing `@` with elementwise multiplies and adds removes the dependence for about three lines.
- **Default integer width.** numpy before 2.0 uses a 32-bit default `int` on Windows [R, high]. `nebula_shape.py` passes `dtype=np.int64` where it matters; make that a rule.

## 3. The deterministic-math rule and its lint

Rules for generation code (the packages `galaxy`, `generation`, `names`, `physics`, `population`, `util`), in order of value:

1. **Randomness only through `draw`.** Built and scanned by `src/tests/test_reproducible_draws.py`. Seed with non-negative ints only.
2. **No `math.hypot`, `math.dist` or builtin `sum()` of floats.** Add `planetgen/util/detmath.py` with `hypot(*c)` (the square root of a left-to-right sum of squares, plain IEEE operations), `dist(a, b)` and `fsum` (`math.fsum`, correctly rounded and identical in every version run). Extend `test_reproducible_draws.py`'s AST scan to flag the three, with an allowlist for integer sums and a comment. `math.fsum` costs more than `sum`; check the performance benchmarks, none of these are hot paths. Switching changes stored values, so it ships with an epoch bump.
3. **Transcendentals are allowed, exact float comparison on their output is not.** Avoid `sinh`, `cosh`, `expm1`, `log1p` where `exp` and `log` do the job. A libm-free normal draw (a ziggurat, or Wichura's AS241) would remove `Stream.gauss`'s dependence on `log` and `cos`; not worth doing unless an exact cross-libm guarantee is wanted.
4. **No BLAS in generation.** Replace `@`, `dot`, `matmul` and `linalg.*` in `nebula_shape.py` and `sector_path.py` with elementwise code, or show the quantity feeds only quantised output. Give numpy an explicit dtype. numpy stays for display meshes.
5. **Freeze astropy-derived constants as literals.**
6. **Never iterate a set or dict built from strings or tuples without sorting**, and never `os.listdir` or `glob` unsorted. The scan finds no `list(set(` today.
7. **Text.** Keep names and word lists ASCII or handle them as bytes: the Unicode database differs from 13.0 (3.9) to 16.0 (3.14), so `.lower()`, `.title()` and `casefold()` can change for non-ASCII letters [R for the effect, C for the versions].
8. **Files and float text.** Write canonical bytes in binary mode with `\n` only and hash the canonical body, not the file. Float `repr` does not call the platform printf [R, high]; the six Linux versions agree, Windows was not run. No locale-aware formatting.
9. **Run the no-database battery (section 7, tier A) on the Windows and macOS CI jobs too.** `ci.yml` has `windows-jobs` and `macos-installers` jobs but the test matrix (3.9, 3.12, 3.13) is Linux only. GitHub's `macos-latest` is Apple Silicon [R], a different CPU family and libm, the strongest single cross-check.
10. **A leg that still differs after rules 2 to 5 is not loosened silently.** The failure report says whether integers and strings or only floats differ; a float-only difference on one platform is a finding to triage, and `epochs.json` can hold a per-platform alternative digest as a last resort.

## 4. Canonical serialization and fingerprints (GEN.58)

### 4.1 Encoder

RFC 8785 (JCS) sorts keys, writes UTF-8 with no whitespace and uses ES6 number serialization [R, high]. Against `json.dumps` (pure-Python `rfc8785` 0.1.4 [C]) it writes `1.0` as `1`, `1e16` as `10000000000000000`, `-0.0` as `0`, and errors on integers of 2^53 or more and on nan.

So JCS is not what `json.dumps` writes, cannot carry 128-bit seeds, 96-bit uids or BIGINT ids as numbers, and erases the sign of zero. Because the fingerprint quantises every float first, it needs a rule stable across Python 3.9 to 3.14, not cross-language portability. A 40-line encoder inside `planetgen`, documented as "JCS-like, not JCS", with no new dependency:

- **Floats:** reject nan and inf; `x = 0.0 if x == 0 else x`; quantise with `float(format(x, ".8e"))` (9 significant digits); `json.dumps` then writes the shortest repr of that double.
- **Integers:** decimal; anything above 2^53 as a hex string.
- **Strings:** UTF-8 as stored, no Unicode normalisation (the tables differ by Python version).
- **Rows:** keys sorted (`sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False`); lists in generation order, never database id order; ids, foreign keys and timestamps skipped.
- **Format tag:** `b"pgfp1\0"` in the hash input, so a rule change (digits, field set) is told apart from a content change. `fp_spec` is kept separate from `generator_epoch` in the golden file.

MySQL 8.4 and MariaDB 11.4 may send DOUBLE values in different text forms (17 significant digits against shortest) [R, medium; not tested]. A digest read back from the database is independent of that only because of the 9-digit quantisation.

### 4.2 Hash function

SHA-256, BLAKE2b and SHA3-256 are in the Python 3.9 standard library; BLAKE3 is not (PyPI `blake3` 1.0.11, released 2026-10-08, needs Python 3.11 or later; 1.0.10 needs 3.8 or later). SHA-256 is already used for unit seeds, uids and the naming key. Measured on the research machine (Python 3.9.25, OpenSSL 3.5.4, no SHA-NI): SHA-256 130 to 310 MB/s, BLAKE2b about 270 MB/s, BLAKE3 2 to 4 GB/s [C]. A 40-system sector generates in about 1.3 s; its canonical text is on the order of 0.5 MB, so hashing is a few milliseconds, under 1 percent of generation. Canonicalisation in Python is the cost, not the hash. **SHA-256.**

### 4.3 Per-sector digests and the region tree

Sector addresses have a total order: `geometry.provisional_sector_designation` and `galaxy/uid.sector_uid` pack ring, biased layer and slot into one integer, ring in the high bits.

- **Leaf** = `SHA256(0x00 || fp_spec || uid (8 bytes) || canonical_sector_bytes)`, stored per sector in the galaxy database (one 32-byte column, or 8 bytes truncated, since it detects change and is not a security feature). DB.9's content checksum is this leaf.
- **Node** = `SHA256(0x01 || for each child in uid order: child key (8 bytes) || child digest)`. The distinct prefix bytes stop a leaf from being read as a node, as in RFC 6962 Merkle trees [R].
- **Levels:** sector leaf, (ring, layer) node, ring node, galaxy root. A region is any ring or layer range or list of sectors; its digest covers the leaves or nodes it spans plus their count, so two regions with different sector sets never coincide. An ungenerated sector is absent and the count says so.
- **Cost** [C]: a flat digest over 1,000,000 stored leaves took 1.0 s in Python (the loop, not the hash). With a stored node table, changing one sector rehashes one 1,000-child group and the 1,000-entry top in 0.55 ms. Compute the leaf at save time from the in-memory sector, which does not depend on the database float round trip; keep a `sector_digest_nodes` table keyed by (ring, layer) marked dirty on save or edit; compute ring and galaxy roots on demand. At 12 billion sectors a full-galaxy digest is a batch job; a region of up to about 10^6 sectors is interactive.

`tests/galaxy_fingerprint.py` is a GEN.39 stand-in using `repr` of full-precision floats and decorated names; GEN.58 replaces it. `tests/reproducible_sector_probe.py` compares runs under one interpreter only, which is why it could not see cross-version drift until it was run under five.

## 5. Version key, epoch and battery digest

### 5.1 What the current key gets wrong

`galaxy/version_key.py` builds 22 hex digits (release, Python major.minor.micro, OS digit, architecture digit; see [reproducible-galaxies.md](reproducible-galaxies.md) section 4).

1. **It tracks releases, not output.** By `changes/README.md`, every `patch` or `minor` note bumps REVISION and BUILD is the sum of TODO category counters, so the key moves when items are added to `docs/TODO.md`. OPS.14 would warn on most updates and OPS.12's "refuse when the running release isn't Y" would refuse almost always, though output is the same.
2. **It over-specifies the environment.** Python micro (3.12.3 against 3.12.4) does not matter. Python minor mattered in 3.10 (hypot) and 3.12 (sum), both fixable in code. OS and architecture matter only through the C maths library and BLAS, which the digits do not name.
3. **It omits what did differ in the tests:** numpy, astropy, scikit-image and scipy versions, which differ per Python leg in the lock. OPS.13's lock hash covers that; the logged key does not.

### 5.2 Recommended record

Keep the 22-digit key as provenance (log line, DB.6 columns) and add:

| Field | Meaning | Changes when |
|---|---|---|
| `generator_epoch` (int) | Algorithm epoch | A PR changes the golden digest, and nothing else |
| `fp_spec` (int) | Fingerprint format version | The canonical encoding or quantisation changes |
| `schema_version` | Galaxy schema (Alembic head) | Migrations (already recorded) |
| `battery_digest` (hex) | Digest of the canned battery run in this environment | Any output change, including platform or library drift |
| `env` | Python x.y.z, OS, machine, `platform.libc_ver()`, numpy, astropy, scipy, scikit-image versions, lock hash | Informational |
| `release` | MAJOR.REVISION.BUILD | Every release |

`src/tests/golden/epochs.json` maps each epoch to its first release and its golden battery digest, so a tool can name the release to check out for any epoch (OPS.12) and OPS.14 can say "epoch 7 was introduced in 7.131.412".

**Enforcing the epoch.** `scripts/bump_version.py --check-pr` already compares a PR with `main`. Add the rule: the golden file changed if and only if the epoch is bumped if and only if the `changes/` note carries a `generation-output: changed` line. This is the mechanism Python's `.pyc` magic number and numpy's stream-compatibility policy (NEP 19) use [R].

**Comparison rule** (OPS.14, OPS.12). A stored galaxy is reproducible here if its epoch matches and the battery digest computed now equals the stored one. Epoch matches, battery differs: "same code, different platform or libraries; the first differing battery case is X" (the case that found the 3.9 hypot drift). Epoch differs: say which release introduced it. Release, Python micro, OS and architecture differences are notes, not errors.

**Per-sector version (DB.7).** The TODO text asks for four columns per sector (`version_key`, version, python, platform). At 12 billion sectors that is repeated text. Recommended: the `generation_runs` table that already exists (`db/schema.sql`) as the reference, and per-sector columns of only `generator_epoch` (SMALLINT) and `run_id` (INT, foreign key), the full key text living on the run row. A mixed-version galaxy is then a sector-by-sector comparison of epochs.

### 5.3 The OPS.13 question: "recalculate the seed value"

Boss (2026-10-02 02:08Z): each update "will have to make a sweet to recalculate the seed value based on the current code version and it'll need to keep the seed/version combo for say the last 10 updates."

**Do not change the seed; recalculate what is derived from the code and record it next to the seed.**

1. The galaxy seed identifies the galaxy. The stored rows came from it; changing it makes the stored seed claim a galaxy that does not exist.
2. Mixing the version into the seed (`H(seed || epoch)`) would re-roll every sector at each epoch bump, including sectors whose generator code did not change, which defeats "a release that changes one formula changes only what that formula touches" ([reproducible-galaxies.md](reproducible-galaxies.md) section 3).
3. What is recalculated on each update and kept in the history row is the new version key and epoch, the dependency fingerprints (lock hash, numpy and astropy versions, libc) and the battery digest under the new code (OPS.15). That is the sweep: does the same seed still produce the same thing on this code.

An OPS.13 row is `(galaxy_seed, version_key, epoch, battery_digest, lock_sha256, env_json, date)`, last 10 per galaxy, plus the corpus, offensive-words and name-list hashes already planned (the filtered list the generator actually uses is hashed, so a corpus download difference that leaves the list unchanged does not alarm). A `derived_run_id = SHA-256(seed || version_key)` could label a row but must never feed generation; not recommended.

### 5.4 The battery (OPS.15)

A pure function `generator_battery_digest()` in the package (not a script), no database and no network, running a fixed canned battery and returning an overall digest plus per-case digests:

- 6 to 8 tiny cases, each a fixed (galaxy seed, address, config) giving one 3-system sector through `generation/run_sector.generate_sector` under `draw.bound(unit_seed(...))`; one rogue planet or phenomenon; one nebula shape (`draw_shape`, `scale`, 50 `contains` points); one bright-star scatter on a tiny region; one population stream; one `unit_seed` known answer. Each is hashed with the section 4 encoder.
- Cost [C]: eight 3-system sectors take 1.2 s; importing the package costs 2 to 3 s of CPU. So 3 to 5 s per update. `update.sh` and `update.ps1` run it with the new code and compare with the previous history row; ADM.17 and API.16 can show "battery OK" or "differs" by running it at start-up or on demand.
- It detects an intentional change (the epoch should have been bumped), an accidental one (the golden test should have failed first) and environment drift with the epoch unchanged. It cannot see code paths it does not exercise; grow it as GEN.42, PERF.18, API.12 and API.13 land.
- The pure part (seed, address, config to sector object) is also what API.17's remote generation needs; keep it free of database reads. The entry points are `generation/run_galaxy.generate_and_save_sector_at` (`galaxySeed.seeded(galaxy_seed, "sector", address)`), `galaxy/nebula_field.py` (`"nebula-cell"`) and `generation/run_sector.generate_sector`.

### 5.5 How other systems handle it (all [R], unconfirmed)

Minecraft (Java) stores a data version per chunk, so old chunks keep their terrain and new chunks use the new generator, with 1.18 blending the seams (medium-high): DB.7's per-sector version in miniature. Factorio (exchange string with settings, seed and game version), Terraria (1.4 changed generation) and Dwarf Fortress (separate seeds for world, history and names; a saved world is a snapshot) tie a seed to a version (medium). Brogue and DCSS show the version beside the seed (medium). No Man's Sky stores player edits on top of deterministic seeds; how updates preserve them is unknown (low, do not rely on it). numpy's NEP 19 freezes `RandomState`'s stream while `Generator`'s may change, and Rust `rand`'s `StdRng` may change algorithm while a named ChaCha generator is fixed [R, high]: "the stream may change" is kept apart from "this named thing never changes", which `draw.py` does and an epoch counter formalises. None promises that a seed rebuilds across versions, which matches Boss's "no automatic migration".

## 6. Net-diff and edit-log format (GEN.59, GEN.61)

### 6.1 Candidates

RFC 6902 JSON Patch addresses arrays by position, is an ordered script, and fails safe on a changed base only with added `test` ops. RFC 7396 JSON Merge Patch replaces arrays whole, cannot set a field to null (null deletes it) and checks nothing. A custom delta addresses objects by stable generation index, holds one final state per object (order independent), marks deletes as tombstones and checks the base through a `was` value.

Demonstrated with `jsonpatch` and a hand-written RFC 7396 [C]. A diff of one planet mass in system 1 is `[{"op":"replace","path":"/systems/1/planets/1/mass","value":99.0}]` (71 bytes). If a later generator version inserts a system at the front, the same patch applies without error to the wrong planet; adding a `test` op first makes it fail loudly. The merge patch for the same edit has to carry the whole `systems` array (130 bytes in a 2-system toy, far larger in real data). The custom form was 92 bytes with no position in it.

### 6.2 The format

A custom document built on merge-patch ideas but keyed by stable object address, one entry per changed object, flat field maps:

```json
{
  "format": "pgdelta/1",
  "epoch": 7,
  "orbit_epoch": "2026-10-09T00:00:00Z",
  "objects": {
    "S12.3.0/y4/p2":    {"set": {"mass_kg": {"was": 5.97e24, "now": 6.1e24}}},
    "S12.3.0/y4/p2/m1": {"deleted": true},
    "S12.3.0/y5":       {"regen": "9F3A07C2E81B44D5A1C06E7B3D2F9081"},
    "S12.3.0/y6/p+3F0A91": {"added": {"class": "ice giant"}}
  }
}
```

Rules, which make the GEN.61 merge deterministic and idempotent:

1. **Address path:** sector `S<ring>.<layer>.<slot>`, then `y<n>` system, `p<n>` planet, `m<n>` moon, `b<n>` belt, `c<n>` comet, where `n` is the **generation index**: the object's position in generation order when first produced, stored and never renumbered. Deleting object 2 leaves 0, 1, 3. Admin-added objects get `+<random 24-bit hex>` so they cannot collide with generated indices. `galaxy/uid.py` and `store.assign_uids` derive a uid from the object's rank among its parent's rows by row id; if an admin deletes a planet and the sector is later re-saved, that rank can shift. The path helper therefore needs a stored generation-index column, not row rank and not the uid.
2. **`set` holds field to `{was, now}`.** `now` is the admin's value; `was` is the generated value captured before the first change. A field whose latest `now` equals its `was` is dropped, and an entry with no fields is dropped, without regenerating (this is what GEN.59's "an edit that puts a value back drops its entry" needs). It also gives rebase checks and a readable log.
3. **`regen` holds the new 128-bit seed.** A regenerate at path P clears every entry at or under P and writes `{"regen": seed}`; later edits are recorded on top and their `was` values refer to the regenerated content. A delete at P clears everything under P and writes the tombstone. Both are prefix operations: the merge sorts by path and processes prefix-first.
4. **Order independence:** the file maps path to final state, so folding pending deltas is "for each pending delta in time order, overwrite by rule 2 or 3", and a different grouping gives the same file. Each pending row is `(time, path, op, payload)`.
5. **Canonical bytes:** sorted keys; `now` and `was` as the exact Python float repr, not quantised (they are stored values; the fingerprint quantises when comparing); seeds and ids above 2^53 as hex strings; binary write, `\n`, no BOM; a `content_sha256` of the canonical body inside the file, so GEN.61's "written and read back" check is one comparison.
6. **When the base changed (different epoch):** by Boss's rule there is no automatic migration, so the default is to refuse; OPS.12 stops at the epoch check and names the release. A separate explicit rebase tool, if ever wanted, is a three-way compare: base is `was`, ours is `now`, theirs is the regenerated value. Equal `was` applies; different `was` is a conflict (keep the admin's value, mark it conflicted); a missing path is an orphan (listed, not applied); a `regen` is kept but flagged because new code gives different content from what the admin saw. Never apply silently.

**Regenerate needs a seed argument.** GEN.59's text says it replaces today's `random.seed()` in `api/edits.py`. That call is gone: a search of `src/planetgen` finds no `random.seed`, and GEN.56 moved the draws to `util/draw.py`. What remains is that `admin/edits.py: regenerate_planet`, `regenerate_moon` and `regenerate_belt` draw from whatever stream is bound, the run stream in a request, with no stored seed. GEN.59 must add a seed argument to those functions, bind `draw.bound(seed)` around them and store the seed.

### 6.3 The settings JSON (ADM.18) in the same style

The same encoder: `seed` as 32 uppercase hex, `version_key` with parts spelled out, `epoch`, `fp_spec`, every creation setting at full float repr (they are generator inputs and must not be quantised), `lock_sha256`, `env`, the naming key and its `CODEC_VERSION` ([object-ids.md](object-ids.md), GEN.70), the filtered word list's SHA-256 with the list itself in gzip and base64 (stars and sectors keep word-salad names per Boss's 2026-10-08 decisions), and `content_sha256`. The file name `<seed>-<key>-<YYYYMMDD>-<HHMMSS>Z.json` has no colons and is valid on Windows.

## 7. Golden-seed test structure (TEST.77)

### 7.1 Tiers

Three tiers, cheapest first, so a failure is localised before the expensive test runs.

- **Tier A, no database, every leg including Windows and macOS (3 to 5 s).** The OPS.15 battery. The golden record per case is a `structural_digest` (ints, strings, booleans, counts, names) and a `float_digest` (floats at 9 digits), each over named sections (system, star, planet, moon, belt, comet, phenomenon, population), plus the case's overall digest. On mismatch the report prints a table: case, section, expected, actual, structural or float. A structural mismatch is a code change; a float-only mismatch on one platform is drift.
- **Tier B, database, 1 against 4 workers.** TEST.77 as written: plan, a few sectors, a scatter and a backfill. Two assertions: 1 worker equals 4 workers (always true, no golden needed), and the per-sector digests equal the golden's. Print differences by address (`S12.3.0: expected 9fa1..., got 02be...`), then the region digest. A few dozen sectors.
- **Tier C:** the existing hash-seed and locale subprocess probe (`test_a_sector_does_not_depend_on_hash_seed_or_locale`), unchanged.

One "Rosetta" case is stored as full canonical text (a 1-system sector, a few KB, `src/tests/golden/rosetta_epoch_N.txt`) so a failing diff reads line by line; the other cases store digests only. On failure the actual text is written next to the golden path and `diff -u` is printed.

### 7.2 Updating goldens

- One tool, `python -m tests.golden_update --reason "..."` (or `pytest --update-golden`), rewrites `src/tests/golden/` and `epochs.json`, bumps `generator_epoch` and records the reason. It refuses to run when the `CI` environment variable is set and when the diff is empty.
- `bump_version.py --check-pr` enforces golden file changed, epoch bumped and note line `generation-output: changed` together (section 5.2).
- The golden file records the environment it was made on (`python`, `platform`, `numpy`, `astropy`) as information only; the digest is the same on all legs. This replaces the TODO's "pinned per release with the version it was made on": per-leg goldens would hide the drift of section 2, and a per-release pin would be rewritten on every release because REVISION moves each time.
- Order of work: the deterministic-math items (section 3, rules 2 to 5) land before TEST.77, or the 3.9 leg is red.

## Decisions already taken

Boss's recorded decisions that this note builds on; the rest of the note is recommendation.

- Seed: 128 bits, stored in the database, logged at the top of every generation (Boss, 2026-10-02 01:40Z).
- A seed reproduces a galaxy only on the exact setup that made it; there is no automatic migration of old galaxies to a new release's output (Boss's rule, recorded in [reproducible-galaxies.md](reproducible-galaxies.md) sections 2 and 10).
- Admin changes are a net difference from the generated galaxy, a regenerate draws a fresh 128-bit seed for that object, and the JSON file changes only with the deltas at the end of the day (Boss, 2026-10-02 02:20Z and 02:28Z).
- Each update keeps the seed/version combination for the last 10 updates (Boss, 2026-10-02 02:08Z).
- Generation draws reach the unit's stream through a `ContextVar`, not an `rng` argument on every function (Boss, 2026-10-09 07:02Z).
- Stars and sectors keep word-salad names; planets, moons and belts keep the "<system> I" pattern (Boss, 2026-10-08 03:57Z and 04:00Z).

## Evidence notes

Computed from scripts in the research run's `exp/` scratchpad. The hypot and sum rates marked as a spot check in section 2.1 were rerun on 2026-10-09.

[R] items to verify when paper and documentation access is allowed (the research environment could only read search-result text, not the papers or documentation):

- That the Python documentation promises reproducibility only for `random()` (it is what `draw.py`'s docstring says).
- That 3.12's `sum()` uses compensated summation and that 3.10 changed `math.hypot`: the behaviour was measured; the attribution to those releases' "What's New" is from memory.
- Windows (UCRT) and macOS (libm, Accelerate) behaviour: nothing was run on them.
- That `macos-latest` on GitHub is arm64; that numpy before 2.0 defaulted to a 32-bit `int` on Windows; that Python's float `repr` and `format` do not call the platform printf; the MySQL against MariaDB DOUBLE text formatting claim.
- All prior-art statements in section 5.5; RFC 6962 domain-separation prefixes; RFC 8785's rules (the library's behaviour was tested, the RFC text was not read).
- The musl `sqrt` misrounding is one build's observation, undiagnosed.

## Sources

- Repo files read: `docs/TODO.md` (the items listed), the two design notes, `util/draw.py`, `galaxy/seed.py`, `version_key.py`, `uid.py`, `geometry.py`, `nebula_shape.py`, `physics/constants.py`, `kepler.py`, `sector_path.py`, `admin/edits.py`, `src/tests/test_reproducible_draws.py`, `reproducible_sector_probe.py`, `requirements.lock`, `.github/workflows/ci.yml`, `scripts/bump_version.py`, `changes/README.md`.
- PyPI JSON API: https://pypi.org/pypi/blake3/json (1.0.11 needs Python 3.11 or later, 1.0.10 needs 3.8 or later), https://pypi.org/pypi/rfc8785/json (0.1.4, pure Python, Trail of Bits); also `canonicaljson` 2.0.0, `jcs` 0.2.1, `xxhash` 4.0.1.
- Runtimes used: CPython 3.9.25, 3.10.20, 3.11.17, 3.12.3, 3.13.16, 3.14.6 (glibc 2.39) and 3.12.13 musl (python-build-standalone via `uv`); mpmath 1.4.1, numpy 2.5.3, rfc8785 0.1.4, jsonpatch.
- No web pages were consulted.
