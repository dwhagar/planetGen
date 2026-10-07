# Open-Source Library Migration and Refactoring Roadmap for planetGen

## Executive Summary and Architectural Principles

The planetGen codebase currently maintains substantial custom implementations for capabilities that are standard in the open-source Python ecosystem. These include custom RFC 6238 TOTP algorithms, QR code matrix generators, Kepler orbital solvers, raw SQL persistence frameworks, custom multi-process work queues, rate limiters, file-backed caches, and regex-based Markdown converters.

Reinventing these domain utilities increases the managed surface area to over 35,000 lines of Python code, introducing security hazards, numerical precision limitations, and elevated maintenance costs. Modernizing planetGen through strategic adoption of mature third-party open-source libraries will significantly reduce code footprint, enhance security compliance, improve numerical accuracy, and boost execution speed via C-optimized native bindings.

Architectural refactoring is structured around four primary criteria:

1. Security and Standard Compliance: Offload cryptography, authentication, and token handling to audited, industry-standard packages.
2. Computational Performance: Leverage C/Fortran bindings (NumPy, SciPy) for intensive mathematical and astrodynamic calculations.
3. Operational Reliability: Replace hand-rolled multi-processing queues and database connection managers with battle-tested frameworks.
4. Maintainability: Delete single-purpose internal helper modules in favor of declarative schema models and standardized utility libraries.

## Domain Analysis, Alternative Library Evaluations, and Code Impact

### 1. Security, Two-Factor Authentication, and QR Code Generation

#### Current Codebase State

* `src/planetgen/admin/totp.py` (102 lines): Custom implementation of HMAC-SHA1 Time-Based One-Time Passwords (RFC 6238).
* `src/planetgen/admin/qrcode.py` (907 lines): Single-file copy of Project Nayuki's QR Code generator, handling low-level bit buffers, Reed-Solomon error correction, and SVG matrix generation.

#### Recommended Stack

* `pyotp`: Standard library for TOTP/HOTP validation.
* `qrcode` (or `segno`): Industry-standard QR code generation.

#### Alternative Libraries Evaluated

* `segno`: Pure Python, zero-dependency QR code encoder supporting micro QR codes and direct SVG rendering.
* `authlib`: Comprehensive suite covering OAuth 1.0/2.0, JOSE, JWT, and TOTP.

#### Pros, Cons, and Code Impact

* `pyotp` Pros: Pure Python, zero external dependencies, strict adherence to RFC 6238 and RFC 4226, fully tested against clock drift and replay attacks.
* `pyotp` Cons: Requires explicit timestamp formatting for customized step intervals.
* `qrcode` Pros: Widespread adoption, native integration with PIL/Pillow for raster rendering and vector SVG output.
* `segno` Pros: Faster SVG generation than `qrcode`, lighter weight, zero third-party dependencies.
* Code Impact: Removes ~1,000 lines of complex bit manipulation and matrix math. Replaces custom TOTP and QR functions with 5 lines of declarative library calls.

### 2. Astrodynamics, Orbital Mechanics, and Physical Constants

#### Current Codebase State

* `src/planetgen/physics/kepler.py` (408 lines): Custom bisection and Newton-Raphson solvers for Kepler's equation (elliptical orbits) and Barker's equation (parabolic orbits).
* `src/planetgen/physics/constants.py` (447 lines): Hardcoded astronomical, physical, and conversion constants.
* `src/planetgen/physics/planets.py` (1,154 lines): Custom calculations for orbital period, mean anomaly, and gravitational parameters.

#### Theoretical Mathematical Models

The orbital propagation logic in `keplerMotion.py` models two-body celestial dynamics using the following governing equations:

Kepler's Equation for Elliptical Orbits ($0 \le e < 1$):
$M = E - e \sin(E)$

True Anomaly ($f$) derived from Eccentric Anomaly ($E$):
$\tan\left(\frac{f}{2}\right) = \sqrt{\frac{1 + e}{1 - e}} \tan\left(\frac{E}{2}\right)$

Barker's Equation for Parabolic Orbits ($e = 1$):
$M_p = \frac{1}{2} \tan\left(\frac{f}{2}\right) + \frac{1}{6} \tan^3\left(\frac{f}{2}\right)$

Two-Body Gravitational Parameter ($\mu$) in astronomical units ($AU$, solar masses $M_\odot$, years $yr$):
$\mu = G (M_1 + M_2) = 4\pi^2 M$

#### Recommended Stack

* `scipy`: Numerical root finding (`scipy.optimize.root_scalar`) and vectorized array operations.
* `astropy`: Standardized physical constants (`astropy.constants`) and unit transformation framework (`astropy.units`).

#### Alternative Libraries Evaluated

* `skyfield`: High-precision astrodynamics library designed for planetary position calculation using NASA JPL ephemerides.
* `poliastro`: Specialized orbital mechanics library built on top of Astropy and Numba.

#### Pros, Cons, and Code Impact

* `scipy` Pros: C-optimized vector root finders solve Kepler's equation across thousands of orbital bodies concurrently; eliminates custom bisection loops.
* `scipy` Cons: Large binary footprint (~30 MB wheel size).
* `astropy` Pros: Guarantees unit correctness across metric, imperial, and astronomical coordinate systems; provides IAU-compliant physical constants.
* `astropy` Cons: Slightly higher import latency; unit conversions add minor execution overhead if not vectorized.
* `skyfield` / `poliastro` Pros: Extremely accurate real-world ephemeris models.
* `skyfield` / `poliastro` Cons: Overkill for procedurally generated fictional stellar systems.
* Code Impact: Reduces numerical iteration time by up to $10\times$ via vectorized SciPy calls; replaces manual constant definitions with `astropy.constants`.

### 3. Database Persistence and Schema Migrations

#### Current Codebase State

* `src/planetgen/db/store.py` (9,793 lines): Hand-rolled database abstraction layer utilizing raw PyMySQL, custom connection pooling (`dbutils`), manual SQL string construction, and custom transactional batching scopes.
* `planetgen.cli.migrate` (172 lines): Custom script tracking migrations via raw SQL script files stored as compressed `.sql.gz` fixtures (`schema_v14.sql.gz` to `schema_v52.sql.gz`).

#### Recommended Stack

* `SQLAlchemy`: Industrial-strength SQL toolkit and Object Relational Mapper (ORM).
* `Alembic`: Database migration framework integrated with SQLAlchemy.

#### Alternative Libraries Evaluated

* `Peewee`: Expressive, lightweight ORM with simple setup and small footprint.
* `Tortoise ORM`: Async-first ORM modeled after Django's syntax.
* `DuckDB` / `Polars`: In-process analytical database for high-speed sector queries.

#### Pros, Cons, and Code Impact

* `SQLAlchemy` Pros: Complete abstraction over MySQL/MariaDB differences, built-in connection pooling, automatic statement preparation, prevention of SQL injection, declarative schema definitions.
* `SQLAlchemy` Cons: Steeper learning curve for complex ORM query optimization.
* `Alembic` Pros: Programmatic Python-based schema migrations, auto-generation of migration scripts via model diffs, complete replacement of compressed raw SQL files.
* `Peewee` Pros: Easy transition for smaller projects, lighter cognitive overhead.
* `Peewee` Cons: Less robust migration engine compared to Alembic when managing complex foreign keys and constraints.
* Code Impact: Deletes up to 8,000 lines of manual query formatting, parameter binding, cursor lifecycle management, and schema detection logic in `store.py`.

### 4. Background Processing and Work Queue Management

#### Current Codebase State

* `src/stellarObjects/workQueue.py` (1,611 lines): Custom process pool supervisor managing independent generation units (sectors, star scattering) with database lease locks (`work_lease` table), heartbeats, stale run reclamation, and custom process IPC.
* `src/jobRunner.py` (305 lines): Custom process tree runner executing background jobs started by the web administration interface.

#### Recommended Stack

* `RQ` (Redis Queue): Simple, lightweight background task processing built on Redis.
* `Huey`: Lightweight alternative supporting SQLite or Redis backends.

#### Alternative Libraries Evaluated

* `Celery`: Full-featured, distributed task queue supporting multiple message brokers (RabbitMQ, Redis) and complex workflows.
* `APScheduler`: Advanced Python scheduler for cron-like background jobs.

#### Pros, Cons, and Code Impact

* `RQ` Pros: Minimal setup, native Python decorator integration (`@task`), built-in process monitoring, job retries, timeout enforcement, and status tracking.
* `RQ` Cons: Requires a Redis server instance.
* `Huey` Pros: Can run directly on SQLite or file-based backends without needing an external service like Redis.
* `Celery` Pros: Multi-broker flexibility, complex task chaining and routing.
* `Celery` Cons: High configuration complexity and heavy operational overhead for single-node deployments.
* Code Impact: Replaces ~1,900 lines of complex IPC, database locking, multi-process signaling, and supervisor polling loops with a standard `@task` queue architecture.

### 5. Rate Limiting, Throttling, and Brute-Force Guards

#### Current Codebase State

* `src/planetgen/admin/throttle.py` (393 lines): Dual-counter lockout tracking (per client IP address and per username) backed by custom database and memory stores.
* `src/html/api/limiter.py` (108 lines): Flask rate-limiter setup wrapper.
* `src/html/api/loginguard.py` (156 lines): Middleware guarding password checks against brute-force attacks.

#### Recommended Stack

* `Flask-Limiter` / `limits`: Production-ready rate limiting extension for Flask applications.

#### Alternative Libraries Evaluated

* `slowapi`: Rate limiting library designed for ASGI/FastAPI applications (useful if migrating web layers).
* `redis-py`: Direct sliding-window implementation using Redis sorted sets.

#### Pros, Cons, and Code Impact

* `Flask-Limiter` Pros: Declarative route decorators (e.g., `@limiter.limit("5 per minute")`), supports multiple storage backends (Redis, Memcached, in-memory), handles rate limit headers (`X-RateLimit-*`) automatically.
* `Flask-Limiter` Cons: Default memory storage is not shared across Gunicorn worker processes unless configured with Redis or Memcached.
* Code Impact: Removes custom database schema tables (`login_throttle`), sliding exponential backoff math, and ~650 lines of custom lockout guard code.

### 6. In-Memory and Disk Caching

#### Current Codebase State

* `src/html/lib/pagecache.py` (183 lines): Custom thread-safe in-memory cache for public API HTTP responses with TTL eviction.
* `src/html/lib/tilecache.py` (516 lines): Custom file-backed disk cache for 3D Galaxy Map JSON tiles with manual directory cleanup, stale file eviction, and database stamp checking.

#### Recommended Stack

* `cachetools`: In-memory LRU, LFU, and TTL caching algorithms.
* `diskcache`: SQLite and file-backed disk cache optimized for multi-process Python applications.

#### Alternative Libraries Evaluated

* `dogpile.cache`: Comprehensive caching API supporting Redis, Memcached, and file backends with stampede protection (lock stealing).
* `Flask-Caching`: High-level wrapper combining in-memory and disk caching for Flask routes.

#### Pros, Cons, and Code Impact

* `cachetools` Pros: Thread-safe, lightweight, declarative function decorators (`@ttl_cache`, `@lru_cache`).
* `diskcache` Pros: POSIX file-locking, thread-safe and process-safe concurrent access, fast SQLite index metadata lookups, auto-managed cache size limits.
* `dogpile.cache` Pros: Excellent protection against cache stampedes on expensive generation routes.
* Code Impact: Deletes custom disk traversal loops, tempfile renames, atomic file replacement routines, and ~700 lines of cache management boilerplate.

### 7. Markdown Rendering and Data Validation

#### Current Codebase State

* `src/html/lib/mdconvert.py` (179 lines): Custom Markdown-to-HTML converter using regex string parsing to handle ATX headings, GFM pipe tables, and scientific superscript tags (`<sup>`).
* `src/planetgen/generation/validation.py` (783 lines): Imperative checking functions validating planet parameters, moon system orbits, and physical boundaries.

#### Recommended Stack

* `markdown` (or `mistune`): Standard Python Markdown parsers.
* `Pydantic`: High-performance data validation using Python type hints.

#### Alternative Libraries Evaluated

* `mistletoe`: Fast, extensible pure-Python CommonMark parser.
* `marshmallow`: Ecosystem for object serialization and payload validation.

#### Pros, Cons, and Code Impact

* `markdown` Pros: Full extension support (tables, superscripts, sanitize HTML), fully compliant with CommonMark/GFM specifications.
* `mistune` Pros: Extremely fast C-compiled parser option, highly customizable renderers.
* `Pydantic` Pros: Rust-backed validation engine (`pydantic-core`), automatic type coercion, JSON schema generation, declarative dataclass definitions.
* Code Impact: Eliminates custom regex HTML escaping/un-escaping rules in `mdconvert.py` and removes ~780 lines of imperative validation functions in `validation.py`.

## Implementation Matrix and Phased Adoption Roadmap

```
Phase 1: Low Risk / High Reduction (Security & Utility Modules)
├── Replace src/planetgen/admin/totp.py       --> pyotp
├── Replace src/planetgen/admin/qrcode.py  --> qrcode / segno
└── Replace src/html/lib/mdconvert.py        --> markdown

Phase 2: Numerical & Cache Refactoring (Performance & Physical Accuracy)
├── Refactor src/planetgen/physics/kepler.py --> scipy.optimize
├── Refactor physical_constants.py              --> astropy.constants & astropy.units
├── Replace src/html/lib/pagecache.py           --> cachetools
└── Replace src/html/lib/tilecache.py           --> diskcache

Phase 3: Core Architecture Refactoring (Persistence & Task Execution)
├── Refactor src/planetgen/db/store.py          --> SQLAlchemy Core / ORM
├── Replace planetgen.cli.migrate                    --> Alembic
├── Replace src/stellarObjects/workQueue.py     --> RQ or Huey
└── Refactor src/planetgen/generation/validation.py  --> Pydantic
```

## Reference List

Astropy Collaboration, Price-Whelan, A. M., Lim, P. L., Earl, N., Starkman, N., Bradley, L., Shupe, D. L., Patil, A. A., Corrales, L., Brasseur, C. E., Nöthe, M., Donath, A., Tollerud, E. J., Morris, B. M., Ginsburg, A., Vaher, E., Weaver, B. A., Tocknell, J., Jamieson, W., ... Sipőcz, B. M. (2022). The Astropy Project: Sustaining and growing a community-oriented open-source project and the latest major release (v5.0) of the core package. The Astrophysical Journal, 935(2), Article 167. https://doi.org/10.3847/1538-4357/ac7c74

Bayer, M. (2012). SQLAlchemy. In A. Brown & G. Wilson (Eds.), The Architecture of Open Source Applications (Vol. 2, pp. 1–20). aosabook.org. http://aosabook.org/en/sqlalchemy.html

Bayer, M., & SQLAlchemy Contributors. (2026). SQLAlchemy: The Python SQL toolkit and Object Relational Mapper (Version 2.0) [Computer software]. https://www.sqlalchemy.org/

Colvin, S., & Pydantic Contributors. (2023). Pydantic: Data validation using Python type hints (Version 2.0) [Computer software]. https://pydantic.dev/

Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R., Virtanen, P., Cournapeau, D., Wieser, E., Taylor, J., Berg, S., Smith, N. J., Kern, R., Picus, M., Hoyer, S., van Kerkwijk, M. H., Brett, M., Haldane, A., del Río, F. A., Wiebe, M., Peterson, P., ... Oliphant, T. E. (2020). Array programming with NumPy. Nature, 585(7825), 357–362. https://doi.org/10.1038/s41586-020-2649-2

Jenks, G. (2021). DiskCache: Disk and file backed persistent cache in Python (Version 5.2) [Computer software]. https://grantjenks.com/docs/diskcache/

Lissner, A., & PyOTP Maintainers. (2023). PyOTP: Python One-Time Password Library (Version 2.9.0) [Computer software]. Python Software Foundation. https://github.com/pyauth/pyotp

MacFarlane, J. (2019). CommonMark specification (Version 0.29). CommonMark. https://spec.commonmark.org/0.29/

M'Raihi, D., Machani, S., Pei, M., & Rydell, J. (2011). TOTP: Time-Based One-Time Password Algorithm (RFC No. 6238). Internet Engineering Task Force. https://doi.org/10.17487/RFC6238

Murray, C. D., & Dermott, S. F. (1999). Solar system dynamics. Cambridge University Press. https://doi.org/10.1017/CBO9781139174817

Nayuki Project. (2022). QR Code generator library [Computer software]. https://www.nayuki.io/page/qr-code-generator-library

Python Markdown Developers. (2023). Python-Markdown (Version 3.5) [Computer software]. Python Software Foundation. https://python-markdown.github.io/

Roy, A. E. (2004). Orbital motion (4th ed.). Institute of Physics Publishing. https://doi.org/10.1201/9781420056884

RQ Contributors. (2024). RQ: Simple job queues for Python (Version 1.16) [Computer software]. GitHub. https://python-rq.org/

Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., Burovski, E., Peterson, P., Weckesser, W., Bright, J., van der Walt, S. J., Brett, M., Wilson, J., Millman, K. J., Mayorov, N., Nelson, A. R. J., Jones, E., Kern, R., Larson, E., ... SciPy 1.0 Contributors. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. Nature Methods, 17(3), 261–272. https://doi.org/10.1038/s41592-019-0686-2
