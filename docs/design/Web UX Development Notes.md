# Architectural Improvements & Simplification via Open-Source UX Libraries

This proposal evaluates refactoring custom infrastructure within the `planetGen` codebase using production-ready open-source libraries across five distinct subsystems.

## 1. General Controls & Paginated Lists

Building custom paginated lists and interactive UI controls in `tabledisplay.py` and `sectormap.js` creates unnecessary maintenance overhead.

### Proposed Architecture

* **Data Grids & Virtualization:** Integrate TanStack Table for headless pagination, sorting, and filtering, paired with TanStack Virtual (or `lit-virtualizer`) for viewport rendering (TanStack, 2024).
* **Control Primitives:** Replace custom web controls with web-component suites like Shoelace / Web Awesome or headless primitive suites such as Radix Primitives or React Aria (Shoelace, 2024).

### Architectural Evaluation

* **Pros:**
  * **Codebase:** Removes hundreds of lines of fragile manual DOM manipulation from `sectormap.js` and Jinja templates.
  * **Server Performance:** Reduces payload sizes by serializing raw JSON arrays rather than rendering heavy HTML tables on the backend.
  * **Client Performance:** Virtualization limits DOM nodes to only visible viewport elements, preventing browser layout reflows and memory leaks when browsing catalogs of 1,000+ celestial bodies.
* **Cons:**
  * **Codebase & Client:** Requires introducing a JavaScript build or bundling step into an existing vanilla JavaScript pipeline, slightly increasing initial client bundle size.

## 2. The 3D Multi-Scale "Galactic-to-System" Viewport

`galaxymap3d.js` and `galaxyprisms.js` encounter hardware single-precision float32 precision limits (z-fighting and vertex jitter) when transitioning across astronomical distances.

### Proposed Architecture

* **3D Tiles Specification:** Adopt the Open Geospatial Consortium (OGC) 3D Tiles specification using libraries such as `3d-tiles-renderer` or `three-mesh-bvh` (Open Geospatial Consortium, 2023).
* **Spatial Partitioning (Octree / BVH):** Partition galaxy sectors and stellar systems into Hierarchical Bounding Volume Hierarchies (BVH) to dynamically stream Level-of-Detail (LOD) nodes based on camera proximity (Open Geospatial Consortium, 2023).

### Architectural Evaluation

* **Pros:**
  * **Codebase:** Offloads manual camera calculations and frustum culling algorithms from custom scripts into standard 3D geospatial specifications.
  * **Server Performance:** Enables streaming spatial bounding chunks on demand rather than serializing multi-megabyte galactic coordinates in a single database payload.
  * **Client Performance:** Resolves GPU floating-point precision jitter by localizing coordinate spaces per spatial node; drastically lowers GPU VRAM usage via dynamic LOD streaming.
* **Cons:**
  * **Codebase:** High refactoring complexity to translate `planetGen` data structures into hierarchical 3D tilesets.

## 3. Log Streaming & Terminal Windows

Piping Python stdout/stderr process output directly into standard browser `<div>` tags in `generatejobs.js` degrades browser performance as log buffers grow.

### Proposed Architecture

* **Web Terminal Emulator:** Replace DOM text containers with Xterm.js for fast, memory-capped log rendering with native ANSI escape code support (Xterm.js, 2024).
* **Communication Layer:** Replace HTTP endpoint polling with Server-Sent Events (SSE) or WebSockets to push output streams from the backend directly to the client (Xterm.js, 2024).

### Architectural Evaluation

* **Pros:**
  * **Codebase:** Eliminates manual text-parsing, string-concatenation, and auto-scrolling logic in `generatejobs.js`.
  * **Server Performance:** SSE eliminates high-frequency HTTP polling requests, freeing server workers and reducing SQLite query overhead.
  * **Client Performance:** Xterm.js utilizes WebGL/Canvas virtual scrolling with configurable buffer limits, ensuring steady 60 FPS rendering regardless of log throughput.
* **Cons:**
  * **Server Performance & Codebase:** Requires an asynchronous server layer (such as Gevent or ASGI) to handle persistent SSE connections and proxy configuration adjustments.

## 4. Progress Bars and Status Feedback

Custom DOM progress updates in `generatefolds.js` and tracking routines in `progressRate.py` increase UI complexity.

### Proposed Architecture

* **Native Elements & Accessibility:** Standardize progress tracking using native HTML `<progress>` elements enhanced with W3C ARIA live region attributes (`role="progressbar"`, `aria-valuenow`).

### Architectural Evaluation

* **Pros:**
  * **Codebase:** Zero external library dependencies; allows complete deletion of bespoke DOM-updating routines in `generatefolds.js`.
  * **Server Performance:** Negligible overhead; server streams small integer status updates.
  * **Client Performance:** Browser executes native hardware-accelerated animations and updates screen readers with zero JavaScript main-thread evaluation overhead.
* **Cons:**
  * **Codebase:** Styling cross-browser default progress elements requires modern CSS custom properties or pseudo-element overrides.

## 5. Background Task Scheduling & Parallel Processing

The custom task queue in `workQueue.py` and `planetgen.cli.job` uses SQLite locking mechanisms that can lead to race conditions, worker process crashes, and incomplete state recovery.

### Proposed Architecture

* **Distributed Task Queues:** Adopt standard Python job schedulers such as Celery, RQ (Redis Queue), or Dramatiq to manage multi-worker execution, serialization, and fault recovery (Celery Project, 2024).
* **Single-Server Alternative:** Utilize Python's built-in `concurrent.futures.ProcessPoolExecutor` paired with lightweight dispatchers like APScheduler or FastAPI BackgroundTasks for broker-less setups (Celery Project, 2024).

### Architectural Evaluation

* **Pros:**
  * **Codebase:** Fully replaces `workQueue.py` and `planetgen.cli.job`, eliminating custom process management and recovery logic.
  * **Server Performance:** Resolves SQLite `database is locked` errors caused by concurrent worker writes, enabling true multi-core CPU scaling for procedural generation algorithms.
  * **Client Performance:** Faster background task execution reduces user wait times and eliminates UI job-status timeouts.
* **Cons:**
  * **Server Performance & Codebase:** Distributed brokers (e.g., Redis or RabbitMQ) introduce an extra system dependency that must be deployed and monitored.

## References

Celery Project. (2024). *Celery: Distributed task queue*. https://docs.celeryq.dev/

Open Geospatial Consortium. (2023). *3D Tiles specification (Version 1.1)*. https://www.ogc.org/standard/3dtiles/

Shoelace. (2024). *Web components for building modern web applications*. https://shoelace.style/

TanStack. (2024). *TanStack table: Headless UI for building powerful tables & datagrids*. https://tanstack.com/table/

Xterm.js. (2024). *Xterm.js: A terminal for the web*. https://xtermjs.org/
