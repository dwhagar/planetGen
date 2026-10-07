# Single-Server Migration Guide & Codebase Modernization Outline

This guide outlines a step-by-step conversion strategy to modernize the planetGen architecture while maintaining a single-server, broker-less deployment model (Celery Project, 2024; Python Software Foundation, 2024). It replaces custom queues and DOM scripts with standard Python concurrency mechanisms and open-source frontend libraries.

## Phased Migration Progression

### Phase 1: Native In-Process Queue & Admin Management Refactoring

* **Target Files:** `src/jobRunner.py`, `src/stellarObjects/workQueue.py`, `src/html/web/admin_pages.py`, `src/html/web/templates/admin_queue.html`.
* **Objective:** Replace custom SQLite process-polling loops with Python's built-in `concurrent.futures.ProcessPoolExecutor` managed in-process, while expanding the admin web interface for queue inspection and control (Python Software Foundation, 2024).
* **Steps:**
  1. Wrap `ProcessPoolExecutor` inside a global task manager singleton within the Flask application process (Python Software Foundation, 2024).
  2. Refactor `workQueue.py` to maintain SQLite job state for persistence while delegating execution futures to the pool.
  3. Expand `admin_pages.py` and `admin_queue.html` to expose active task futures, worker utilization, job cancellations, and queue re-prioritization directly in the web UI.

### Phase 2: Real-Time SSE Log Streaming & Terminal Integration

* **Target Files:** `src/html/static/generatejobs.js`, `src/html/web/jobs.py`, `src/planetgen/util/log.py`.
* **Objective:** Transition process log rendering from client HTTP polling to Server-Sent Events (SSE) and Xterm.js (Xterm.js, 2024).
* **Steps:**
  1. Modify `log.py` to publish log outputs into thread-safe in-memory queues per job ID.
  2. Add an SSE route in `jobs.py` that streams log buffers to the client as text streams.
  3. Replace the `<div>` log containers in `generatejobs.js` with an Xterm.js terminal instance that listens to the SSE stream (Xterm.js, 2024).

### Phase 3: Data Grids & Controls Modernization

* **Target Files:** `src/html/lib/tabledisplay.py`, `src/html/static/sectormap.js`, Jinja templates in `src/html/web/templates/`.
* **Objective:** Eliminate server-rendered table DOM bloat and manual DOM construction in JavaScript.
* **Steps:**
  1. Refactor `tabledisplay.py` and Jinja templates to serve lightweight container shell markup rather than large pre-rendered HTML tables.
  2. Integrate TanStack Table and TanStack Virtual in `sectormap.js` to render catalog lists headlessly with viewport virtualization (TanStack, 2024).
  3. Standardize form controls, dropdowns, and modal dialogs using Shoelace web components (Shoelace, 2024).

### Phase 4: Spatial Partitioning & 3D Viewport Optimization

* **Target Files:** `src/html/static/galaxymap3d.js`, `src/html/static/galaxyprisms.js`, `src/planetgen/galaxy/geometry.py`.
* **Objective:** Eliminate GPU coordinate jitter and memory limits across astronomical scale jumps.
* **Steps:**
  1. Update `galaxyGeometry.py` to export sector bounding volumes into a Hierarchical Bounding Volume Hierarchy (BVH) structure (Open Geospatial Consortium, 2023).
  2. Refactor `galaxymap3d.js` to stream scene nodes dynamically based on camera distance using `three-mesh-bvh` or an OGC 3D Tiles renderer (Open Geospatial Consortium, 2023).
  3. Deprecate custom frustum culling and raw prism mathematics in `galaxyprisms.js`.

## Potential Code Conversion Pain Points

### SQLite Concurrency & Database Locking

* **Issue:** `ProcessPoolExecutor` worker processes writing job completions while the Flask web app reads/writes worker statuses can trigger SQLite database lock errors (Python Software Foundation, 2024).
* **Mitigation:** Enable Write-Ahead Logging (WAL) mode in SQLite during startup, shorten transaction lock durations in `workQueue.py`, and funnel status updates through a dedicated thread-safe queue in the master application process.

### WSGI Worker Starvation During Log Streaming

* **Issue:** Standard WSGI servers (e.g., Gunicorn with synchronous workers) block an entire HTTP worker process for each open SSE log stream connection.
* **Mitigation:** Run the single server instance using an asynchronous or eventlet/gevent worker class, or use thread-pooled streaming responses to prevent worker pool depletion.

### Asset Bundling in a Vanilla JS / Jinja Pipeline

* **Issue:** Introducing modern npm modules (TanStack, Xterm.js, Shoelace) into a codebase built around Jinja templates and static vanilla JS files (`src/html/static/`) can break simple static file serving.
* **Mitigation:** Use standalone ES Module (ESM) builds or integrate a lightweight bundler (such as Vite or ESBuild) configured to output compiled ESM bundles directly into `src/html/static/`.

### Coordinate Frame Shift for 3D Rendering

* **Issue:** Translating absolute galactic coordinates from `galaxyGeometry.py` directly to WebGL causes single-precision float32 vertex jitter in Three.js (Open Geospatial Consortium, 2023).
* **Mitigation:** Compute visual coordinates relative to local sector origins (Camera-Relative Rendering) before passing geometry buffers to the GPU.

## Integration & Implementation Guide

### How should single-server queue management be structured without Redis or Celery?

* **Architecture:** Use a hybrid model consisting of SQLite for persistent state and Python's `concurrent.futures.ProcessPoolExecutor` for execution (Python Software Foundation, 2024).
* **Why:** SQLite preserves task history across server restarts, while `ProcessPoolExecutor` utilizes multi-core CPU performance without external broker dependencies (Celery Project, 2024; Python Software Foundation, 2024).
* **Admin Web Control:** Flask admin endpoints in `admin_pages.py` interact with the executor's `Future` objects. Canceling a job calls `future.cancel()` or terminates the worker process, while pausing queues holds new dispatches in SQLite without worker pool execution.

### How can SSE log streaming be implemented efficiently on a single Python backend?

* **Architecture:** Implement an in-memory pub-sub channel using Python's `queue.Queue` per active job ID in `log.py`.
* **Why:** Avoids polling SQLite for log lines, reducing disk I/O and server load.
* **Client Handshake:** `generatejobs.js` initiates an `EventSource('/api/jobs/<id>/stream')` connection. The server yields log chunks as text streams, which Xterm.js appends directly into its terminal canvas buffer (Xterm.js, 2024).

### How should modern open-source UX components be integrated into Jinja2 templates?

* **Architecture:** Adopt a progressive web component model using Shoelace and ESM-imported JavaScript libraries (Shoelace, 2024).
* **Why:** Preserves the existing Jinja layout architecture in `src/html/web/templates/` while replacing client-side rendering with fast, accessible web components (Shoelace, 2024).
* **Data Flow:** Jinja renders lightweight HTML structural shells containing `data-*` attributes; static scripts in `src/html/static/` hydrate these containers using TanStack Table and Xterm.js (TanStack, 2024; Xterm.js, 2024).

## References

Celery Project. (2024). *Celery: Distributed task queue*. https://docs.celeryq.dev/

Open Geospatial Consortium. (2023). *3D Tiles specification (Version 1.1)*. https://www.ogc.org/standard/3dtiles/

Python Software Foundation. (2024). *concurrent.futures — Managing pools of concurrent tasks*. https://docs.python.org/3/library/concurrent.futures.html

Shoelace. (2024). *Web components for building modern web applications*. https://shoelace.style/

TanStack. (2024). *TanStack table: Headless UI for building powerful tables & datagrids*. https://tanstack.com/table/

Xterm.js. (2024). *Xterm.js: A terminal for the web*. https://xtermjs.org/
