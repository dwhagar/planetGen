# Progress ledger: the phase each finished item belonged to

Boss (2026-10-16, 11:09 AM EDT) asked for a progress bar on every tier of the
tier plan. A phase file drops an item's row when the item ships, so the phase of a
finished item has to be kept somewhere; this file keeps it. The rows were
read back from the git history of `docs/plan/phase-*.md` on 2026-10-16 (the phase an
item last sat in), and every later retirement adds its row in the same PR
(the retire step copies the row from the phase file before deleting it).

`scripts/build_todo_docs.py` counts, per phase, the items here whose
number-map status is done (finished), the items still open in the phase files,
and the dropped ones (closed, superseded, folded or not reproduced; they are
not counted either way). Progress is finished / (finished + open). A finished item
with no row here was built before the phase plan began on 2026-10-01 and is
counted under "built before the plan".

The build fails if an item is done in the number map with a PR above the
baseline below and has no row here, so a retire PR cannot forget it.

Baseline: PR #1817.

| ID | Phase | Last thread |
|---|---|---|
| ADM.14 | 0 | Groundwork: web components |
| ADM.15 | 1 | Queue |
| ADM.16 | 0 | Bugfixes: prevalence |
| ADM.17 | 2 | Reproducible galaxies |
| ADM.18 | 1 | Reproducible galaxies |
| ADM.19 | 2 | Daily maintenance |
| ADM.20 | 3 | Daily maintenance |
| ADM.21 | 0 | Groundwork: data model |
| ADM.22 | 0 | Groundwork: queue, logs and caches |
| ADM.23 | 0 | Groundwork: queue, logs and caches |
| ADM.24 | 0 | Groundwork: queue, logs and caches |
| ADM.25 | 0 | Groundwork: queue, logs and caches |
| ADM.26 | 0 | Groundwork: queue, logs and caches |
| ADM.27 | 0 | Bugfixes: generation |
| ADM.28 | 1 | Generate page |
| ADM.29 | 1 | Generate page |
| ADM.30 | 1 | Generate page |
| ADM.31 | 1 | Generate page |
| ADM.33 | 0 | Bugfixes: console and progress |
| ADM.34 | 0 | Groundwork: web components |
| ADM.37 | 0 | Bugfixes: prevalence |
| ADM.38 | 0 | Bugfixes: console and progress |
| ADM.39 | 1 | Bugs from the GitHub issues |
| ADM.40 | 1 | Bugs from the GitHub issues |
| ADM.41 | 1 | Bugs from the GitHub issues |
| ADM.42 | 1 | Foundations for the issue features |
| ADM.45 | 1 | Generate page |
| ADM.46 | 1 | Generate page |
| ADM.47 | 1 | Generate page |
| ADM.48 | 1 | Foundations for the issue features |
| ADM.49 | 1 | Generate page |
| ADM.51 | 1 | Bugs from the GitHub issues |
| ADM.53 | 1 | Bugs from the GitHub issues |
| ADM.54 | 1 | Foundations for the issue features |
| API.9 | 1 | API |
| API.16 | 2 | API |
| API.20 | 2 | Recipes |
| API.21 | 2 | Recipes |
| API.22 | 1 | API |
| API.23 | 1 | API |
| API.25 | 1 | Maps |
| DB.2 | 0 | Database fixes |
| DB.3 | 0 | Database fixes |
| DB.4 | 0 | Database fixes |
| DB.5 | 0 | Database fixes |
| DB.6 | 0 | Parallel path (top priority) |
| DB.7 | 1 | Reproducible galaxies |
| DB.8 | 1 | Database consistency check |
| DB.10 | 3 | Database consistency check |
| DB.11 | 0 | Groundwork: data model |
| DB.12 | 0 | Bugfixes: CI red |
| DB.13 | 0 | Groundwork: data model |
| DB.14 | 1 | Maps |
| DB.15 | 1 | Foundations for the issue features |
| DB.19 | 1 | Foundations for the issue features |
| DB.20 | 1 | Database consistency check |
| DB.21 | 1 | Database consistency check |
| DB.22 | 1 | API |
| DB.23 | 1 | Foundations for the issue features |
| DB.24 | 1 | Maps |
| DB.25 | 1 | Bugs from the GitHub issues |
| DB.26 | 1 | Foundations for the issue features |
| DB.27 | 2 | Features from the GitHub issues |
| DB.28 | 1 | Tech levels and facilities |
| DB.32 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| DB.34 | 2 | Picker |
| DB.35 | 2 | Picker |
| DOC.5 | 1 | Database consistency check |
| DOC.17 | 1 | API |
| DOC.21 | 1 | Bugs from the GitHub issues |
| DOC.22 | 1 | Bugs from the GitHub issues |
| DOC.23 | 1 | Bugs from the GitHub issues |
| DOC.25 | 1 | Bugs from the GitHub issues |
| DOC.26 | 1 | Bugs from the GitHub issues |
| DOC.27 | 1 | Bugs from the GitHub issues |
| DOC.28 | 1 | Bugs from the GitHub issues |
| DOC.29 | 1 | Bugs from the GitHub issues |
| DOC.30 | 1 | Bugs from the GitHub issues |
| DOC.31 | 1 | Bugs from the GitHub issues |
| DOC.32 | 1 | Bugs from the GitHub issues |
| DOC.33 | 1 | Bugs from the GitHub issues |
| DOC.34 | 1 | Bugs from the GitHub issues |
| DOC.35 | 1 | Bugs from the GitHub issues |
| DOC.36 | 1 | Bugs from the GitHub issues |
| DOC.37 | 1 | Bugs from the GitHub issues |
| DOC.38 | 1 | Bugs from the GitHub issues |
| DOC.39 | 1 | Bugs from the GitHub issues |
| DOC.40 | 1 | Bugs from the GitHub issues |
| DOC.41 | 1 | Bugs from the GitHub issues |
| DOC.42 | 1 | Bugs from the GitHub issues |
| DOC.43 | 1 | Bugs from the GitHub issues |
| DOC.44 | 1 | Bugs from the GitHub issues |
| DOC.46 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.25 | 0 | Physics bugs |
| GEN.31 | 0 | Galaxy geometry and names |
| GEN.32 | 0 | Parallel path (top priority) |
| GEN.34 | 0 | Physics bugs |
| GEN.35 | 0 | Physics bugs |
| GEN.36 | 0 | Physics bugs |
| GEN.37 | 0 | Physics bugs |
| GEN.38 | 0 | Generation bugs |
| GEN.39 | 0 | Parallel path (top priority) |
| GEN.41 | 1 | Galaxy gen |
| GEN.43 | 2 | Backfill |
| GEN.44 | 0 | Sector stats and colors |
| GEN.45 | 0 | Physics bugs |
| GEN.46 | 0 | Galaxy geometry and names |
| GEN.47 | 0 | Generation bugs |
| GEN.48 | 0 | Bugfixes: prevalence |
| GEN.49 | 0 | Binary pairs and single-system forcing |
| GEN.50 | 0 | Binary pairs and single-system forcing |
| GEN.51 | 0 | Binary pairs and single-system forcing |
| GEN.52 | 0 | Bugfixes: prevalence |
| GEN.53 | 0 | Binary pairs and single-system forcing |
| GEN.54 | 0 | Binary pairs and single-system forcing |
| GEN.56 | 1 | Reproducible galaxies |
| GEN.58 | 1 | Reproducible galaxies |
| GEN.59 | 1 | Reproducible galaxies |
| GEN.60 | 0 | Generation bugs |
| GEN.61 | 2 | Daily maintenance |
| GEN.62 | 0 | Binary pairs and single-system forcing |
| GEN.63 | 1 | Reproducible galaxies |
| GEN.64 | 1 | Interstellar IDs |
| GEN.65 | 0 | Bugfixes: generation |
| GEN.66 | 0 | Groundwork: data model |
| GEN.67 | 0 | Groundwork: names from IDs |
| GEN.68 | 0 | Groundwork: names from IDs |
| GEN.69 | 0 | Groundwork: names from IDs |
| GEN.70 | 0 | Groundwork: names from IDs |
| GEN.71 | 0 | Groundwork: names from IDs |
| GEN.72 | 0 | Groundwork: names from IDs |
| GEN.73 | 0 | Groundwork: names from IDs |
| GEN.74 | 0 | Groundwork: data model |
| GEN.75 | 0 | Groundwork: nebula shapes |
| GEN.76 | 0 | Bugfixes: generation |
| GEN.77 | 0 | Bugfixes: generation |
| GEN.78 | 0 | Bugfixes: generation |
| GEN.79 | 0 | Bugfixes: generation |
| GEN.80 | 0 | Bugfixes: generation |
| GEN.81 | 0 | Bugfixes: console and progress |
| GEN.82 | 0 | Bugfixes: generation |
| GEN.83 | 1 | Habitability |
| GEN.84 | 1 | Habitability |
| GEN.85 | 1 | Habitability |
| GEN.86 | 1 | Habitability |
| GEN.87 | 1 | Habitability |
| GEN.88 | 1 | Habitability |
| GEN.89 | 1 | Habitability |
| GEN.94 | 1 | Nebula planets |
| GEN.96 | 1 | Generate page |
| GEN.97 | 1 | Generate page |
| GEN.98 | 1 | Galaxy gen |
| GEN.100 | 1 | Galaxy gen |
| GEN.101 | 1 | Galaxy gen |
| GEN.104 | 1 | Orbital updates |
| GEN.106 | 1 | Orbital updates |
| GEN.107 | 1 | Orbital updates |
| GEN.108 | 1 | Orbital updates |
| GEN.113 | 2 | Anomalies |
| GEN.114 | 3 | Anomalies |
| GEN.116 | 0 | Bugfixes: generation |
| GEN.117 | 0 | Groundwork: map engine |
| GEN.118 | 0 | Bugfixes: generation |
| GEN.119 | 0 | Bugfixes: generation |
| GEN.120 | 1 | Galaxy gen |
| GEN.121 | 1 | Orbital updates |
| GEN.122 | 1 | Orbital updates |
| GEN.123 | 1 | Orbital updates |
| GEN.124 | 1 | Orbital updates |
| GEN.125 | 1 | Orbital updates |
| GEN.126 | 0 | Groundwork: data model |
| GEN.127 | 1 | Bugs from the GitHub issues |
| GEN.128 | 1 | Foundations for the issue features |
| GEN.131 | 1 | Foundations for the issue features |
| GEN.132 | 1 | Foundations for the issue features |
| GEN.133 | 1 | Foundations for the issue features |
| GEN.137 | 1 | Routing |
| GEN.138 | 2 | Orbital updates |
| GEN.147 | 2 | Classes |
| GEN.148 | 1 | Habitability |
| GEN.149 | 2 | Nebula planets |
| GEN.150 | 1 | Galaxy gen |
| GEN.151 | 2 | Nebula planets |
| GEN.152 | 1 | Foundations for the issue features |
| GEN.153 | 2 | Anomalies |
| GEN.154 | 2 | Anomalies |
| GEN.155 | 1 | Galaxy gen |
| GEN.156 | 1 | Reproducible galaxies |
| GEN.158 | 3 | Anomalies |
| GEN.159 | 3 | Anomalies |
| GEN.160 | 3 | Anomalies |
| GEN.165 | 2 | Features from the GitHub issues |
| GEN.166 | 1 | Foundations for the issue features |
| GEN.167 | 1 | Foundations for the issue features |
| GEN.168 | 1 | Foundations for the issue features |
| GEN.169 | 1 | Foundations for the issue features |
| GEN.170 | 1 | Database consistency check |
| GEN.171 | 1 | Database consistency check |
| GEN.172 | 1 | Database consistency check |
| GEN.173 | 1 | Database consistency check |
| GEN.174 | 1 | Database consistency check |
| GEN.175 | 1 | Database consistency check |
| GEN.176 | 1 | Database consistency check |
| GEN.177 | 1 | Habitability |
| GEN.178 | 1 | Habitability |
| GEN.182 | 1 | Bugs from the GitHub issues |
| GEN.183 | 1 | Foundations for the issue features |
| GEN.184 | 1 | Foundations for the issue features |
| GEN.185 | 1 | Foundations for the issue features |
| GEN.186 | 1 | Generate page |
| GEN.187 | 1 | Foundations for the issue features |
| GEN.188 | 1 | Foundations for the issue features |
| GEN.189 | 2 | Maps |
| GEN.190 | 1 | Generate page |
| GEN.191 | 1 | Generate page |
| GEN.192 | 1 | Foundations for the issue features |
| GEN.194 | 1 | Maps |
| GEN.195 | 1 | Maps |
| GEN.196 | 1 | Maps |
| GEN.200 | 1 | Foundations for the issue features |
| GEN.201 | 1 | Foundations for the issue features |
| GEN.203 | 1 | Foundations for the issue features |
| GEN.204 | 1 | Foundations for the issue features |
| GEN.205 | 1 | Foundations for the issue features |
| GEN.206 | 1 | Bugs from the GitHub issues |
| GEN.207 | 2 | Maps |
| GEN.208 | 2 | Maps |
| GEN.211 | 1 | Bugs from the GitHub issues |
| GEN.213 | 1 | Habitability |
| GEN.214 | 1 | Maps |
| GEN.215 | 1 | Maps |
| GEN.216 | 1 | Maps |
| GEN.217 | 1 | Maps |
| GEN.222 | 1 | Maps |
| GEN.223 | 1 | Maps |
| GEN.224 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.225 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.226 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.227 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.228 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.230 | 1 | Maps |
| GEN.232 | 1 | Maps |
| GEN.233 | 2 | Orbital updates |
| GEN.234 | 2 | Orbital updates |
| GEN.235 | 2 | Orbital updates |
| GEN.250 | 2 | Orbital updates |
| GEN.251 | 2 | Orbital updates |
| GEN.254 | 2 | Orbital updates |
| GEN.260 | 2 | Orbital updates |
| GEN.261 | 2 | Orbital updates |
| GEN.262 | 2 | Orbital updates |
| GEN.265 | 2 | Orbital updates |
| GEN.266 | 2 | Orbital updates |
| GEN.268 | 2 | Orbital updates |
| GEN.277 | 2 | Anomalies |
| GEN.278 | 2 | Anomalies |
| GEN.279 | 2 | Anomalies |
| GEN.288 | 1 | Maps |
| MAP.52 | 0 | Galaxy map picker and arc |
| MAP.53 | 0 | Galaxy Map drill-down |
| MAP.54 | 0 | Galaxy Map drill-down |
| MAP.55 | 0 | Galaxy map picker and arc |
| MAP.56 | 0 | Galaxy Map drill-down |
| MAP.57 | 0 | System Map |
| MAP.58 | 2 | Galaxy Map |
| MAP.60 | 0 | Galaxy map picker and arc |
| MAP.61 | 0 | Groundwork: map engine |
| MAP.62 | 3 | 3D system |
| MAP.63 | 0 | Map groundwork |
| MAP.64 | 0 | Map groundwork |
| MAP.65 | 0 | Groundwork: map engine |
| MAP.66 | 0 | Groundwork: map engine |
| MAP.67 | 0 | Groundwork: map engine |
| MAP.68 | 0 | Groundwork: map engine |
| MAP.69 | 2 | 3D system |
| MAP.70 | 2 | 3D system |
| MAP.71 | 3 | 3D system |
| MAP.72 | 3 | 3D system |
| MAP.73 | 3 | 3D system |
| MAP.74 | 3 | 3D system |
| MAP.76 | 0 | Galaxy Map drill-down |
| MAP.77 | 0 | Galaxy Map drill-down |
| MAP.78 | 0 | Galaxy Map drill-down |
| MAP.79 | 0 | Groundwork: map engine |
| MAP.80 | 0 | Sector stats and colors |
| MAP.81 | 0 | Map groundwork |
| MAP.82 | 0 | Map groundwork |
| MAP.83 | 0 | Map groundwork |
| MAP.84 | 0 | Map groundwork |
| MAP.85 | 0 | Galaxy map picker and arc |
| MAP.86 | 0 | Sector stats and colors |
| MAP.87 | 0 | Map groundwork |
| MAP.88 | 0 | System Map |
| MAP.89 | 1 | System Map |
| MAP.90 | 0 | Galaxy geometry and names |
| MAP.91 | 0 | Galaxy Map follow-ups |
| MAP.92 | 0 | System Map |
| MAP.93 | 0 | Galaxy Map follow-ups |
| MAP.94 | 0 | Galaxy Map follow-ups |
| MAP.95 | 0 | Groundwork: map engine |
| MAP.96 | 0 | Galaxy Map drill-down |
| MAP.97 | 0 | Galaxy Map drill-down |
| MAP.98 | 0 | Galaxy Map drill-down |
| MAP.99 | 0 | Galaxy Map drill-down |
| MAP.100 | 0 | Galaxy Map drill-down |
| MAP.101 | 0 | Galaxy Map drill-down |
| MAP.102 | 0 | Groundwork: map engine |
| MAP.103 | 0 | Groundwork: nebula shapes |
| MAP.104 | 0 | Groundwork: nebula shapes |
| MAP.105 | 0 | Groundwork: nebula shapes |
| MAP.106 | 0 | Groundwork: map engine |
| MAP.107 | 0 | Groundwork: map engine |
| MAP.108 | 0 | Groundwork: map engine |
| MAP.109 | 0 | Groundwork: map engine |
| MAP.110 | 0 | Groundwork: map engine |
| MAP.111 | 0 | Groundwork: map engine |
| MAP.112 | 0 | Groundwork: map engine |
| MAP.113 | 0 | Groundwork: map engine |
| MAP.114 | 0 | Bugfixes: CI red |
| MAP.115 | 0 | Bugfixes: maps and pages |
| MAP.116 | 0 | Groundwork: map engine |
| MAP.117 | 0 | Bugfixes: maps and pages |
| MAP.118 | 0 | Bugfixes: maps and pages |
| MAP.123 | 1 | Maps |
| MAP.124 | 1 | Maps |
| MAP.125 | 3 | Infinite zoom |
| MAP.126 | 3 | 3D system |
| MAP.127 | 0 | Groundwork: map engine |
| MAP.128 | 1 | Maps |
| MAP.129 | 1 | Maps |
| MAP.130 | 1 | Maps |
| MAP.131 | 2 | Galaxy Map |
| MAP.133 | 0 | Bugfixes: maps and pages |
| MAP.134 | 1 | Maps |
| MAP.135 | 1 | Bugs from the GitHub issues |
| MAP.136 | 1 | Bugs from the GitHub issues |
| MAP.137 | 1 | Bugs from the GitHub issues |
| MAP.138 | 1 | Bugs from the GitHub issues |
| MAP.144 | 2 | Galaxy Map |
| MAP.148 | 1 | Maps |
| MAP.149 | 1 | Maps |
| MAP.150 | 1 | Maps |
| MAP.151 | 1 | Maps |
| MAP.153 | 1 | Maps |
| MAP.154 | 1 | Maps |
| MAP.155 | 1 | Maps |
| MAP.157 | 1 | Maps |
| MAP.158 | 1 | Maps |
| MAP.159 | 1 | Maps |
| MAP.161 | 1 | Maps |
| MAP.162 | 1 | Maps |
| MAP.163 | 1 | Maps |
| MAP.164 | 1 | Maps |
| MAP.165 | 1 | Maps |
| MAP.166 | 1 | Maps |
| MAP.175 | 1 | Maps |
| MAP.182 | 2 | Orbital updates |
| MAP.184 | 1 | Bugs from the GitHub issues |
| MAP.185 | 1 | Bugs from the GitHub issues |
| MAP.186 | 1 | Bugs from the GitHub issues |
| MAP.187 | 1 | Bugs from the GitHub issues |
| MAP.188 | 1 | Bugs from the GitHub issues |
| MAP.196 | 1 | Bugs from the GitHub issues |
| MAP.197 | 1 | Maps |
| MAP.198 | 1 | Maps |
| MAP.199 | 1 | Maps |
| MAP.200 | 1 | Maps |
| MAP.201 | 1 | Maps |
| MAP.202 | 1 | Maps |
| MAP.203 | 1 | Maps |
| MAP.204 | 1 | Maps |
| MAP.205 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| MAP.206 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| MAP.208 | 2 | Orbital updates |
| MAP.211 | 2 | Orbital updates |
| NAV.3 | 1 | Picker |
| NAV.7 | 0 | Groundwork: map engine |
| NAV.8 | 1 | References |
| NAV.9 | 1 | References |
| NAV.10 | 1 | Routing |
| NAV.11 | 1 | Routing |
| NAV.12 | 1 | Routing |
| NAV.13 | 0 | Groundwork: map engine |
| NAV.14 | 0 | Groundwork: map engine |
| NAV.15 | 0 | Groundwork: map engine |
| NAV.16 | 1 | References |
| NAV.20 | 2 | Courses |
| NAV.24 | 2 | Courses |
| NAV.29 | 0 | Groundwork: map engine |
| NAV.30 | 0 | Map groundwork |
| NAV.31 | 0 | Galaxy Map follow-ups |
| NAV.32 | 0 | Groundwork: map engine |
| NAV.33 | 0 | Groundwork: map engine |
| NAV.34 | 0 | Routing groundwork |
| NAV.35 | 0 | Routing guarantee |
| NAV.37 | 2 | Courses |
| NAV.38 | 0 | Galaxy geometry and names |
| NAV.40 | 0 | Galaxy Map follow-ups |
| NAV.41 | 0 | Bugfixes: maps and pages |
| NAV.42 | 1 | Routing |
| NAV.43 | 1 | Nearby search |
| NAV.44 | 1 | Nearby search |
| NAV.46 | 0 | Groundwork: map engine |
| NAV.47 | 1 | Routing |
| NAV.48 | 1 | Routing |
| NAV.50 | 1 | References |
| NAV.52 | 1 | Routing |
| NAV.53 | 1 | Routing |
| OPS.6 | 0 | Bugfixes: ops and flakes |
| OPS.7 | 0 | Bugfixes: ops and flakes |
| OPS.8 | 1 | Reproducible galaxies |
| OPS.9 | 0 | Bugfixes: console and progress |
| OPS.10 | 0 | Parallel path (top priority) |
| OPS.11 | 1 | Reproducible galaxies |
| OPS.12 | 3+ | Reproducible galaxies |
| OPS.13 | 1 | Reproducible galaxies |
| OPS.14 | 1 | Reproducible galaxies |
| OPS.18 | 2 | Daily maintenance |
| OPS.19 | 0 | Groundwork: queue, logs and caches |
| OPS.20 | 0 | Groundwork: layout and libraries |
| OPS.21 | 0 | Groundwork: layout and libraries |
| OPS.22 | 0 | Groundwork: layout and libraries |
| OPS.23 | 0 | Groundwork: layout and libraries |
| OPS.24 | 0 | Groundwork: layout and libraries |
| OPS.25 | 0 | Bugfixes: CI red |
| OPS.26 | 0 | Bugfixes: ops and flakes |
| OPS.27 | 0 | Groundwork: queue, logs and caches |
| OPS.31 | 2 | Daily maintenance |
| OPS.32 | 2 | Daily maintenance |
| OPS.33 | 1 | Reproducible galaxies |
| OPS.34 | 2 | Daily maintenance |
| OPS.36 | 1 | Foundations for the issue features |
| OPS.38 | 1 | Bugs from the GitHub issues |
| OPS.39 | 1 | Generate page |
| OPS.40 | 1 | Bugs from the GitHub issues |
| OPS.41 | 1 | Bookmarks |
| OPS.42 | 1 | Bugs from the GitHub issues |
| OPS.43 | 1 | Bugs from the GitHub issues |
| OPS.44 | 1 | Bugs from the GitHub issues |
| OPS.45 | 1 | Bugs from the GitHub issues |
| OPS.46 | 1 | Foundations for the issue features |
| OPS.47 | 1 | Foundations for the issue features |
| OPS.48 | 1 | Foundations for the issue features |
| OPS.49 | 1 | Foundations for the issue features |
| OPS.50 | 1 | Foundations for the issue features |
| OPS.51 | 1 | Bookmarks |
| OPS.52 | 1 | Bookmarks |
| OPS.53 | 1 | Bookmarks |
| OPS.54 | 1 | Bookmarks |
| PERF.1 | 0 | Sector stats and colors |
| PERF.11 | 0 | Sector stats and colors |
| PERF.19 | 0 | Groundwork: queue, logs and caches |
| PERF.21 | 0 | Parallel path (top priority) |
| PERF.22 | 0 | Parallel path (top priority) |
| PERF.23 | 0 | Parallel path (top priority) |
| PERF.24 | 0 | Groundwork: queue, logs and caches |
| PERF.25 | 0 | Groundwork: queue, logs and caches |
| PERF.26 | 0 | Bugfixes: generation |
| PERF.27 | 0 | Bugfixes: CI red |
| PERF.28 | 0 | Bugfixes: console and progress |
| PERF.31 | 1 | Foundations for the issue features |
| PERF.32 | 1 | Foundations for the issue features |
| PERF.33 | 1 | Foundations for the issue features |
| PERF.34 | 1 | Bugs from the GitHub issues |
| PERF.37 | 2 | Features from the GitHub issues |
| PERF.38 | 1 | Foundations for the issue features |
| PERF.39 | 1 | Foundations for the issue features |
| PERF.40 | 1 | Foundations for the issue features |
| PERF.42 | 1 | Foundations for the issue features |
| PERF.43 | 1 | Foundations for the issue features |
| PERF.44 | 1 | Foundations for the issue features |
| PERF.45 | 1 | Foundations for the issue features |
| PERF.46 | 1 | Foundations for the issue features |
| PERF.49 | 1 | Foundations for the issue features |
| PERF.50 | 1 | Bugs from the GitHub issues |
| PERF.51 | 1 | Foundations for the issue features |
| PERF.52 | 1 | Bugs from the GitHub issues |
| PERF.53 | 1 | Bugs from the GitHub issues |
| PERF.54 | 1 | Foundations for the issue features |
| PERF.55 | 1 | Foundations for the issue features |
| PERF.56 | 1 | Maps |
| PERF.57 | 1 | Maps |
| PERF.58 | 1 | Maps |
| PERF.59 | 1 | Maps |
| PERF.60 | 1 | Maps |
| PERF.61 | 1 | Maps |
| PERF.62 | 1 | Foundations for the issue features |
| PERF.63 | 1 | Maps |
| PERF.64 | 1 | Foundations for the issue features |
| PERF.65 | 1 | Foundations for the issue features |
| PERF.66 | 1 | Foundations for the issue features |
| PERF.67 | 1 | Foundations for the issue features |
| PERF.68 | 1 | Foundations for the issue features |
| PERF.69 | 1 | Foundations for the issue features |
| PERF.70 | 1 | Foundations for the issue features |
| PERF.71 | 1 | Foundations for the issue features |
| PERF.73 | 1 | Foundations for the issue features |
| PERF.74 | 1 | Foundations for the issue features |
| PERF.76 | 1 | Foundations for the issue features |
| PERF.77 | 1 | Foundations for the issue features |
| PERF.78 | 1 | Foundations for the issue features |
| PERF.79 | 1 | Foundations for the issue features |
| PERF.80 | 1 | Foundations for the issue features |
| PERF.81 | 1 | Foundations for the issue features |
| PERF.82 | 1 | Foundations for the issue features |
| PERF.83 | 1 | Bugs from the GitHub issues |
| PERF.84 | 1 | Bugs from the GitHub issues |
| PERF.85 | 1 | Bugs from the GitHub issues |
| PERF.87 | 1 | Bugs from the GitHub issues |
| PERF.89 | 1 | Bugs from the GitHub issues |
| PERF.90 | 1 | Bugs from the GitHub issues |
| PERF.92 | 1 | Bugs from the GitHub issues |
| PERF.93 | 1 | Foundations for the issue features |
| PERF.94 | 1 | Foundations for the issue features |
| PERF.95 | 1 | Foundations for the issue features |
| PERF.96 | 1 | Foundations for the issue features |
| PERF.101 | 1 | Maps |
| PERF.102 | 1 | Maps |
| PERF.103 | 2 | Orbital updates |
| POP.7 | 1 | Tech levels and facilities |
| POP.8 | 1 | Tech levels and facilities |
| POP.9 | 1 | Tech levels and facilities |
| POP.10 | 1 | Tech levels and facilities |
| POP.12 | 1 | Tech levels and facilities |
| SEC.29 | 0 | Groundwork: layout and libraries |
| SEC.30 | 0 | Groundwork: layout and libraries |
| SEC.31 | 0 | Bugfixes: maps and pages |
| TEST.70 | 0 | Map groundwork |
| TEST.71 | 0 | Bugfixes: ops and flakes |
| TEST.72 | 0 | Groundwork: layout and libraries |
| TEST.73 | 0 | Parallel path (top priority) |
| TEST.74 | 0 | Parallel path (top priority) |
| TEST.75 | 0 | Bugfixes: prevalence |
| TEST.76 | 0 | Parallel path (top priority) |
| TEST.78 | 0 | Bugfixes: ops and flakes |
| TEST.79 | 0 | Routing groundwork |
| TEST.80 | 0 | Bugfixes: ops and flakes |
| TEST.81 | 0 | Bugfixes: CI red |
| TEST.82 | 0 | Bugfixes: ops and flakes |
| TEST.83 | 0 | Groundwork: layout and libraries |
| TEST.84 | 0 | Bugfixes: ops and flakes |
| TEST.85 | 0 | Binary pairs and single-system forcing |
| TEST.86 | 0 | Bugfixes: ops and flakes |
| TEST.87 | 0 | Bugfixes: CI red |
| TEST.88 | 0 | Bugfixes: ops and flakes |
| TEST.89 | 0 | Bugfixes: ops and flakes |
| TEST.90 | 0 | Bugfixes: ops and flakes |
| TEST.91 | 0 | Bugfixes: ops and flakes |
| TEST.92 | 0 | Bugfixes: ops and flakes |
| TEST.93 | 0 | Bugfixes: ops and flakes |
| TEST.94 | 0 | Bugfixes: ops and flakes |
| TEST.95 | 0 | Bugfixes: ops and flakes |
| TEST.96 | 0 | Bugfixes: ops and flakes |
| TEST.97 | 0 | Bugfixes: ops and flakes |
| TEST.98 | 0 | Bugfixes: ops and flakes |
| TEST.99 | 0 | Bugfixes: ops and flakes |
| TEST.100 | 0 | Bugfixes: ops and flakes |
| TEST.101 | 0 | Bugfixes: ops and flakes |
| TEST.102 | 0 | Bugfixes: generation |
| TEST.103 | 0 | Bugfixes: maps and pages |
| TEST.104 | 0 | Bugfixes: ops and flakes |
| TEST.105 | 0 | Bugfixes: maps and pages |
| TEST.106 | 0 | Bugfixes: maps and pages |
| TEST.107 | 0 | Bugfixes: maps and pages |
| TEST.108 | 0 | Bugfixes: ops and flakes |
| TEST.109 | 0 | Bugfixes: ops and flakes |
| TEST.110 | 1 | Database consistency check |
| TEST.112 | 1 | Bugs from the GitHub issues |
| TEST.113 | 1 | Bugs from the GitHub issues |
| TEST.114 | 1 | Bugs from the GitHub issues |
| TEST.115 | 1 | Bugs from the GitHub issues |
| TEST.117 | 1 | Bugs from the GitHub issues |
| TEST.118 | 1 | Bugs from the GitHub issues |
| TEST.119 | 1 | Bugs from the GitHub issues |
| TEST.120 | 1 | Bugs from the GitHub issues |
| TEST.121 | 1 | Bugs from the GitHub issues |
| TEST.122 | 1 | Generate page |
| TEST.123 | 1 | Bugs from the GitHub issues |
| TEST.124 | 1 | Maps |
| TEST.125 | 1 | Bugs from the GitHub issues |
| TEST.132 | 1 | Bugs from the GitHub issues |
| TEST.133 | 1 | Foundations for the issue features |
| TEST.136 | 1 | Foundations for the issue features |
| TEST.139 | 1 | Bugs from the GitHub issues |
| TEST.142 | 1 | Bugs from the GitHub issues |
| TEST.143 | 1 | Bugs from the GitHub issues |
| TEST.145 | 1 | Bugs from the GitHub issues |
| TEST.146 | 1 | Bugs from the GitHub issues |
| TEST.147 | 1 | Bugs from the GitHub issues |
| TEST.149 | 1 | Bugs from the GitHub issues |
| TEST.150 | 1 | Bugs from the GitHub issues |
| TEST.151 | 1 | Bugs from the GitHub issues |
| TEST.152 | 1 | Bugs from the GitHub issues |
| TEST.153 | 1 | Bugs from the GitHub issues |
| TEST.154 | 1 | Bugs from the GitHub issues |
| TEST.157 | 1 | Bugs from the GitHub issues |
| TEST.158 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| TEST.163 | 1 | Bugs from the GitHub issues |
| UX.2 | 0 | Groundwork: web components |
| UX.21 | 0 | Groundwork: UX sweep |
| UX.22 | 1 | Pages |
| UX.23 | 1 | Pages |
| UX.24 | 0 | Bugfixes: maps and pages |
| UX.25 | 0 | Bugfixes: maps and pages |
| UX.26 | 0 | Groundwork: web components |
| UX.27 | 0 | Groundwork: web components |
| UX.28 | 0 | Bugfixes: maps and pages |
| UX.29 | 0 | Bugfixes: maps and pages |
| UX.30 | 2 | Pages |
| UX.31 | 0 | Groundwork: web components |
| UX.33 | 0 | Groundwork: web components |
| UX.34 | 0 | Bugfixes: console and progress |
| UX.35 | 1 | Routing |
| UX.36 | 0 | Bugfixes: maps and pages |
| UX.37 | 0 | Groundwork: UX sweep |
| UX.38 | 0 | Groundwork: nebula shapes |
| UX.39 | 0 | Groundwork: layout and libraries |
| UX.40 | 0 | Groundwork: web components |
| UX.41 | 0 | Groundwork: web components |
| UX.42 | 1 | Pages |
| UX.44 | 0 | Bugfixes: maps and pages |
| UX.50 | 0 | Groundwork: UX sweep |
| UX.51 | 0 | Groundwork: UX sweep |
| UX.52 | 0 | Groundwork: UX sweep |
| UX.53 | 0 | Groundwork: UX sweep |
| UX.54 | 0 | Groundwork: UX sweep |
| UX.55 | 0 | Groundwork: UX sweep |
| UX.56 | 0 | Groundwork: UX sweep |
| UX.57 | 0 | Groundwork: UX sweep |
| UX.58 | 0 | Groundwork: UX sweep |
| UX.59 | 0 | Groundwork: UX sweep |
| UX.60 | 0 | Groundwork: UX sweep |
| UX.61 | 0 | Groundwork: UX sweep |
| UX.62 | 0 | Groundwork: UX sweep |
| UX.63 | 0 | Groundwork: UX sweep |
| UX.64 | 0 | Groundwork: UX sweep |
| UX.65 | 0 | Groundwork: UX sweep |
| UX.66 | 0 | Groundwork: UX sweep |
| UX.67 | 0 | Groundwork: UX sweep |
| UX.68 | 0 | Groundwork: UX sweep |
| UX.69 | 0 | Groundwork: UX sweep |
| UX.70 | 0 | Groundwork: UX sweep |
| UX.71 | 0 | Groundwork: UX sweep |
| UX.72 | 0 | Groundwork: UX sweep |
| UX.73 | 0 | Groundwork: UX sweep |
| UX.74 | 0 | Groundwork: UX sweep |
| UX.75 | 0 | Groundwork: UX sweep |
| UX.76 | 1 | Bugs from the GitHub issues |
| UX.77 | 1 | Bugs from the GitHub issues |
| UX.78 | 1 | Pages |
| UX.79 | 1 | Pages |
| UX.80 | 1 | Pages |
| UX.81 | 1 | Pages |
| UX.83 | 1 | Bugs from the GitHub issues |
| UX.84 | 1 | Bugs from the GitHub issues |
| UX.85 | 1 | Bugs from the GitHub issues |
| UX.86 | 1 | Bugs from the GitHub issues |
| UX.87 | 1 | Maps |
| UX.88 | 1 | Maps |
| UX.89 | 1 | Maps |
| UX.90 | 1 | Maps |
| UX.91 | 1 | Maps |
| UX.93 | 1 | Routing |
| UX.95 | 1 | Foundations for the issue features |
| UX.96 | 1 | Foundations for the issue features |
| UX.97 | 1 | Foundations for the issue features |
| UX.98 | 1 | Bugs from the GitHub issues |
| UX.99 | 1 | Bugs from the GitHub issues |
| UX.103 | 1 | Maps |
| UX.104 | 1 | Maps |
| UX.106 | 1 | Maps |
| UX.111 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| UX.112 | 1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| UX.117 | 1 | Maps |
| VIEW.5 | 2 | View |
| GEN.281 | 2 | Anomalies |
