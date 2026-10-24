# Progress ledger: the version each finished item belonged to

Boss (2026-10-16, 11:09 AM EDT) asked for a progress bar for every group of the plan; the same day (1:07 PM EDT) the plan became a version map, so the bars are one per version. A version file drops an item's row when the item ships, so the version of a finished item has to be kept somewhere; this file keeps it. The rows were read back from the git history on 2026-10-16 and sorted into the version the item would sit in on the version map (the old Phase 0 is 8.0; the old threads Orbital updates, Daily maintenance, Routing, Courses, Classes, Habitability, Nebula planets, Recipes and View went to the version that carries that feature, everything else to 8.1). Every later retirement adds its row in the same PR (the retire step copies the version from the version file row before deleting it).

`scripts/build_todo_docs.py` counts, per version, the items here whose number-map status is done (finished), the items still open in the version files, and the dropped ones (closed, superseded, folded or not reproduced; they are not counted either way). Progress is finished / (finished + open). A finished item with no row here was built before the plan began on 2026-10-01 and is counted under "built before the plan".

The build fails if an item is done in the number map with a PR above the baseline below and has no row here, so a retire PR cannot forget it.

Baseline: PR #1817.

| ID | Version | Last thread |
|---|---|---|
| ADM.14 | 8.0 | Groundwork: web components |
| ADM.15 | 8.1 | Queue |
| ADM.16 | 8.0 | Bugfixes: prevalence |
| ADM.17 | 8.1 | Reproducible galaxies |
| ADM.18 | 8.1 | Reproducible galaxies |
| ADM.19 | 8.3 | Daily maintenance |
| ADM.20 | 8.3 | Daily maintenance |
| ADM.21 | 8.0 | Groundwork: data model |
| ADM.22 | 8.0 | Groundwork: queue, logs and caches |
| ADM.23 | 8.0 | Groundwork: queue, logs and caches |
| ADM.24 | 8.0 | Groundwork: queue, logs and caches |
| ADM.25 | 8.0 | Groundwork: queue, logs and caches |
| ADM.26 | 8.0 | Groundwork: queue, logs and caches |
| ADM.27 | 8.0 | Bugfixes: generation |
| ADM.28 | 8.1 | Generate page |
| ADM.29 | 8.1 | Generate page |
| ADM.30 | 8.1 | Generate page |
| ADM.31 | 8.1 | Generate page |
| ADM.33 | 8.0 | Bugfixes: console and progress |
| ADM.34 | 8.0 | Groundwork: web components |
| ADM.37 | 8.0 | Bugfixes: prevalence |
| ADM.38 | 8.0 | Bugfixes: console and progress |
| ADM.39 | 8.1 | Bugs from the GitHub issues |
| ADM.40 | 8.1 | Bugs from the GitHub issues |
| ADM.41 | 8.1 | Bugs from the GitHub issues |
| ADM.42 | 8.1 | Foundations for the issue features |
| ADM.45 | 8.1 | Generate page |
| ADM.46 | 8.1 | Generate page |
| ADM.47 | 8.1 | Generate page |
| ADM.48 | 8.1 | Foundations for the issue features |
| ADM.49 | 8.1 | Generate page |
| ADM.51 | 8.1 | Bugs from the GitHub issues |
| ADM.53 | 8.1 | Bugs from the GitHub issues |
| ADM.54 | 8.1 | Foundations for the issue features |
| API.9 | 8.1 | API |
| API.16 | 8.1 | API |
| API.20 | 8.9 | Recipes |
| API.21 | 8.9 | Recipes |
| API.22 | 8.1 | API |
| API.23 | 8.1 | API |
| API.25 | 8.1 | Maps |
| DB.2 | 8.0 | Database fixes |
| DB.3 | 8.0 | Database fixes |
| DB.4 | 8.0 | Database fixes |
| DB.5 | 8.0 | Database fixes |
| DB.6 | 8.0 | Parallel path (top priority) |
| DB.7 | 8.1 | Reproducible galaxies |
| DB.8 | 8.1 | Database consistency check |
| DB.10 | 8.1 | Database consistency check |
| DB.11 | 8.0 | Groundwork: data model |
| DB.12 | 8.0 | Bugfixes: CI red |
| DB.13 | 8.0 | Groundwork: data model |
| DB.14 | 8.1 | Maps |
| DB.15 | 8.1 | Foundations for the issue features |
| DB.19 | 8.1 | Foundations for the issue features |
| DB.20 | 8.1 | Database consistency check |
| DB.21 | 8.1 | Database consistency check |
| DB.22 | 8.1 | API |
| DB.23 | 8.1 | Foundations for the issue features |
| DB.24 | 8.1 | Maps |
| DB.25 | 8.1 | Bugs from the GitHub issues |
| DB.26 | 8.1 | Foundations for the issue features |
| DB.27 | 8.1 | Features from the GitHub issues |
| DB.28 | 8.1 | Tech levels and facilities |
| DB.32 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| DB.34 | 8.1 | Picker |
| DB.35 | 8.1 | Picker |
| DOC.5 | 8.1 | Database consistency check |
| DOC.17 | 8.1 | API |
| DOC.21 | 8.1 | Bugs from the GitHub issues |
| DOC.22 | 8.1 | Bugs from the GitHub issues |
| DOC.23 | 8.1 | Bugs from the GitHub issues |
| DOC.25 | 8.1 | Bugs from the GitHub issues |
| DOC.26 | 8.1 | Bugs from the GitHub issues |
| DOC.27 | 8.1 | Bugs from the GitHub issues |
| DOC.28 | 8.1 | Bugs from the GitHub issues |
| DOC.29 | 8.1 | Bugs from the GitHub issues |
| DOC.30 | 8.1 | Bugs from the GitHub issues |
| DOC.31 | 8.1 | Bugs from the GitHub issues |
| DOC.32 | 8.1 | Bugs from the GitHub issues |
| DOC.33 | 8.1 | Bugs from the GitHub issues |
| DOC.34 | 8.1 | Bugs from the GitHub issues |
| DOC.35 | 8.1 | Bugs from the GitHub issues |
| DOC.36 | 8.1 | Bugs from the GitHub issues |
| DOC.37 | 8.1 | Bugs from the GitHub issues |
| DOC.38 | 8.1 | Bugs from the GitHub issues |
| DOC.39 | 8.1 | Bugs from the GitHub issues |
| DOC.40 | 8.1 | Bugs from the GitHub issues |
| DOC.41 | 8.1 | Bugs from the GitHub issues |
| DOC.42 | 8.1 | Bugs from the GitHub issues |
| DOC.43 | 8.1 | Bugs from the GitHub issues |
| DOC.44 | 8.1 | Bugs from the GitHub issues |
| DOC.46 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.25 | 8.0 | Physics bugs |
| GEN.31 | 8.0 | Galaxy geometry and names |
| GEN.32 | 8.0 | Parallel path (top priority) |
| GEN.34 | 8.0 | Physics bugs |
| GEN.35 | 8.0 | Physics bugs |
| GEN.36 | 8.0 | Physics bugs |
| GEN.37 | 8.0 | Physics bugs |
| GEN.38 | 8.0 | Generation bugs |
| GEN.39 | 8.0 | Parallel path (top priority) |
| GEN.41 | 8.1 | Galaxy gen |
| GEN.43 | 8.1 | Backfill |
| GEN.44 | 8.0 | Sector stats and colors |
| GEN.45 | 8.0 | Physics bugs |
| GEN.46 | 8.0 | Galaxy geometry and names |
| GEN.47 | 8.0 | Generation bugs |
| GEN.48 | 8.0 | Bugfixes: prevalence |
| GEN.49 | 8.0 | Binary pairs and single-system forcing |
| GEN.50 | 8.0 | Binary pairs and single-system forcing |
| GEN.51 | 8.0 | Binary pairs and single-system forcing |
| GEN.52 | 8.0 | Bugfixes: prevalence |
| GEN.53 | 8.0 | Binary pairs and single-system forcing |
| GEN.54 | 8.0 | Binary pairs and single-system forcing |
| GEN.56 | 8.1 | Reproducible galaxies |
| GEN.58 | 8.1 | Reproducible galaxies |
| GEN.59 | 8.1 | Reproducible galaxies |
| GEN.60 | 8.0 | Generation bugs |
| GEN.61 | 8.3 | Daily maintenance |
| GEN.62 | 8.0 | Binary pairs and single-system forcing |
| GEN.63 | 8.1 | Reproducible galaxies |
| GEN.64 | 8.1 | Interstellar IDs |
| GEN.65 | 8.0 | Bugfixes: generation |
| GEN.66 | 8.0 | Groundwork: data model |
| GEN.67 | 8.0 | Groundwork: names from IDs |
| GEN.68 | 8.0 | Groundwork: names from IDs |
| GEN.69 | 8.0 | Groundwork: names from IDs |
| GEN.70 | 8.0 | Groundwork: names from IDs |
| GEN.71 | 8.0 | Groundwork: names from IDs |
| GEN.72 | 8.0 | Groundwork: names from IDs |
| GEN.73 | 8.0 | Groundwork: names from IDs |
| GEN.74 | 8.0 | Groundwork: data model |
| GEN.75 | 8.0 | Groundwork: nebula shapes |
| GEN.76 | 8.0 | Bugfixes: generation |
| GEN.77 | 8.0 | Bugfixes: generation |
| GEN.78 | 8.0 | Bugfixes: generation |
| GEN.79 | 8.0 | Bugfixes: generation |
| GEN.80 | 8.0 | Bugfixes: generation |
| GEN.81 | 8.0 | Bugfixes: console and progress |
| GEN.82 | 8.0 | Bugfixes: generation |
| GEN.83 | 8.4 | Habitability |
| GEN.84 | 8.4 | Habitability |
| GEN.85 | 8.4 | Habitability |
| GEN.86 | 8.4 | Habitability |
| GEN.87 | 8.4 | Habitability |
| GEN.88 | 8.4 | Habitability |
| GEN.89 | 8.4 | Habitability |
| GEN.94 | 8.4 | Nebula planets |
| GEN.96 | 8.1 | Generate page |
| GEN.97 | 8.1 | Generate page |
| GEN.98 | 8.1 | Galaxy gen |
| GEN.100 | 8.1 | Galaxy gen |
| GEN.101 | 8.1 | Galaxy gen |
| GEN.104 | 8.3 | Orbital updates |
| GEN.106 | 8.3 | Orbital updates |
| GEN.107 | 8.3 | Orbital updates |
| GEN.108 | 8.3 | Orbital updates |
| GEN.113 | 8.1 | Anomalies |
| GEN.114 | 8.1 | Anomalies |
| GEN.116 | 8.0 | Bugfixes: generation |
| GEN.117 | 8.0 | Groundwork: map engine |
| GEN.118 | 8.0 | Bugfixes: generation |
| GEN.119 | 8.0 | Bugfixes: generation |
| GEN.120 | 8.1 | Galaxy gen |
| GEN.121 | 8.3 | Orbital updates |
| GEN.122 | 8.3 | Orbital updates |
| GEN.123 | 8.3 | Orbital updates |
| GEN.124 | 8.3 | Orbital updates |
| GEN.125 | 8.3 | Orbital updates |
| GEN.126 | 8.0 | Groundwork: data model |
| GEN.127 | 8.1 | Bugs from the GitHub issues |
| GEN.128 | 8.1 | Foundations for the issue features |
| GEN.131 | 8.1 | Foundations for the issue features |
| GEN.132 | 8.1 | Foundations for the issue features |
| GEN.133 | 8.1 | Foundations for the issue features |
| GEN.137 | 8.5 | Routing |
| GEN.138 | 8.3 | Orbital updates |
| GEN.147 | 8.4 | Classes |
| GEN.148 | 8.4 | Habitability |
| GEN.149 | 8.4 | Nebula planets |
| GEN.150 | 8.1 | Galaxy gen |
| GEN.151 | 8.4 | Nebula planets |
| GEN.152 | 8.1 | Foundations for the issue features |
| GEN.153 | 8.1 | Anomalies |
| GEN.154 | 8.1 | Anomalies |
| GEN.155 | 8.1 | Galaxy gen |
| GEN.156 | 8.1 | Reproducible galaxies |
| GEN.158 | 8.1 | Anomalies |
| GEN.159 | 8.1 | Anomalies |
| GEN.160 | 8.1 | Anomalies |
| GEN.165 | 8.1 | Features from the GitHub issues |
| GEN.166 | 8.1 | Foundations for the issue features |
| GEN.167 | 8.1 | Foundations for the issue features |
| GEN.168 | 8.1 | Foundations for the issue features |
| GEN.169 | 8.1 | Foundations for the issue features |
| GEN.170 | 8.1 | Database consistency check |
| GEN.171 | 8.1 | Database consistency check |
| GEN.172 | 8.1 | Database consistency check |
| GEN.173 | 8.1 | Database consistency check |
| GEN.174 | 8.1 | Database consistency check |
| GEN.175 | 8.1 | Database consistency check |
| GEN.176 | 8.1 | Database consistency check |
| GEN.177 | 8.4 | Habitability |
| GEN.178 | 8.4 | Habitability |
| GEN.182 | 8.1 | Bugs from the GitHub issues |
| GEN.183 | 8.1 | Foundations for the issue features |
| GEN.184 | 8.1 | Foundations for the issue features |
| GEN.185 | 8.1 | Foundations for the issue features |
| GEN.186 | 8.1 | Generate page |
| GEN.187 | 8.1 | Foundations for the issue features |
| GEN.188 | 8.1 | Foundations for the issue features |
| GEN.189 | 8.1 | Maps |
| GEN.190 | 8.1 | Generate page |
| GEN.191 | 8.1 | Generate page |
| GEN.192 | 8.1 | Foundations for the issue features |
| GEN.194 | 8.1 | Maps |
| GEN.195 | 8.1 | Maps |
| GEN.196 | 8.1 | Maps |
| GEN.200 | 8.1 | Foundations for the issue features |
| GEN.201 | 8.1 | Foundations for the issue features |
| GEN.203 | 8.1 | Foundations for the issue features |
| GEN.204 | 8.1 | Foundations for the issue features |
| GEN.205 | 8.1 | Foundations for the issue features |
| GEN.206 | 8.1 | Bugs from the GitHub issues |
| GEN.207 | 8.1 | Maps |
| GEN.208 | 8.1 | Maps |
| GEN.211 | 8.1 | Bugs from the GitHub issues |
| GEN.213 | 8.4 | Habitability |
| GEN.214 | 8.1 | Maps |
| GEN.215 | 8.1 | Maps |
| GEN.216 | 8.1 | Maps |
| GEN.217 | 8.1 | Maps |
| GEN.222 | 8.1 | Maps |
| GEN.223 | 8.1 | Maps |
| GEN.224 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.225 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.226 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.227 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.228 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| GEN.230 | 8.1 | Maps |
| GEN.232 | 8.1 | Maps |
| GEN.233 | 8.1 | Orbital updates |
| GEN.234 | 8.1 | Orbital updates |
| GEN.235 | 8.1 | Orbital updates |
| GEN.250 | 8.1 | Orbital updates |
| GEN.251 | 8.1 | Orbital updates |
| GEN.254 | 8.1 | Orbital updates |
| GEN.260 | 8.1 | Orbital updates |
| GEN.261 | 8.1 | Orbital updates |
| GEN.262 | 8.1 | Orbital updates |
| GEN.265 | 8.1 | Orbital updates |
| GEN.266 | 8.1 | Orbital updates |
| GEN.268 | 8.1 | Orbital updates |
| GEN.277 | 8.1 | Anomalies |
| GEN.278 | 8.1 | Anomalies |
| GEN.279 | 8.1 | Anomalies |
| GEN.288 | 8.1 | Maps |
| MAP.52 | 8.0 | Galaxy map picker and arc |
| MAP.53 | 8.0 | Galaxy Map drill-down |
| MAP.54 | 8.0 | Galaxy Map drill-down |
| MAP.55 | 8.0 | Galaxy map picker and arc |
| MAP.56 | 8.0 | Galaxy Map drill-down |
| MAP.57 | 8.0 | System Map |
| MAP.58 | 8.1 | Galaxy Map |
| MAP.60 | 8.0 | Galaxy map picker and arc |
| MAP.61 | 8.0 | Groundwork: map engine |
| MAP.62 | 8.1 | 3D system |
| MAP.63 | 8.0 | Map groundwork |
| MAP.64 | 8.0 | Map groundwork |
| MAP.65 | 8.0 | Groundwork: map engine |
| MAP.66 | 8.0 | Groundwork: map engine |
| MAP.67 | 8.0 | Groundwork: map engine |
| MAP.68 | 8.0 | Groundwork: map engine |
| MAP.69 | 8.1 | 3D system |
| MAP.70 | 8.1 | 3D system |
| MAP.71 | 8.1 | 3D system |
| MAP.72 | 8.1 | 3D system |
| MAP.73 | 8.1 | 3D system |
| MAP.74 | 8.1 | 3D system |
| MAP.76 | 8.0 | Galaxy Map drill-down |
| MAP.77 | 8.0 | Galaxy Map drill-down |
| MAP.78 | 8.0 | Galaxy Map drill-down |
| MAP.79 | 8.0 | Groundwork: map engine |
| MAP.80 | 8.0 | Sector stats and colors |
| MAP.81 | 8.0 | Map groundwork |
| MAP.82 | 8.0 | Map groundwork |
| MAP.83 | 8.0 | Map groundwork |
| MAP.84 | 8.0 | Map groundwork |
| MAP.85 | 8.0 | Galaxy map picker and arc |
| MAP.86 | 8.0 | Sector stats and colors |
| MAP.87 | 8.0 | Map groundwork |
| MAP.88 | 8.0 | System Map |
| MAP.89 | 8.1 | System Map |
| MAP.90 | 8.0 | Galaxy geometry and names |
| MAP.91 | 8.0 | Galaxy Map follow-ups |
| MAP.92 | 8.0 | System Map |
| MAP.93 | 8.0 | Galaxy Map follow-ups |
| MAP.94 | 8.0 | Galaxy Map follow-ups |
| MAP.95 | 8.0 | Groundwork: map engine |
| MAP.96 | 8.0 | Galaxy Map drill-down |
| MAP.97 | 8.0 | Galaxy Map drill-down |
| MAP.98 | 8.0 | Galaxy Map drill-down |
| MAP.99 | 8.0 | Galaxy Map drill-down |
| MAP.100 | 8.0 | Galaxy Map drill-down |
| MAP.101 | 8.0 | Galaxy Map drill-down |
| MAP.102 | 8.0 | Groundwork: map engine |
| MAP.103 | 8.0 | Groundwork: nebula shapes |
| MAP.104 | 8.0 | Groundwork: nebula shapes |
| MAP.105 | 8.0 | Groundwork: nebula shapes |
| MAP.106 | 8.0 | Groundwork: map engine |
| MAP.107 | 8.0 | Groundwork: map engine |
| MAP.108 | 8.0 | Groundwork: map engine |
| MAP.109 | 8.0 | Groundwork: map engine |
| MAP.110 | 8.0 | Groundwork: map engine |
| MAP.111 | 8.0 | Groundwork: map engine |
| MAP.112 | 8.0 | Groundwork: map engine |
| MAP.113 | 8.0 | Groundwork: map engine |
| MAP.114 | 8.0 | Bugfixes: CI red |
| MAP.115 | 8.0 | Bugfixes: maps and pages |
| MAP.116 | 8.0 | Groundwork: map engine |
| MAP.117 | 8.0 | Bugfixes: maps and pages |
| MAP.118 | 8.0 | Bugfixes: maps and pages |
| MAP.123 | 8.1 | Maps |
| MAP.124 | 8.1 | Maps |
| MAP.125 | 8.1 | Infinite zoom |
| MAP.126 | 8.1 | 3D system |
| MAP.127 | 8.0 | Groundwork: map engine |
| MAP.128 | 8.1 | Maps |
| MAP.129 | 8.1 | Maps |
| MAP.130 | 8.1 | Maps |
| MAP.131 | 8.1 | Galaxy Map |
| MAP.133 | 8.0 | Bugfixes: maps and pages |
| MAP.134 | 8.1 | Maps |
| MAP.135 | 8.1 | Bugs from the GitHub issues |
| MAP.136 | 8.1 | Bugs from the GitHub issues |
| MAP.137 | 8.1 | Bugs from the GitHub issues |
| MAP.138 | 8.1 | Bugs from the GitHub issues |
| MAP.144 | 8.1 | Galaxy Map |
| MAP.148 | 8.1 | Maps |
| MAP.149 | 8.1 | Maps |
| MAP.150 | 8.1 | Maps |
| MAP.151 | 8.1 | Maps |
| MAP.153 | 8.1 | Maps |
| MAP.154 | 8.1 | Maps |
| MAP.155 | 8.1 | Maps |
| MAP.157 | 8.1 | Maps |
| MAP.158 | 8.1 | Maps |
| MAP.159 | 8.1 | Maps |
| MAP.161 | 8.1 | Maps |
| MAP.162 | 8.1 | Maps |
| MAP.163 | 8.1 | Maps |
| MAP.164 | 8.1 | Maps |
| MAP.165 | 8.1 | Maps |
| MAP.166 | 8.1 | Maps |
| MAP.175 | 8.1 | Maps |
| MAP.182 | 8.1 | Orbital updates |
| MAP.184 | 8.1 | Bugs from the GitHub issues |
| MAP.185 | 8.1 | Bugs from the GitHub issues |
| MAP.186 | 8.1 | Bugs from the GitHub issues |
| MAP.187 | 8.1 | Bugs from the GitHub issues |
| MAP.188 | 8.1 | Bugs from the GitHub issues |
| MAP.196 | 8.1 | Bugs from the GitHub issues |
| MAP.197 | 8.1 | Maps |
| MAP.198 | 8.1 | Maps |
| MAP.199 | 8.1 | Maps |
| MAP.200 | 8.1 | Maps |
| MAP.201 | 8.1 | Maps |
| MAP.202 | 8.1 | Maps |
| MAP.203 | 8.1 | Maps |
| MAP.204 | 8.1 | Maps |
| MAP.205 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| MAP.206 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| MAP.208 | 8.1 | Orbital updates |
| MAP.211 | 8.1 | Orbital updates |
| NAV.3 | 8.1 | Picker |
| NAV.7 | 8.0 | Groundwork: map engine |
| NAV.8 | 8.1 | References |
| NAV.9 | 8.1 | References |
| NAV.10 | 8.5 | Routing |
| NAV.11 | 8.5 | Routing |
| NAV.12 | 8.5 | Routing |
| NAV.13 | 8.0 | Groundwork: map engine |
| NAV.14 | 8.0 | Groundwork: map engine |
| NAV.15 | 8.0 | Groundwork: map engine |
| NAV.16 | 8.1 | References |
| NAV.20 | 8.5 | Courses |
| NAV.24 | 8.5 | Courses |
| NAV.29 | 8.0 | Groundwork: map engine |
| NAV.30 | 8.0 | Map groundwork |
| NAV.31 | 8.0 | Galaxy Map follow-ups |
| NAV.32 | 8.0 | Groundwork: map engine |
| NAV.33 | 8.0 | Groundwork: map engine |
| NAV.34 | 8.0 | Routing groundwork |
| NAV.35 | 8.0 | Routing guarantee |
| NAV.37 | 8.5 | Courses |
| NAV.38 | 8.0 | Galaxy geometry and names |
| NAV.40 | 8.0 | Galaxy Map follow-ups |
| NAV.41 | 8.0 | Bugfixes: maps and pages |
| NAV.42 | 8.5 | Routing |
| NAV.43 | 8.1 | Nearby search |
| NAV.44 | 8.1 | Nearby search |
| NAV.46 | 8.0 | Groundwork: map engine |
| NAV.47 | 8.5 | Routing |
| NAV.48 | 8.5 | Routing |
| NAV.50 | 8.1 | References |
| NAV.52 | 8.5 | Routing |
| NAV.53 | 8.5 | Routing |
| OPS.6 | 8.0 | Bugfixes: ops and flakes |
| OPS.7 | 8.0 | Bugfixes: ops and flakes |
| OPS.8 | 8.1 | Reproducible galaxies |
| OPS.9 | 8.0 | Bugfixes: console and progress |
| OPS.10 | 8.0 | Parallel path (top priority) |
| OPS.11 | 8.1 | Reproducible galaxies |
| OPS.12 | 8.1 | Reproducible galaxies |
| OPS.13 | 8.1 | Reproducible galaxies |
| OPS.14 | 8.1 | Reproducible galaxies |
| OPS.18 | 8.3 | Daily maintenance |
| OPS.19 | 8.0 | Groundwork: queue, logs and caches |
| OPS.20 | 8.0 | Groundwork: layout and libraries |
| OPS.21 | 8.0 | Groundwork: layout and libraries |
| OPS.22 | 8.0 | Groundwork: layout and libraries |
| OPS.23 | 8.0 | Groundwork: layout and libraries |
| OPS.24 | 8.0 | Groundwork: layout and libraries |
| OPS.25 | 8.0 | Bugfixes: CI red |
| OPS.26 | 8.0 | Bugfixes: ops and flakes |
| OPS.27 | 8.0 | Groundwork: queue, logs and caches |
| OPS.31 | 8.3 | Daily maintenance |
| OPS.32 | 8.3 | Daily maintenance |
| OPS.33 | 8.1 | Reproducible galaxies |
| OPS.34 | 8.3 | Daily maintenance |
| OPS.36 | 8.1 | Foundations for the issue features |
| OPS.38 | 8.1 | Bugs from the GitHub issues |
| OPS.39 | 8.1 | Generate page |
| OPS.40 | 8.1 | Bugs from the GitHub issues |
| OPS.41 | 8.1 | Bookmarks |
| OPS.42 | 8.1 | Bugs from the GitHub issues |
| OPS.43 | 8.1 | Bugs from the GitHub issues |
| OPS.44 | 8.1 | Bugs from the GitHub issues |
| OPS.45 | 8.1 | Bugs from the GitHub issues |
| OPS.46 | 8.1 | Foundations for the issue features |
| OPS.47 | 8.1 | Foundations for the issue features |
| OPS.48 | 8.1 | Foundations for the issue features |
| OPS.49 | 8.1 | Foundations for the issue features |
| OPS.50 | 8.1 | Foundations for the issue features |
| OPS.51 | 8.1 | Bookmarks |
| OPS.52 | 8.1 | Bookmarks |
| OPS.53 | 8.1 | Bookmarks |
| OPS.54 | 8.1 | Bookmarks |
| PERF.1 | 8.0 | Sector stats and colors |
| PERF.11 | 8.0 | Sector stats and colors |
| PERF.19 | 8.0 | Groundwork: queue, logs and caches |
| PERF.21 | 8.0 | Parallel path (top priority) |
| PERF.22 | 8.0 | Parallel path (top priority) |
| PERF.23 | 8.0 | Parallel path (top priority) |
| PERF.24 | 8.0 | Groundwork: queue, logs and caches |
| PERF.25 | 8.0 | Groundwork: queue, logs and caches |
| PERF.26 | 8.0 | Bugfixes: generation |
| PERF.27 | 8.0 | Bugfixes: CI red |
| PERF.28 | 8.0 | Bugfixes: console and progress |
| PERF.31 | 8.1 | Foundations for the issue features |
| PERF.32 | 8.1 | Foundations for the issue features |
| PERF.33 | 8.1 | Foundations for the issue features |
| PERF.34 | 8.1 | Bugs from the GitHub issues |
| PERF.37 | 8.1 | Features from the GitHub issues |
| PERF.38 | 8.1 | Foundations for the issue features |
| PERF.39 | 8.1 | Foundations for the issue features |
| PERF.40 | 8.1 | Foundations for the issue features |
| PERF.42 | 8.1 | Foundations for the issue features |
| PERF.43 | 8.1 | Foundations for the issue features |
| PERF.44 | 8.1 | Foundations for the issue features |
| PERF.45 | 8.1 | Foundations for the issue features |
| PERF.46 | 8.1 | Foundations for the issue features |
| PERF.49 | 8.1 | Foundations for the issue features |
| PERF.50 | 8.1 | Bugs from the GitHub issues |
| PERF.51 | 8.1 | Foundations for the issue features |
| PERF.52 | 8.1 | Bugs from the GitHub issues |
| PERF.53 | 8.1 | Bugs from the GitHub issues |
| PERF.54 | 8.1 | Foundations for the issue features |
| PERF.55 | 8.1 | Foundations for the issue features |
| PERF.56 | 8.1 | Maps |
| PERF.57 | 8.1 | Maps |
| PERF.58 | 8.1 | Maps |
| PERF.59 | 8.1 | Maps |
| PERF.60 | 8.1 | Maps |
| PERF.61 | 8.1 | Maps |
| PERF.62 | 8.1 | Foundations for the issue features |
| PERF.63 | 8.1 | Maps |
| PERF.64 | 8.1 | Foundations for the issue features |
| PERF.65 | 8.1 | Foundations for the issue features |
| PERF.66 | 8.1 | Foundations for the issue features |
| PERF.67 | 8.1 | Foundations for the issue features |
| PERF.68 | 8.1 | Foundations for the issue features |
| PERF.69 | 8.1 | Foundations for the issue features |
| PERF.70 | 8.1 | Foundations for the issue features |
| PERF.71 | 8.1 | Foundations for the issue features |
| PERF.73 | 8.1 | Foundations for the issue features |
| PERF.74 | 8.1 | Foundations for the issue features |
| PERF.76 | 8.1 | Foundations for the issue features |
| PERF.77 | 8.1 | Foundations for the issue features |
| PERF.78 | 8.1 | Foundations for the issue features |
| PERF.79 | 8.1 | Foundations for the issue features |
| PERF.80 | 8.1 | Foundations for the issue features |
| PERF.81 | 8.1 | Foundations for the issue features |
| PERF.82 | 8.1 | Foundations for the issue features |
| PERF.83 | 8.1 | Bugs from the GitHub issues |
| PERF.84 | 8.1 | Bugs from the GitHub issues |
| PERF.85 | 8.1 | Bugs from the GitHub issues |
| PERF.87 | 8.1 | Bugs from the GitHub issues |
| PERF.89 | 8.1 | Bugs from the GitHub issues |
| PERF.90 | 8.1 | Bugs from the GitHub issues |
| PERF.92 | 8.1 | Bugs from the GitHub issues |
| PERF.93 | 8.1 | Foundations for the issue features |
| PERF.94 | 8.1 | Foundations for the issue features |
| PERF.95 | 8.1 | Foundations for the issue features |
| PERF.96 | 8.1 | Foundations for the issue features |
| PERF.101 | 8.1 | Maps |
| PERF.102 | 8.1 | Maps |
| PERF.103 | 8.1 | Orbital updates |
| POP.7 | 8.1 | Tech levels and facilities |
| POP.8 | 8.1 | Tech levels and facilities |
| POP.9 | 8.1 | Tech levels and facilities |
| POP.10 | 8.1 | Tech levels and facilities |
| POP.12 | 8.1 | Tech levels and facilities |
| SEC.29 | 8.0 | Groundwork: layout and libraries |
| SEC.30 | 8.0 | Groundwork: layout and libraries |
| SEC.31 | 8.0 | Bugfixes: maps and pages |
| TEST.70 | 8.0 | Map groundwork |
| TEST.71 | 8.0 | Bugfixes: ops and flakes |
| TEST.72 | 8.0 | Groundwork: layout and libraries |
| TEST.73 | 8.0 | Parallel path (top priority) |
| TEST.74 | 8.0 | Parallel path (top priority) |
| TEST.75 | 8.0 | Bugfixes: prevalence |
| TEST.76 | 8.0 | Parallel path (top priority) |
| TEST.78 | 8.0 | Bugfixes: ops and flakes |
| TEST.79 | 8.0 | Routing groundwork |
| TEST.80 | 8.0 | Bugfixes: ops and flakes |
| TEST.81 | 8.0 | Bugfixes: CI red |
| TEST.82 | 8.0 | Bugfixes: ops and flakes |
| TEST.83 | 8.0 | Groundwork: layout and libraries |
| TEST.84 | 8.0 | Bugfixes: ops and flakes |
| TEST.85 | 8.0 | Binary pairs and single-system forcing |
| TEST.86 | 8.0 | Bugfixes: ops and flakes |
| TEST.87 | 8.0 | Bugfixes: CI red |
| TEST.88 | 8.0 | Bugfixes: ops and flakes |
| TEST.89 | 8.0 | Bugfixes: ops and flakes |
| TEST.90 | 8.0 | Bugfixes: ops and flakes |
| TEST.91 | 8.0 | Bugfixes: ops and flakes |
| TEST.92 | 8.0 | Bugfixes: ops and flakes |
| TEST.93 | 8.0 | Bugfixes: ops and flakes |
| TEST.94 | 8.0 | Bugfixes: ops and flakes |
| TEST.95 | 8.0 | Bugfixes: ops and flakes |
| TEST.96 | 8.0 | Bugfixes: ops and flakes |
| TEST.97 | 8.0 | Bugfixes: ops and flakes |
| TEST.98 | 8.0 | Bugfixes: ops and flakes |
| TEST.99 | 8.0 | Bugfixes: ops and flakes |
| TEST.100 | 8.0 | Bugfixes: ops and flakes |
| TEST.101 | 8.0 | Bugfixes: ops and flakes |
| TEST.102 | 8.0 | Bugfixes: generation |
| TEST.103 | 8.0 | Bugfixes: maps and pages |
| TEST.104 | 8.0 | Bugfixes: ops and flakes |
| TEST.105 | 8.0 | Bugfixes: maps and pages |
| TEST.106 | 8.0 | Bugfixes: maps and pages |
| TEST.107 | 8.0 | Bugfixes: maps and pages |
| TEST.108 | 8.0 | Bugfixes: ops and flakes |
| TEST.109 | 8.0 | Bugfixes: ops and flakes |
| TEST.110 | 8.1 | Database consistency check |
| TEST.112 | 8.1 | Bugs from the GitHub issues |
| TEST.113 | 8.1 | Bugs from the GitHub issues |
| TEST.114 | 8.1 | Bugs from the GitHub issues |
| TEST.115 | 8.1 | Bugs from the GitHub issues |
| TEST.117 | 8.1 | Bugs from the GitHub issues |
| TEST.118 | 8.1 | Bugs from the GitHub issues |
| TEST.119 | 8.1 | Bugs from the GitHub issues |
| TEST.120 | 8.1 | Bugs from the GitHub issues |
| TEST.121 | 8.1 | Bugs from the GitHub issues |
| TEST.122 | 8.1 | Generate page |
| TEST.123 | 8.1 | Bugs from the GitHub issues |
| TEST.124 | 8.1 | Maps |
| TEST.125 | 8.1 | Bugs from the GitHub issues |
| TEST.132 | 8.1 | Bugs from the GitHub issues |
| TEST.133 | 8.1 | Foundations for the issue features |
| TEST.136 | 8.1 | Foundations for the issue features |
| TEST.139 | 8.1 | Bugs from the GitHub issues |
| TEST.142 | 8.1 | Bugs from the GitHub issues |
| TEST.143 | 8.1 | Bugs from the GitHub issues |
| TEST.145 | 8.1 | Bugs from the GitHub issues |
| TEST.146 | 8.1 | Bugs from the GitHub issues |
| TEST.147 | 8.1 | Bugs from the GitHub issues |
| TEST.149 | 8.1 | Bugs from the GitHub issues |
| TEST.150 | 8.1 | Bugs from the GitHub issues |
| TEST.151 | 8.1 | Bugs from the GitHub issues |
| TEST.152 | 8.1 | Bugs from the GitHub issues |
| TEST.153 | 8.1 | Bugs from the GitHub issues |
| TEST.154 | 8.1 | Bugs from the GitHub issues |
| TEST.157 | 8.1 | Bugs from the GitHub issues |
| TEST.158 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| TEST.163 | 8.1 | Bugs from the GitHub issues |
| UX.2 | 8.0 | Groundwork: web components |
| UX.21 | 8.0 | Groundwork: UX sweep |
| UX.22 | 8.1 | Pages |
| UX.23 | 8.1 | Pages |
| UX.24 | 8.0 | Bugfixes: maps and pages |
| UX.25 | 8.0 | Bugfixes: maps and pages |
| UX.26 | 8.0 | Groundwork: web components |
| UX.27 | 8.0 | Groundwork: web components |
| UX.28 | 8.0 | Bugfixes: maps and pages |
| UX.29 | 8.0 | Bugfixes: maps and pages |
| UX.30 | 8.1 | Pages |
| UX.31 | 8.0 | Groundwork: web components |
| UX.33 | 8.0 | Groundwork: web components |
| UX.34 | 8.0 | Bugfixes: console and progress |
| UX.35 | 8.5 | Routing |
| UX.36 | 8.0 | Bugfixes: maps and pages |
| UX.37 | 8.0 | Groundwork: UX sweep |
| UX.38 | 8.0 | Groundwork: nebula shapes |
| UX.39 | 8.0 | Groundwork: layout and libraries |
| UX.40 | 8.0 | Groundwork: web components |
| UX.41 | 8.0 | Groundwork: web components |
| UX.42 | 8.1 | Pages |
| UX.44 | 8.0 | Bugfixes: maps and pages |
| UX.50 | 8.0 | Groundwork: UX sweep |
| UX.51 | 8.0 | Groundwork: UX sweep |
| UX.52 | 8.0 | Groundwork: UX sweep |
| UX.53 | 8.0 | Groundwork: UX sweep |
| UX.54 | 8.0 | Groundwork: UX sweep |
| UX.55 | 8.0 | Groundwork: UX sweep |
| UX.56 | 8.0 | Groundwork: UX sweep |
| UX.57 | 8.0 | Groundwork: UX sweep |
| UX.58 | 8.0 | Groundwork: UX sweep |
| UX.59 | 8.0 | Groundwork: UX sweep |
| UX.60 | 8.0 | Groundwork: UX sweep |
| UX.61 | 8.0 | Groundwork: UX sweep |
| UX.62 | 8.0 | Groundwork: UX sweep |
| UX.63 | 8.0 | Groundwork: UX sweep |
| UX.64 | 8.0 | Groundwork: UX sweep |
| UX.65 | 8.0 | Groundwork: UX sweep |
| UX.66 | 8.0 | Groundwork: UX sweep |
| UX.67 | 8.0 | Groundwork: UX sweep |
| UX.68 | 8.0 | Groundwork: UX sweep |
| UX.69 | 8.0 | Groundwork: UX sweep |
| UX.70 | 8.0 | Groundwork: UX sweep |
| UX.71 | 8.0 | Groundwork: UX sweep |
| UX.72 | 8.0 | Groundwork: UX sweep |
| UX.73 | 8.0 | Groundwork: UX sweep |
| UX.74 | 8.0 | Groundwork: UX sweep |
| UX.75 | 8.0 | Groundwork: UX sweep |
| UX.76 | 8.1 | Bugs from the GitHub issues |
| UX.77 | 8.1 | Bugs from the GitHub issues |
| UX.78 | 8.1 | Pages |
| UX.79 | 8.1 | Pages |
| UX.80 | 8.1 | Pages |
| UX.81 | 8.1 | Pages |
| UX.83 | 8.1 | Bugs from the GitHub issues |
| UX.84 | 8.1 | Bugs from the GitHub issues |
| UX.85 | 8.1 | Bugs from the GitHub issues |
| UX.86 | 8.1 | Bugs from the GitHub issues |
| UX.87 | 8.1 | Maps |
| UX.88 | 8.1 | Maps |
| UX.89 | 8.1 | Maps |
| UX.90 | 8.1 | Maps |
| UX.91 | 8.1 | Maps |
| UX.93 | 8.5 | Routing |
| UX.95 | 8.1 | Foundations for the issue features |
| UX.96 | 8.1 | Foundations for the issue features |
| UX.97 | 8.1 | Foundations for the issue features |
| UX.98 | 8.1 | Bugs from the GitHub issues |
| UX.99 | 8.1 | Bugs from the GitHub issues |
| UX.103 | 8.1 | Maps |
| UX.104 | 8.1 | Maps |
| UX.106 | 8.1 | Maps |
| UX.111 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| UX.112 | 8.1 | Spherical Sector Coordinates (TOP PRIORITY, Boss 2026-10-15 06:39Z) |
| UX.117 | 8.1 | Maps |
| VIEW.5 | 9.1 | View |
| GEN.281 | 8.1 | Anomalies |
| GEN.236 | 8.1 | Orbital updates |
| PERF.104 | 8.1 | Orbital updates |
| MAP.152 | 8.1 | Fly-through map |
| GEN.99 | 8.1 | Components and wanderers |
| PERF.29 | 8.1 | Queue, caches and performance |
| PERF.30 | 8.1 | Queue, caches and performance |
| API.7 | 8.1 | Admin and API housekeeping |
| MAP.59 | 8.1 | Fly-through map |
| MAP.147 | 8.1 | Fly-through map |
| GEN.24 | 8.1 | Galactic core and halo |
| GEN.210 | 8.1 | Galactic core and halo |
| DOC.4 | 8.1 | Docs |
| GEN.112 | 8.4 | Asteroid fields |
| GEN.9 | later | Multiple galaxies |
| GEN.237 | 8.1 | Components and wanderers |
| GEN.280 | 8.1 | Galactic core and halo |
| ADM.96 | 8.1 | Settings page |
| ADM.87 | 8.1 | Settings page |
| ADM.83 | 8.1 | Settings page |
| MAP.226 | 8.1 | Fly-through map |
| GEN.364 | 8.1 | Galactic core and halo |
| PERF.88 | 8.1 | Queue, caches and performance |
| ADM.84 | 8.1 | Settings page |
| OPS.56 | 8.1 | Reproducible galaxies |
| GEN.102 | 8.1 | Components and wanderers |
| DB.59 | 8.1 | Nebula data foundations |
| TEST.166 | 8.1 | Whole-sphere audits |
| TEST.192 | 8.1 | Bugs and flaky tests |
| GEN.335 | 8.3 | Positions as one object |
| DB.54 | 8.1 | Galactic core and halo |
| ADM.50 | 8.1 | Reproducible galaxies |
| DB.18 | 8.1 | Reproducible galaxies |
| GEN.333 | 8.1 | Galactic core and halo |
| GEN.283 | 8.1 | Galactic core and halo |
| ADM.85 | 8.1 | Settings page |
| UX.116 | 8.1 | Galactic core and halo |
| API.27 | 8.1 | Region query |
| UX.113 | 8.1 | Region query |
| ADM.43 | 8.1 | Settings page |
| ADM.44 | 8.1 | Settings page |
| GEN.135 | 8.1 | Reproducible galaxies |
| TEST.205 | 8.1 | Bugs and flaky tests |
| GEN.286 | 8.1 | Galactic core and halo |
| GEN.287 | 8.1 | Galactic core and halo |
| GEN.306 | 8.1 | Reproducible galaxies |
| GEN.408 | 8.4 | Atmosphere sources and design |
| DB.55 | 8.3 | Positions as one object |
| GEN.336 | 8.3 | Positions as one object |
| GEN.337 | 8.3 | Positions as one object |
| GEN.103 | 8.1 | Components and wanderers |
| GEN.425 | 8.1 | Composition research and design |
| GEN.372 | 8.3 | Positions as one object |
| GEN.368 | 8.1 | Components and wanderers |
| GEN.276 | 8.1 | Region query |
| GEN.431 | 8.1 | Composition research and design |
| GEN.432 | 8.1 | Composition research and design |
| GEN.433 | 8.1 | Composition research and design |
| GEN.434 | 8.1 | Composition research and design |
| GEN.134 | 8.1 | Components and wanderers |
| NAV.60 | 8.1 | Region query |
| GEN.407 | 8.1 | Bugs and flaky tests |
| UX.114 | 8.1 | Region query |
| TEST.217 | 8.1 | Bugs and flaky tests |
| GEN.322 | 8.1 | Components and wanderers |
| GEN.323 | 8.1 | Components and wanderers |
| GEN.324 | 8.1 | Components and wanderers |
| GEN.240 | 8.1 | Components and wanderers |
| TEST.193 | 8.1 | Bugs and flaky tests |
| PERF.97 | 8.1 | Queue, caches and performance |
| NAV.45 | 8.1 | Region query |
| GEN.319 | 8.1 | Nebula and young stars |
| MAP.217 | 8.3 | Movement on the maps |
| GEN.242 | 8.1 | Nebula and young stars |
| GEN.400 | 8.1 | Galaxy laws build |
| GEN.320 | 8.1 | Nebula and young stars |
| GEN.406 | 8.1 | Galaxy laws build |
| DB.65 | 8.3 | Positions as one object |
| MAP.225 | 8.1 | Fly-through map |
| MAP.227 | 8.1 | Fly-through map |
| MAP.228 | 8.1 | Fly-through map |
| DOC.47 | 8.1 | Docs |
| TEST.186 | 8.1 | Whole-sphere audits |
| GEN.526 | 8.1 | Components and wanderers |
| PERF.109 | 8.1 | Queue, caches and performance |
| PERF.110 | 8.1 | Queue, caches and performance |
| PERF.47 | 8.1 | Queue, caches and performance |
| PERF.48 | 8.1 | Queue, caches and performance |
| TEST.156 | 8.1 | Bugs and flaky tests |
| TEST.187 | 8.1 | Bugs and flaky tests |
| TEST.188 | 8.1 | Bugs and flaky tests |
| GEN.401 | 8.1 | Galaxy laws build |
| GEN.472 | 8.1 | Composition tables |
| GEN.409 | 8.1 | Composition research and design |
| GEN.426 | 8.1 | Composition research and design |
| GEN.435 | 8.1 | Composition research and design |
| GEN.436 | 8.1 | Composition research and design |
| GEN.437 | 8.1 | Composition research and design |
| DB.72 | 8.1 | Bugs and flaky tests |
| PERF.105 | 8.1 | Queue, caches and performance |
| PERF.36 | 8.1 | Queue, caches and performance |
| PERF.91 | 8.1 | Queue, caches and performance |
| GEN.308 | 8.1 | Queue, caches and performance |
| TEST.191 | 8.1 | Bugs and flaky tests |
| TEST.221 | 8.1 | Bugs and flaky tests |
| DOC.62 | 8.1 | Bugs and flaky tests |
| DOC.53 | 8.1 | Docs |
| DOC.60 | 8.1 | Docs |
| GEN.527 | 8.1 | Components and wanderers |
| GEN.528 | 8.1 | Components and wanderers |
| DB.53 | 8.1 | Star-system storage |
| GEN.294 | 8.3 | Positions as one object |
| UX.3 | 8.1 | Generate page and directives |
| DOC.24 | 8.1 | Docs |
| UX.118 | 8.1 | Queue, caches and performance |
| GEN.449 | 8.1 | Multi-star and compact systems |
| GEN.448 | 8.1 | Multi-star and compact systems |
| GEN.447 | 8.1 | Multi-star and compact systems |
| MAP.245 | 8.1 | Fly-through map |
| GEN.289 | 8.3 | Positions as one object |
| GEN.463 | 8.1 | Bugs and flaky tests |
| GEN.467 | 8.1 | Bugs and flaky tests |
| GEN.502 | 8.1 | Bugs and flaky tests |
| GEN.529 | 8.1 | Nebula and young stars |
| GEN.259 | 8.1 | Galactic core and halo |
| TEST.213 | 8.1 | Bugs and flaky tests |
| GEN.314 | 8.1 | Multi-star and compact systems |
| TEST.230 | 8.1 | Bugs and flaky tests |
| GEN.238 | 8.1 | Components and wanderers |
| GEN.298 | 8.1 | Components and wanderers |
| GEN.469 | 8.1 | Components and wanderers |
| GEN.321 | 8.1 | Components and wanderers |
| GEN.241 | 8.1 | Components and wanderers |
| TEST.231 | 8.1 | Bugs and flaky tests |
| TEST.232 | 8.1 | Bugs and flaky tests |
| MAP.239 | 8.1 | Fly-through map |
| PERF.20 | 8.1 | Queue, caches and performance |
| PERF.75 | 8.1 | Queue, caches and performance |
| PERF.112 | 8.1 | Queue, caches and performance |
| DOC.52 | 8.1 | Docs |
| DOC.67 | 8.1 | Docs |
| MAP.250 | 8.1 | Fly-through map |
| DOC.68 | 8.1 | Docs |
| PERF.118 | 8.1 | Queue, caches and performance |
| GEN.252 | 8.1 | Multi-star and compact systems |
| GEN.541 | 8.1 | Components and wanderers |
| PERF.135 | 8.1 | Queue, caches and performance |
| PERF.127 | 8.1 | Queue, caches and performance |
| GEN.551 | 8.4 | Bulk composition |
| GEN.243 | 8.1 | Nebula and young stars |
| GEN.438 | 8.1 | Composition tables |
| GEN.532 | 8.1 | Galaxy laws build |
| GEN.533 | 8.1 | Galaxy laws build |
| GEN.567 | 8.1 | Composition research and design |
| MAP.243 | 8.1 | Interface update |
| DOC.69 | 8.1 | Docs |
| OPS.59 | 8.1 | Queue, caches and performance |
| UX.110 | 8.1 | Description layer |
| TEST.134 | 8.1 | Bugs and flaky tests |
| TEST.170 | 8.1 | Bugs and flaky tests |
| TEST.171 | 8.1 | Bugs and flaky tests |
| TEST.172 | 8.1 | Bugs and flaky tests |
| TEST.174 | 8.1 | Bugs and flaky tests |
| UX.133 | 8.1 | Bugs and flaky tests |
| TEST.189 | 8.1 | Bugs and flaky tests |
| TEST.216 | 8.1 | Bugs and flaky tests |
| TEST.203 | 8.1 | Bugs and flaky tests |
| TEST.199 | 8.1 | Bugs and flaky tests |
| GEN.397 | 8.1 | Bugs and flaky tests |
| TEST.208 | 8.1 | Bugs and flaky tests |
| GEN.218 | 8.1 | Description layer |
| GEN.219 | 8.1 | Description layer |
| DB.50 | 8.1 | Description layer |
| DB.51 | 8.1 | Description layer |
| UX.119 | 8.1 | Description layer |
| UX.120 | 8.1 | Description layer |
| UX.105 | 8.1 | Description layer |
| GEN.220 | 8.1 | Description layer |
| UX.108 | 8.1 | Description layer |
| GEN.221 | 8.1 | Description layer |
| DB.52 | 8.1 | Description layer |
| TEST.155 | 8.1 | Description layer |
| GEN.439 | 8.1 | Composition tables |
| GEN.440 | 8.1 | Star composition |
| PERF.126 | 8.1 | Queue, caches and performance |
| GEN.552 | 8.1 | Multi-star and compact systems |
| TEST.220 | 8.1 | Bugs and flaky tests |
| TEST.235 | 8.1 | Bugs and flaky tests |
| GEN.562 | 8.1 | Multi-star and compact systems |
| GEN.402 | 8.1 | Galaxy laws build |
| OPS.60 | 8.1 | Queue, caches and performance |
| OPS.61 | 8.1 | Queue, caches and performance |
| TEST.236 | 8.1 | Bugs and flaky tests |
| TEST.225 | 8.1 | Bugs and flaky tests |
| TEST.245 | 8.1 | Bugs and flaky tests |
| MAP.247 | 8.1 | Fly-through map |
| UX.109 | 8.1 | Description layer |
| GEN.548 | 8.1 | Composition tables |
| GEN.543 | 8.1 | Multi-star and compact systems |
| GEN.561 | 8.1 | Description layer |
| UX.107 | 8.1 | Description layer |
| DOC.45 | 8.1 | Description layer |
| GEN.572 | 8.1 | Composition research and design |
| GEN.549 | 8.1 | Galaxy laws build |
| TEST.223 | 8.1 | Bugs and flaky tests |
| TEST.238 | 8.1 | Bugs and flaky tests |
| MAP.252 | 8.1 | Fly-through map |
| GEN.244 | 8.1 | Components and wanderers |
| GEN.282 | 8.1 | Galactic core and halo |
| GEN.563 | 8.1 | Galactic core and halo |
| TEST.248 | 8.1 | Bugs and flaky tests |
| GEN.317 | 8.1 | Multi-star and compact systems |
| GEN.547 | 8.1 | Galaxy laws build |
| GEN.546 | 8.1 | Galaxy laws build |
| TEST.250 | 8.1 | Bugs and flaky tests |
| GEN.404 | 8.1 | Galaxy laws build |
| TEST.247 | 8.1 | Bugs and flaky tests |
| TEST.249 | 8.1 | Bugs and flaky tests |
| TEST.251 | 8.1 | Bugs and flaky tests |
| PERF.130 | 8.2 | Robust tile caching |
| PERF.106 | 8.1 | Queue, caches and performance |
| GEN.590 | 8.1 | Multi-star and compact systems |
| TEST.246 | 8.1 | Bugs and flaky tests |
| TEST.252 | 8.1 | Bugs and flaky tests |
| TEST.224 | 8.1 | Bugs and flaky tests |
| TEST.227 | 8.1 | Bugs and flaky tests |
| TEST.243 | 8.1 | Bugs and flaky tests |
| PERF.138 | 8.1 | Multi-star and compact systems |
| PERF.132 | 8.2 | Robust tile caching |
| TEST.129 | 8.1 | Bugs and flaky tests |
| PERF.107 | 8.1 | Queue, caches and performance |
| PERF.108 | 8.1 | Queue, caches and performance |
| GEN.591 | 8.1 | Multi-star and compact systems |
| GEN.318 | 8.1 | Multi-star and compact systems |
| GEN.130 | 8.1 | Multi-star and compact systems |
| PERF.133 | 8.2 | Robust tile caching |
| PERF.139 | 8.2 | Robust tile caching |
| GEN.544 | 8.1 | Galaxy laws build |
| TEST.256 | 8.1 | Bugs and flaky tests |
| GEN.592 | 8.1 | Multi-star and compact systems |
| TEST.253 | 8.1 | Bugs and flaky tests |
| TEST.255 | 8.1 | Bugs and flaky tests |
| TEST.259 | 8.1 | Bugs and flaky tests |
| PERF.143 | 8.1 | Queue, caches and performance |
| PERF.142 | 8.1 | Queue, caches and performance |
| DB.79 | 8.2 | Robust tile caching |
| PERF.144 | 8.1 | Queue, caches and performance |
| GEN.378 | 8.1 | Galaxy types groundwork |
| PERF.140 | 8.2 | Robust tile caching |
| MAP.251 | 8.1 | Fly-through map |
| GEN.457 | 8.1 | Composition research and design |
| GEN.530 | 8.1 | Components and wanderers |
| DB.82 | 8.2 | Robust tile caching |
| DB.31 | 8.1 | Description layer |
| GEN.601 | 8.1 | Admin control suite |
| GEN.608 | 8.1 | Admin control suite |
| PERF.148 | 8.1 | Queue, caches and performance |
| GEN.597 | 8.1 | Galaxy laws build |
| GEN.602 | 8.1 | Admin control suite |
| GEN.603 | 8.1 | Admin control suite |
| GEN.604 | 8.1 | Admin control suite |
| GEN.605 | 8.1 | Admin control suite |
| GEN.606 | 8.1 | Admin control suite |
| GEN.609 | 8.1 | Admin control suite |
| GEN.610 | 8.1 | Admin control suite |
| DB.81 | 8.1 | Admin control suite |
| GEN.611 | 8.1 | Admin control suite |
| DB.80 | 8.1 | Bugs and flaky tests |
| UX.141 | 8.2 | Robust tile caching |
| ADM.98 | 8.1 | Admin control suite |
| ADM.101 | 8.1 | Admin control suite |
| PERF.141 | 8.2 | Robust tile caching |
| PERF.145 | 8.1 | Queue, caches and performance |
| TEST.263 | 8.1 | Bugs and flaky tests |
| TEST.257 | 8.1 | Bugs and flaky tests |
| UX.145 | 8.1 | Admin control suite |
| GEN.621 | 8.1 | Admin control suite |
| DB.44 | 8.1 | Reproducible galaxies |
| TEST.260 | 8.1 | Bugs and flaky tests |
| TEST.264 | 8.1 | Bugs and flaky tests |
| GEN.596 | 8.1 | Galaxy laws build |
| MAP.172 | 8.1 | Nebula data foundations |
| GEN.545 | 8.1 | Galaxy laws build |
| DB.33 | 8.1 | Bugs and flaky tests |
| DB.66 | 8.3 | Galactic gravity |
| TEST.265 | 8.1 | Bugs and flaky tests |
| GEN.607 | 8.1 | Galaxy laws build |
| TEST.241 | 8.1 | Multi-star and compact systems |
| ADM.102 | 8.1 | Admin control suite |
| ADM.103 | 8.1 | Admin control suite |
| ADM.104 | 8.1 | Admin control suite |
| ADM.97 | 8.1 | Admin control suite |
| ADM.105 | 8.1 | Admin control suite |
| ADM.106 | 8.1 | Admin control suite |
| ADM.100 | 8.1 | Admin control suite |
| GEN.370 | 8.3 | Galactic gravity |
| TEST.258 | 8.1 | Bugs and flaky tests |
| GEN.600 | 8.1 | Queue, caches and performance |
| TEST.266 | 8.1 | Bugs and flaky tests |
| MAP.139 | 8.1 | Fly-through map |
| GEN.617 | 8.1 | Admin control suite |
| GEN.374 | 8.3 | Positions as one object |
| GEN.140 | 8.3 | Positions as one object |
| GEN.144 | 8.3 | Galactic gravity |
| GEN.334 | 8.3 | Orbital update core |
| GEN.341 | 8.3 | Orbital update core |
| GEN.345 | 8.3 | Orbital update core |
| GEN.568 | 8.4 | Bulk composition |
| GEN.598 | 8.1 | Galaxy laws build |
| GEN.311 | 8.1 | Multi-star and compact systems |
| GEN.269 | 8.1 | Multi-star and compact systems |
| PERF.121 | 8.1 | Queue, caches and performance |
| PERF.122 | 8.1 | Queue, caches and performance |
| TEST.267 | 8.1 | Bugs and flaky tests |
| DB.76 | 8.4 | Bulk composition |
| DB.43 | 8.1 | Reproducible galaxies |
| GEN.569 | 8.4 | Bulk composition |
| UX.139 | 8.4 | Bulk composition |
| GEN.573 | 8.4 | Bulk composition |
| GEN.299 | 8.3 | The Quantity class |
| TEST.274 | 8.1 | Bugs and flaky tests |
| TEST.269 | 8.1 | Bugs and flaky tests |
| ADM.118 | 8.1 | Admin control suite |
| ADM.92 | 8.1 | Admin control suite |
| ADM.119 | 8.1 | Admin control suite |
| ADM.120 | 8.1 | Admin control suite |
| ADM.122 | 8.1 | Admin control suite |
| ADM.123 | 8.1 | Admin control suite |
| ADM.107 | 8.1 | Admin control suite |
| ADM.108 | 8.1 | Admin control suite |
| ADM.109 | 8.1 | Admin control suite |
| ADM.124 | 8.1 | Admin control suite |
| ADM.113 | 8.1 | Admin control suite |
| ADM.115 | 8.1 | Admin control suite |
| ADM.121 | 8.1 | Admin control suite |
| ADM.93 | 8.1 | Admin control suite |
| ADM.55 | 8.1 | Admin control suite |
| ADM.110 | 8.1 | Admin control suite |
| ADM.111 | 8.1 | Admin control suite |
| ADM.112 | 8.1 | Admin control suite |
| ADM.114 | 8.1 | Admin control suite |
| ADM.116 | 8.1 | Admin control suite |
| ADM.117 | 8.1 | Admin control suite |
| ADM.56 | 8.1 | Admin control suite |
| ADM.99 | 8.1 | Admin control suite |
| TEST.268 | 8.1 | Bugs and flaky tests |
| DB.56 | 8.3 | Orbital update core |
| TEST.184 | 8.1 | Whole-sphere audits |
| TEST.185 | 8.1 | The headline closes |
| TEST.159 | 8.1 | The headline closes |
| TEST.160 | 8.1 | The headline closes |
| GEN.410 | 8.1 | Composition tables |
| TEST.278 | 8.1 | Bugs and flaky tests |
| MAP.254 | 8.2 | Viewers and object pages |
| DB.30 | 8.1 | Tech levels and species |
| DB.49 | 8.1 | Reproducible galaxies |
| DB.16 | 8.1 | Reproducible galaxies |
| DB.36 | 8.6 | Change plans |
| DB.37 | 8.6 | Aspects |
| DB.38 | 8.6 | Regions |
| DB.60 | 8.6 | Regions |
| DB.39 | 10.0 | Request tables |
| DB.61 | 8.8 | Clusters and satellites on the map |
| DB.64 | 8.9 | Sky data |
| DB.63 | 9.3 | Ships and stations |
| GEN.379 | 8.3 | Galactic gravity |
| DB.77 | 8.4 | Bulk composition |
| DB.58 | 8.2 | Accounts and sign-up |
| DB.62 | 8.9 | Upload groundwork |
| GEN.348 | 8.4 | The Quantity call-site sweep |
| TEST.279 | 8.1 | Bugs and flaky tests |
| DB.87 | 8.2 | Viewers and object pages |
| GEN.585 | 8.4 | Bulk composition |
| GEN.586 | 8.4 | Bulk composition |
| UX.129 | 8.2 | Bookmarks |
| DB.88 | 8.2 | Viewers and object pages |
| TEST.282 | 8.1 | Bugs and flaky tests |
| NAV.17 | 8.5 | Saved courses |
| PERF.146 | 8.1 | The headline closes |
| GEN.599 | 8.1 | Queue, caches and performance |
| PERF.153 | 8.1 | Queue, caches and performance |
| PERF.152 | 8.1 | Queue, caches and performance |
| PERF.151 | 8.1 | Queue, caches and performance |
| GEN.416 | 9.0 | Atmosphere colour |
| TEST.283 | 8.1 | The headline closes |
| UX.140 | 8.4 | Bulk composition |
| TEST.270 | 8.1 | Bugs and flaky tests |
| TEST.272 | 8.1 | Bugs and flaky tests |
| TEST.280 | 8.1 | Bugs and flaky tests |
| TEST.271 | 8.1 | Bugs and flaky tests |
| PERF.163 | 8.1 | The headline closes |
| GEN.304 | 8.1 | Reproducible galaxies |
| GEN.427 | 8.1 | Composition tables |
| DB.85 | 8.2 | Viewers and object pages |
| GEN.498 | 8.1 | Composition tables |
| GEN.496 | 8.1 | Composition tables |
| PERF.150 | 8.1 | Queue, caches and performance |
| GEN.470 | 8.1 | Composition tables |
| GEN.478 | 8.4 | Atmosphere tables |
| GEN.307 | 8.1 | Reproducible galaxies |
| GEN.623 | 8.1 | Composition tables |
| GEN.495 | 8.4 | Surface and small bodies |
| GEN.626 | 8.1 | Reproducible galaxies |
| TEST.277 | 8.1 | Bugs and flaky tests |
| TEST.275 | 8.1 | Bugs and flaky tests |
| TEST.285 | 8.1 | Bugs and flaky tests |
| ADM.127 | 8.1 | Admin control suite |
| UX.146 | 8.1 | Admin control suite |
| TEST.286 | 8.1 | Bugs and flaky tests |
| GEN.305 | 8.1 | Reproducible galaxies |
| GEN.136 | 8.1 | Reproducible galaxies |
| TEST.281 | 8.1 | Bugs and flaky tests |
| GEN.628 | 8.1 | Reproducible galaxies |
| MAP.146 | 8.1 | Fly-through map |
| PERF.164 | 8.1 | The headline closes |
| GEN.325 | 8.1 | Globular clusters |
| TEST.288 | 8.1 | Bugs and flaky tests |
| GEN.625 | 8.1 | Reproducible galaxies |
| TEST.181 | 8.1 | Reproducible galaxies |
| OPS.57 | 8.1 | Reproducible galaxies |
| PERF.168 | 8.1 | The headline closes |
| GEN.329 | 8.1 | Multi-star and compact systems |
| PERF.165 | 8.1 | The headline closes |
| GEN.627 | 8.1 | Reproducible galaxies |
| GEN.369 | 8.1 | Components and wanderers |
| GEN.624 | 8.1 | Components and wanderers |
| GEN.403 | 8.1 | Galaxy laws build |
| GEN.284 | 8.1 | Galactic core and halo |
| TEST.196 | 8.1 | The headline closes |
| GEN.468 | 8.1 | Composition research and design |
| GEN.465 | 8.1 | Composition research and design |
| ADM.125 | 8.1 | Admin control suite |
| PERF.169 | 8.1 | The headline closes |
| POP.14 | 8.1 | Tech levels and species |
| GEN.636 | 8.1 | The headline closes |
| OPS.58 | 8.1 | Reproducible galaxies |
| OPS.28 | 8.1 | Reproducible galaxies |
| GEN.458 | 8.4 | Bulk composition |
| API.32 | 8.1 | Bugs and flaky tests |
| TEST.293 | 8.1 | Bugs and flaky tests |
| TEST.294 | 8.1 | Bugs and flaky tests |
| SEC.32 | 8.2 | Accounts and sign-up |
| TEST.291 | 8.1 | Bugs and flaky tests |
| GEN.497 | 8.1 | Composition tables |
| GEN.258 | 8.1 | Nebula and young stars |
| ADM.126 | 8.1 | Admin control suite |
| PERF.167 | 8.1 | The headline closes |
| GEN.253 | 8.1 | Multi-star and compact systems |
| GEN.146 | 8.4 | Planet classes |
| TEST.296 | 8.1 | Bugs and flaky tests |
| GEN.630 | 8.1 | Multi-star and compact systems |
| MAP.174 | 8.1 | Nebula data foundations |
| GEN.638 | 8.1 | The headline closes |
| PERF.170 | 8.1 | The headline closes |
| GEN.267 | 8.1 | Multi-star and compact systems |
| UX.124 | 8.1 | Admin control suite |
| OPS.63 | 8.1 | Bugs and flaky tests |
| GEN.417 | 8.1 | Radiation chain for habitability (GEN.445) |
| GEN.453 | 8.1 | Radiation chain for habitability (GEN.445) |
| UX.147 | 8.1 | Generate page and directives |
| PERF.160 | 8.1 | The headline closes |
| PERF.158 | 8.1 | The headline closes |
| PERF.157 | 8.1 | The headline closes |
| PERF.156 | 8.1 | The headline closes |
| UX.101 | 8.1 | Habitability system |
| GEN.622 | 8.1 | Bulk composition |
| GEN.629 | 8.1 | Surface and small bodies |
| PERF.159 | 8.1 | The headline closes |
| GEN.503 | 8.1 | Tuned draws and selection models |
| GEN.504 | 8.1 | Tuned draws and selection models |
| GEN.505 | 8.1 | Tuned draws and selection models |
| GEN.506 | 8.1 | Tuned draws and selection models |
| GEN.507 | 8.1 | Tuned draws and selection models |
| GEN.508 | 8.1 | Tuned draws and selection models |
| GEN.511 | 8.1 | Tuned draws and selection models |
| GEN.513 | 8.1 | Tuned draws and selection models |
| GEN.514 | 8.1 | Tuned draws and selection models |
| GEN.516 | 8.1 | Tuned draws and selection models |
| GEN.517 | 8.1 | Tuned draws and selection models |
| GEN.525 | 8.1 | Tuned draws and selection models |
| GEN.413 | 8.1 | Radiation chain for habitability (GEN.445) |
| PERF.161 | 8.1 | The headline closes |
| PERF.171 | 8.1 | The headline closes |
| UX.151 | 8.1 | Habitability system |
| UX.149 | 8.1 | Generate page and directives |
| GEN.450 | 8.1 | Components and wanderers |
| GEN.418 | 8.1 | Atmosphere colour |
| GEN.419 | 8.1 | Atmosphere colour |
| TEST.309 | 8.1 | Bugs and flaky tests |
| MAP.271 | 8.1 | Fly-through map |
| MAP.272 | 8.1 | Fly-through map |
| MAP.273 | 8.1 | Fly-through map |
| GEN.420 | 8.1 | Descriptions and maps |
| PERF.115 | 8.1 | Queue, caches and performance |
| PERF.116 | 8.1 | Queue, caches and performance |
| PERF.117 | 8.1 | Queue, caches and performance |
| MAP.276 | 8.1 | Fly-through map |
| OPS.68 | 8.1 | Bugs and flaky tests |
| MAP.278 | 8.1 | Fly-through map |
| GEN.566 | 8.1 | Bulk composition |
| GEN.580 | 8.1 | Bulk composition |
| GEN.579 | 8.1 | Bulk composition |
| DOC.70 | 8.1 | Bulk composition |
| GEN.641 | 8.1 | Bulk composition |
| GEN.659 | 8.1 | Bulk composition |
| DB.86 | 8.1 | Bulk composition |
| UX.128 | 8.1 | The Quantity call-site sweep |
| MAP.244 | 8.1 | Descriptions and maps |
| DOC.64 | 8.1 | Descriptions and maps |
| TEST.312 | 8.1 | Bugs and flaky tests |
| UX.100 | 8.1 | Habitability system |
| UX.155 | 8.1 | Habitability system |
| UX.156 | 8.1 | Habitability system |
| DOC.78 | 8.1 | Habitability system |
| GEN.655 | 8.1 | Habitability system |
| DB.93 | 8.1 | Habitability system |
| POP.15 | 8.1 | Habitability system |
| GEN.665 | 8.1 | Tuned draws and selection models |
| GEN.664 | 8.1 | Tuned draws and selection models |
| DOC.74 | 8.1 | Mixture and stored values |
| GEN.451 | 8.1 | Mixture and stored values |
| GEN.454 | 8.1 | Mixture and stored values |
| GEN.499 | 8.1 | Haze, dust and cloud |
| GEN.456 | 8.1 | Haze, dust and cloud |
| GEN.501 | 8.1 | Weight research (atmosphere) |
| GEN.663 | 8.1 | Tuned draws and selection models |
| GEN.92 | 8.4 | Planet classes |
| TEST.310 | 8.1 | Bugs and flaky tests |
| UX.165 | 8.1 | Consistency |
| GEN.680 | 8.1 | Consistency |
| GEN.679 | 8.1 | Consistency |
| UX.148 | 8.1 | Generate page and directives |
| MAP.270 | 8.1 | Fly-through map |
| GEN.656 | 8.1 | Habitability system |
| UX.154 | 8.1 | Habitability system |
| GEN.518 | 8.1 | Tuned draws and selection models |
| GEN.520 | 8.1 | Tuned draws and selection models |
| PERF.113 | 8.1 | Queue, caches and performance |
| PERF.114 | 8.1 | Queue, caches and performance |
| PERF.177 | 8.1 | The headline closes |
| TEST.316 | 8.1 | Bugs and flaky tests |
| MAP.256 | 8.1 | Fly-through map |
| GEN.414 | 8.1 | Radiation chain for habitability (GEN.445) |
| DB.68 | 8.1 | Consistency |
| MAP.143 | 8.1 | Habitability system |
| GEN.694 | 8.1 | Mixture and stored values |
| GEN.693 | 8.1 | Mixture and stored values |
| DB.73 | 8.1 | Radiation chain for habitability (GEN.445) |
| PERF.123 | 8.1 | Queue, caches and performance |
| UX.182 | 8.1 | Habitability system |
| UX.181 | 8.1 | Habitability system |
| UX.180 | 8.1 | Habitability system |
| GEN.459 | 8.1 | Radiation chain for habitability (GEN.445) |
| UX.159 | 8.1 | Admin control suite |
| GEN.412 | 8.1 | Radiation chain for habitability (GEN.445) |
| MAP.288 | 8.1 | Fly-through map |
| PERF.176 | 8.1 | The headline closes |
| UX.162 | 8.1 | The headline closes |
| GEN.270 | 8.1 | Multi-star and compact systems |
| TEST.322 | 8.1 | Bugs and flaky tests |
| TEST.325 | 8.1 | Bugs and flaky tests |
| DB.94 | 8.1 | Radiation chain for habitability (GEN.445) |
| GEN.671 | 8.1 | Habitability system |
| UX.163 | 8.1 | Habitability system |
| UX.160 | 8.1 | Habitability system |
| UX.161 | 8.1 | Habitability system |
| TEST.306 | 8.1 | Bugs and flaky tests |
| MAP.219 | 8.1 | Habitability system |
| MAP.263 | 8.1 | Fly-through map |
| UX.157 | 8.1 | Habitability system |
| UX.166 | 8.1 | Habitability system |
| TEST.327 | 8.1 | Bugs and flaky tests |
| GEN.574 | 8.1 | Bulk composition |
| TEST.273 | 8.1 | Bugs and flaky tests |
| TEST.233 | 8.1 | Bugs and flaky tests |
| TEST.218 | 8.1 | Bugs and flaky tests |
| GEN.666 | 8.1 | Habitability system |
| GEN.692 | 8.1 | Radiation chain for habitability (GEN.445) |
| MAP.264 | 8.1 | Fly-through map |
| DB.97 | 8.1 | Reproducible galaxies |
| UX.167 | 8.1 | Habitability system |
| GEN.428 | 8.1 | Radiation chain for habitability (GEN.445) |
| GEN.429 | 8.1 | Radiation chain for habitability (GEN.445) |
| DB.69 | 8.1 | Bulk composition |
| GEN.443 | 8.1 | Radiation chain for habitability (GEN.445) |
| MAP.258 | 8.1 | Fly-through map |
| GEN.631 | 8.1 | Multi-star and compact systems |
| GEN.643 | 8.1 | Habitability system |
| TEST.313 | 8.1 | Bugs and flaky tests |
| GEN.702 | 8.1 | Habitability system |
| PERF.155 | 8.1 | The headline closes |
| GEN.560 | 8.1 | Composition research and design |
| OPS.65 | 8.1 | Bugs and flaky tests |
| GEN.271 | 8.1 | Multi-star and compact systems |
