# Generated pages

Every page in this folder is built by a script (the plan pages and the code map are HTML; the two code inventories are Markdown, so the documentation website carries them too). Do not edit them by hand: change the source and run the script.

| Page | Built by | From |
|---|---|---|
| [architecture-map.html](architecture-map.html) | `scripts/build_architecture_map.py` | the code under `src/planetgen/` |
| [execution-plan.html](execution-plan.html) | `scripts/build_todo_docs.py` | [../plan/lanes.md](../plan/lanes.md) and [../plan/execution-notes.md](../plan/execution-notes.md) |
| [release-roadmap.html](release-roadmap.html) | `scripts/build_todo_docs.py` | [../plan/release-roadmap.md](../plan/release-roadmap.md) and [../TODO.md](../TODO.md) |
| [todo-reference.html](todo-reference.html) | `scripts/build_todo_docs.py` | [../TODO.md](../TODO.md) and [../plan/todo-number-map.md](../plan/todo-number-map.md) |
| [todo-tree.html](todo-tree.html) | `scripts/build_todo_docs.py` | [../TODO.md](../TODO.md) |
| [version-plan.html](version-plan.html) | `scripts/build_todo_docs.py` | [../TODO.md](../TODO.md) and the files in [../plan/releases/](../plan/releases/) |
| [database-columns.md](database-columns.md) | `scripts/build_reference_inventory.py` | `src/planetgen/db/schema.sql`: every table and column with its comment |
| [api-routes.md](api-routes.md) | `scripts/build_reference_inventory.py` | the route decorators under `src/planetgen/api/`: every HTTP route |

`docs/index.html` (one level up, because GitHub Pages serves it from there) is built by `scripts/build_docs_index.py`. `python scripts/build_global_docs.py` rebuilds all of them; `--check` fails when one is out of date.

`docs/site/` (the documentation website: every Markdown document as a web page, with a menu, search and rendered formulas) is built by `scripts/build_docs_site.py`, which `scripts/build_global_docs.py` also runs. To serve it, set GitHub Pages to deploy from the `main` branch and the `/docs` folder; the entry points are `docs/index.html` and `docs/site/index.html`.
