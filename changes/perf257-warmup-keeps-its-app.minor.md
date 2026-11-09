**Headline:** An unchanged warm-up now takes minutes instead of 28: the workers keep the app and a node stage skips the nodes that are already stored.

### Changed
- **The warm-up app is built once per worker (PERF.257).** The five warm commands whose units visit pages (`warm-galaxy-nodes`, `warm-galaxy-deep-nodes`, `warm-sectors`, `warm-systems`, `warm-phenomena`) build the app in the worker before the work horse is forked, so no unit pays the 7 s the app took to import. Measured: 7.1 s a horse before, 0.0 s after.
- **A node stage opens with a stamp check.** The nodes the disk cache already holds in the current generation are left out of the units (`tilecache.stored_nodes`); an unchanged rerun is one empty unit that reports them as already stored, and a changed sector cuts only its own nodes. A failed check cuts every node as before.
- Generated output is unchanged.

generation-affects: none
