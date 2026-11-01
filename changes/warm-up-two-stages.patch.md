**Headline:** update.sh now builds only what the Apache reload and the health check need before the reload, then runs the full cache warm-up afterwards and stays up, with a progress bar, until it finishes.

### Changed
- Stage 1 runs `planetgen warm opening --progress` in the foreground before the reload (the cluster field, the object index, the opening views and the first page of each list), so the reload can no longer interrupt a warm-up that is still running. After the reload and the /api/health check (at config.json's api_base_url only) the update releases its guard, so the site opens, and stage 2 runs the full warm-up with a progress bar; update.sh ends when it does.
- `planetgen warm` has a new `--progress` option (a bar after each family) and no longer has `--opening-done`.

generation-output: unchanged
TODO items resolved: none (new item, filed by the TODO thread)
