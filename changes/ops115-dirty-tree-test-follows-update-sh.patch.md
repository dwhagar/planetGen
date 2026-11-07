**Headline:** The test of update.sh's dirty-tree check names the git command it now runs.

### Changed

- `src/tests/test_deploy_scripts.py` (OPS.115 follow-up): `test_the_dirty_tree_warning_leaves_out_the_state_folder_the_update_owns` looked for the old text of the check, `status --porcelain -- . ':(exclude)state'`. OPS.115 (PR 5562) added `--untracked-files=no` to that command, so the test went red on main; it now looks for the new text. The check itself is unchanged.

generation-output: unchanged
