**Headline:** Every example plist, XML and systemd unit file in the repository is checked by a test, so a malformed one cannot ship unnoticed.

### Changed
- `test_example_files_lint.py` parses every tracked `.plist` (`plistlib`) and `.xml` file, and runs `systemd-analyze verify` on every tracked `.service`, `.timer`, `.socket`, `.path`, `.mount` and `.target` (skipped where `systemd-analyze` is missing; a "command not installed" notice for the server paths is ignored, anything else fails). New examples are found with `git ls-files`. It replaces the macOS-only plist check.

TODO items resolved: OPS.31 / Added: none
