**Headline:** The updater and the data downloader now print the URL of every file they load.

### Changed
- **`planetgen.cli.datasets fetch` prints each file's URL (GEN.600).** Under every dataset's line it prints one line per file with its https URL and what happened to it: `fetched` (requested and kept), `kept` (already there with the right size and sha256, so not requested again) or `failed` with the reason. `update.sh` and `install.sh` run it, so their output shows the same lines.
- **`planetgen.cli.datasets status` prints the URL of every missing file,** with why it counts as missing (not downloaded, or the size or sha256 does not match), and nothing for a dataset that is present.
