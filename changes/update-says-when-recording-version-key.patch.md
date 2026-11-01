**Headline:** The update says what it is doing while it records the version key, so the quiet stretch after the migration no longer looks hung.

### Fixed

After the databases were migrated, `update.sh` printed nothing while it recorded the version key, which runs the generator's fixed test cases once to fingerprint the code (about 5 seconds, longer on a slow host). It now prints a line saying so before it starts. No generated output changes.
