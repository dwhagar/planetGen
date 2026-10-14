### Fixed
- No internal TODO code (such as GEN.96 or PERF.56) shows in the command-line help, the admin Stats page hint or any other text a user reads, and the help no longer names an internal constant. A test now scans the templates, static files, the strings the code shows and every command's `--help`, and fails if one comes back.
