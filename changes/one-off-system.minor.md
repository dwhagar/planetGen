### Added
- **One-off star systems from the admin site.** A new page,
  `/admin/generate/system` (linked from Generate), offers every
  `generate.py system` option, including a pasted system file and the
  debug narration, and shows the result as Markdown or wikitext with Copy,
  Download and a rendered preview. Nothing is saved to the database.
- **`generate.py system --output FILE`** (`-o`, `-` for stdout) writes
  the system's page instead of saving the system to the database.
