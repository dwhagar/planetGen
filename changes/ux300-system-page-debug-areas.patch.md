**Headline:** The one-off system page can now pick which debug areas its debug log shows.

### Added
- **A "Debug areas" box on the one-off system page (UX.300).** `/admin/generate/system` takes the same area list `--debug-areas` takes (names separated by commas, `all`, or a leading `-` to leave one out); it is checked the way the command line checks it and passed on when the debug log is ticked.
