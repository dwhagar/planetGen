### Fixed

- Log lines printed while a progress bar is up are no longer hard-wrapped at the console's width (80 columns in a web job's log), and text in square brackets is printed as written instead of being read as `rich` markup, which could drop or crash a line (ADM.23).
