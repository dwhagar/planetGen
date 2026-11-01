**Headline:** The update shows a progress bar while it records the version key, so the quiet stretch after the migration is no longer a blank wait.

### Changed

Recording the version key runs the generator's fixed test cases once to fingerprint the code, and the update printed nothing while it did. It now shows a bar with the cases done and the time left on a terminal, and one line per case without one. No generated output changes.
