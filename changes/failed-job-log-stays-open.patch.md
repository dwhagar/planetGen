### Fixed
- A failed Generate-page job no longer reloads away from its error. The log stays open with the error in view until you click **Continue** (ADM.24). At an interactive terminal a failed `planetgen` run waits for Enter before it exits, so the output isn't lost with the console window; runs with redirected input or output exit at once as before.
