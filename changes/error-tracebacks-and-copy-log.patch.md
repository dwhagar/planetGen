### Fixed
- An unexpected error, a database error or a file error in a running action now prints its full traceback to the console and into the job's web log (and the debug log), not only a one-line message. A worker's own traceback comes back with its error. The job panel has a **Copy log** button that puts the whole output on the clipboard (ADM.25).
