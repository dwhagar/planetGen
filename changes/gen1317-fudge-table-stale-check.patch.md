**Headline:** The fudge table now says what it was built from, and `update.sh` says when it is out of date, without ever running the matrix.

### Added

- The header of `fudge_table.json` holds `inputs`: the SHA-256 of the in-scope stand-in rows, of the class table and of the matrix, and the generator epoch the matrix was measured at (stamped in the matrix as `sample.epoch`; the table holds no epoch of its own). `scripts/build_fudge_table.py --check-inputs` compares them with the files; `update.sh` step 3 runs it after the References pages. Silent when current; when a hash differs it names the file and the build command; it also says how many generator epochs, and how many that change the fill stage, have passed since the matrix was measured. It only warns, and the hours-long matrix run stays on a lane machine (GEN.1317).

generation-output: unchanged
