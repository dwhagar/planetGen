**Headline:** The train PR's body now shows how an exact Earth comes out of the generator against a stored baseline, so the items that move the air, water, life or class can quote a before and after.

### Changed

- TEST.726 (stacked draft; proposed, to be confirmed by TODO before it merges): the Earth analog of the Class M audit is now `scripts/earth_analog.py` (`run`, `baseline`, `compare`), with its baseline of 800 paired trials in `scripts/data/earth_analog_baseline.json` (made on main at the commit it names). `python scripts/train.py pr-body` prints a "## Earth analog" section: the classes, hydrospheres, temperature bins, medians and failing Class M conditions of the first 400 trials against the baseline, and the trials whose class changed. The draws are paired by seed, so a train that changes nothing shows no movement. It is for information and never blocks a train (any failure is one sentence in the body; `--no-earth-analog` leaves it out). The analysis-only measurement hooks of the audit copy (scaling the volcanic nitrogen, the carbon dioxide gate, the wet-bulb edge) are not carried over.

generation-output: unchanged
generation-affects: none
