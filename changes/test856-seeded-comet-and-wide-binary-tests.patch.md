**Headline:** The comet test and the wide-binary tests draw from fixed streams, so they take the same time and the same systems in every run.

### Fixed

- TEST.856: `test_mass_object.py::test_a_comet_advances_on_its_hosts_pull` drew unseeded systems until one had a comet, so it took 25 to 46 minutes by luck (and one run passed with a different system than the next). `_first` now draws each attempt from the stream `mass-object:<name>:<n>`, and the comet test forces the comet cloud on and runs three cases (three named streams) instead of one: the whole file takes 1 min 20 s. The loops of `test_wide_binary.py` (the prune trials, the bounds trials, the star-limit trials, the 60 systems of `test_wide_binaries_never_overlap`) draw from fixed streams the way the age and ceiling tests there already do, and the 60-system test is six chunks of ten, so each is under the gate's 1,500 s cap (the single test took 30 minutes) and the chunks run side by side. The same 60 systems are checked; all 43 tests of the file pass in 38 min 42 s with three workers.

generation-output: unchanged
