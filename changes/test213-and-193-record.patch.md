**Headline:** The TEST.193 investigation record and a new loaded-browser flake are filed.

### Added

- TEST.213 (`test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star` fails once in two loaded runs: the wheel loop opens the system before the disc reaches full share), found by Core 1, for Bug Hunt.

### Changed

- TEST.193 carries Core 1's reproduction and fix (a request thread still running when the module fixture drops the database; a bounded drop that names the holders), and says TEST.195's stall was not reproduced and its fetch was removed.

TODO items resolved: none / Added: TEST.213
