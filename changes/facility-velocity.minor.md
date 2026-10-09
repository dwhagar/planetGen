### Added
- A stand-alone facility stores a galactic velocity (`facilities.velocity_x_kms`, `_y_kms`, `_z_kms`, schema v64): the rotation curve's tangent at its place, filled when it is added and turned with its position at each orbit update, so every object in space carries a vector. Run `sudo ./update.sh` to migrate; existing stand-alone facilities get theirs.
