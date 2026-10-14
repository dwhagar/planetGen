### Changed
- The bright and massive star scatters work out each ring's population densities as one NumPy array instead of one call per angle bin (PERF.89). The numbers agree with the one-point form to 1e-12 relative (NumPy's exp, log and cos differ from the maths library's in the last digit or two), so a progress-bar estimate or a star's bin weight can move in the 16th digit.
