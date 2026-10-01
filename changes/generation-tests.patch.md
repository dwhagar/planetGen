### Added
- **Generation tests (TEST.4, 5, 23-26, 28-36).** New tests for resuming an
  interrupted ring, shell or block fill, bright-star scatter edge cases and
  interrupted scatters, `--force` scatter then fill, every `generate.py`
  argument error by its message with each limit at its maximum, the limits
  staying consistent, grid seams and the nucleus, sector placement running
  out of room, and direct tests of the system builder, moon, Kepler, star
  and phenomenon helpers. The old known-bug sweeps run over hundreds of
  seeds, the Tier 2 reports are now hard asserts or strict xfails, and the
  grid boundary tests also run at the real 4 pc sector edge.

### Fixed
- **The Kepler solver could return an unconverged answer.** When Newton's
  method runs out of iterations it now falls back to bisection, and it
  refuses an eccentricity outside [0, 1).
- **A float rounding slip could put a bright star in a zero-weight bin.**
- **A fill after a failed bright-star scatter could place the same bright
  stars twice.** The backfill now clears unbuilt leftovers in the cells it
  draws for when no scatter threshold was recorded.
- **Two colony and population tests failed on about 2% of draws (TEST.69).**
  One searched the page for a planet name with an apostrophe without
  escaping it; the other let a second capital keep random ages, which
  sometimes founded a polity with no systems.
