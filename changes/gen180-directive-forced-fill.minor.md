**Headline:** A sector that drawing cannot make to a habitable-worlds minimum is now met by a forced fill, and an impossible request is refused up front.

### Added
- **Forced fill for a habitable-worlds directive (GEN.180).** After K draws (`tuning.DIRECTIVE_ATTEMPTS`, 200) a
  `--directive habitable>=N` still missed builds the sector once more with a habitable world forced into N systems, and the run
  says `met_forced`. Before, the closest draw was kept and the run only warned.
- **Impossible directives are refused before any drawing.** A minimum a sector reaches with a probability under one in a
  billion, judged from the Poisson tail of the system count and the compound-Poisson table of habitable worlds, and that the
  forced fill cannot reach, stops the command line (and the Generate page's estimate) with the reason.
