**Headline:** A server with no timing history now shows a rough time left, labelled as a guess, instead of "estimating".

### Fixed
- **Time estimates broke on a server that had never timed a stage (PERF.184).** A run's time left, the banner's
  hours and the Generate page's whole-job bar said "estimating" until a first unit finished, and left a later stage out
  of the total. They now start from a rough guess: the stage's units over the workers the host can run at once (its
  cores and memory) at a default cost per unit. The guess is labelled "rough guess" with a wide range (0.4 to 3 times the
  figure), the sector estimate says "anywhere from A to B" while its speed is the default, and each figure is replaced
  by the run's own pace as units finish.
- **A stage with no stored history now learns.** When a stage of a queue run ends, its time and its cost per worker second
  are folded into a running average (`queue/stagetimes.py`), so the next run of the warm-up and the fills starts from
  what the last one took, and the estimate tightens run after run.
