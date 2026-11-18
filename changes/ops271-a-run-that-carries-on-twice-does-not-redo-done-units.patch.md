**Headline:** A run that carries on from a run that itself carried on no longer does its finished units again.
### Fixed
- The units of a stage a run carries on from get the mark `again` in their arguments when the command cannot resume. The next run that carried on from that one compared its cut with those arguments, found them different, and counted none of the units as done, so every unit of the stage ran again. The comparison now leaves the mark out (OPS.271).
