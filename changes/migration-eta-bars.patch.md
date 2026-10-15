**Headline:** The migration's progress bars now show a time remaining.

### Changed
- The migration's two progress bars show a time remaining (PERF.102). A step that names each phase of its batch (the spherical conversion, `Converting to shells: ...`) used to open a new inner bar for every name, so no bar ever finished a unit and neither had a rate; the inner bar now stays across names (it starts over only when its total changes or its count goes backwards). The step bar moves by the share of the step done, so a run with one long step (one migration pending) has a rate on the top bar too (`stage_progress.StageProgress.detail`).

TODO items resolved: PERF.102 / Added: none
