### Changed

Headline: a worker's memory is counted with its parent, and the collector no longer copies its tables

The largest-worker figure in the run report is now the horse's and its worker's memory together (shared pages counted once), which is what one worker takes of the host (a sector worker measured 2.2 GB, not the 1.5 GB the horse alone showed). A worker freezes what it holds before it forks a horse, so the garbage collector does not write to every page of its warm tables in each horse. No generated output changes.
