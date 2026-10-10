### Changed
- A job's log and page count its stages across the whole job: a New galaxy job reads "Step 1 of 12" to "Steps 4 to 12 of 12", not "Step 1 of 4", and each step's own stage lines carry on from the earlier steps' ("Stage 4 of 12"). Every staged process uses the one total.
- The whole-job bar's time left no longer equals the running stage's own until the last stage: each stage still to run counts what earlier runs of it took, else the average of the stages finished so far (in the first stage, that stage's own expected total). The command-line bar does the same without stored times.
