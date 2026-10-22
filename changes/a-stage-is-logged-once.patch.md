**Headline:** A run's log names each stage once: the stage list and the repeated "Stage 4 of 13" line under the job's banner are gone.

### Changed
- **The run output no longer repeats the stage under the banner (PERF.145).** A command used to print `1 stage (stages 4 to 4 of 13): 4. Plan the galaxy.` and then `Stage 4 of 13: Plan the galaxy.` straight after the job's own `=== Stage 4 of 13: Plan the galaxy ===`. The list of stages is gone from every `galaxy` and `plan` run; a step of one stage in a web job prints no line of its own (the banner is the mention), a step of several stages still prints one `Stage n of N: label.` line for each as it enters it, and a command run on its own prints that line for each stage. A stage that is skipped still says why. The Generate page's stage text and the progress file are unchanged.

TODO items resolved / Added: PERF.145 / Added: none
