**Headline:** The data check in front of every Generate job no longer fails because another check holds the check place, and a refusal now names the run that holds it.

### Fixed
- **The data check failed with "run ... is running".** The check only reads files, so when the check place is taken it now logs a warning naming the holder and checks the data without a run record, instead of stopping the whole job.
- **A refused run now says who blocks it.** The message gives the holder's run id, title, kind, start time (UTC) and the process that follows it, instead of the id alone.
- **A job waiting for the machine now says who it waits for.** Under "Waiting at position N for run ..." the job log gives that run's title, kind, start time and process, and says so when the run belongs to the same job.
- **A run started inside a run of the same process no longer waits for itself.** It fails at once, naming the holder and logging where it was started from, instead of printing "Waiting at position 1" for ever (ADM.136).
