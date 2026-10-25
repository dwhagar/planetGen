**Headline:** The words "rough guess" are gone from every time estimate; a time that has not been timed yet is shown as a range and says so.

### Changed
- **Time estimates no longer say "rough guess" (PERF.199, bug).** Where this server has not timed the work yet, the Generate page's estimate and review page now show the range with its units and "(not timed yet on this server)", for example "40 s to 5 m 00 s (not timed yet on this server)", instead of one figure and "a rough guess, anywhere from ..." or "a rough default". The Generate job's whole-job bar reads "roughly A to B left (not timed yet)", the Queue page's time left reads "A to B (not timed yet)", the Galaxy banner says "(nothing has been timed yet)", the deep database check's estimate says "not timed yet", and the estimate's source in the run record API is now "not timed yet" instead of "rough guess".

TODO items resolved: PERF.199. Added: none.
