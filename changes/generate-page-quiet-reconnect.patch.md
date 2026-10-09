### Fixed
- The Generate page no longer says "Connection lost; reconnecting..." every time the job log stream rolls over. The server ends each stream every 40 seconds on purpose and the browser resumes at the right line; the note now appears only if the stream stays down for 8 seconds, and clears when it is back.
