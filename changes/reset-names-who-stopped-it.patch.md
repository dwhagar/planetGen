**Headline:** A reset, data check or database check that is ended by SIGTERM now says who sent it, and the reset says that running it again finishes the wipe.

### Changed

- `planetgen reset` and the data and database checks end on a `Stopped` interrupt that carries the reason the canceller left (update guard, stop-all), or says the signal left none, in place of a bare traceback. The reset prints one line and exits 143.
