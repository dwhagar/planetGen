**Headline:** Pressing "Stop waiting" on a slow page's loading state now stays stopped, even when a progress answer was already on its way.

### Changed
- `static/slowload.js`: a status request that was in flight when "Stop waiting" was pressed used to come back and draw its progress over "Stopped waiting.", so the page read as still waiting while the Check again button showed. The answer is now dropped when the page has been stopped (checked after the request and after its body is read), so the stop holds and "Check again" starts polling again as before.

### Added
- `tests/js/slowload.test.mjs`: "an answer that arrives after Stop waiting does not undo the stop (TEST.235)" holds the status answer, stops, lets the answer arrive, and checks the stop message, the buttons, that no more polls follow, and that Check again still works. It fails without the change. The browser test `test_stop_waiting_stops_the_polling_and_check_again_picks_it_up` failed on this in two gates (`[light]` and `[dark]`, box load 10) and passed alone.

TODO items resolved: TEST.235 / Added: none
