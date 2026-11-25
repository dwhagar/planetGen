**Headline:** Tests: the every-button test looks at each page in a page of its own, so no load is refused

### Fixed

- TEST.836: `test_every_plain_button_does_something` sent one browser page to the address again for every button, and after about ten loads the browser refused the next load's scripts (`ERR_INSUFFICIENT_RESOURCES`: each load makes about 250 requests, the Shoelace modules among them, and the requests of the loads before it stay outstanding in that page's renderer until the server has answered them). The Shoelace elements were then never defined and the click waited for a button that was not there (`web.generate_task`, `web.admin_queue`, `web.admin_queue_jobs`, `web.admin_queue_run`, `web.polities`, `web.sector`, `web.system` came and went red with the load). Each look is now a page of its own; those pages pass (12 passed) and a new page opened right after a refused load always loaded, which is what pointed at the page and not the server.
