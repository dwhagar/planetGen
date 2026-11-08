**Headline:** The Generate page's size and time estimate, which runs in the work queue since ADM.141, is on the queue's allowlist, so the worker no longer refuses it.

### Changed
- `planetgen.web.lib.estimatejob:work` is added to `allowlist.API_FUNCTIONS`. ADM.141 moved the estimate into the work queue but left the function off the list, so `api_jobs.execute` refused it and `test_every_function_the_api_queues_is_on_the_list` failed on main.
