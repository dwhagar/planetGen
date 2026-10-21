**Headline:** Every admin page now sits in one shell with a rail, a search box that finds any page, task or setting (Ctrl+K), and a thin strip showing the running job.

### Changed
- **One admin shell.** The tab row is replaced by a rail of six destinations (Overview, Generate, Settings, Jobs, Data, Account), a search box and a Compact switch (compact by default, remembered in this browser). Every page under `/admin`, and Account, wears it; the registry (`web/adminregistry.py`) lists each page, task and setting once, and a test fails when an admin page is left out of it.
- **Search and palette.** `/admin/search` finds a page or setting by what it does and works without script. Ctrl+K (or `/`) opens the same list as a palette; Enter opens the result, and nothing runs from it.
- **Keys.** `g` then `g`, `s`, `j` or `d` goes to Generate, Settings, Jobs or Data; `?` lists the keys; Esc closes. They are off while typing in a field.
- **Job strip.** The running job, its stage, a bar, the elapsed time, a link and Cancel show above every admin page and follow the job live.
- **Queue one next job.** Starting a second job while one runs no longer fails: it waits as the one job behind it, shown under the strip with Discard. It starts by itself when the running job succeeds or is cancelled; if that job failed it is held, and the strip offers Start now. A third is refused until the waiting one starts or is discarded. (The Galaxy Map's one-click buttons still say a job is running.)
