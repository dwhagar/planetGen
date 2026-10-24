# Admin control suite redesign: Generate, Settings and the rest (research report and proposal)

*Times in this note are the admin redesign thread's own clock as it reported them.*
Written 2026-10-21, 6:20 PM EDT, extended 6:45 PM EDT, from a read-only pass over dwhagar/planetGen at main 48e47b5 (Release 8.0.1864 plus TODO PRs through #2348). No code was changed and nothing is built until Boss answers the cards in section 7.

Boss's asks: (5:57 PM EDT) the Generate page is horrible to use and has grown too complicated for its simple format; the Settings page is incomprehensible and reads like a Linux config file; list only what can be changed, or hide the rest behind an option that is on by default. (6:01 PM EDT) do the redesign right away in its own thread; research ways to pack a lot of options into one admin control suite that includes the Generate control suite; "the interface needs to get out of the way"; research, report and cards before any building.

Sections 1 and 2 describe what exists. Section 0 is the research and the suite design. Sections 3 and 4 detail the Settings and Generate pages inside that suite. Section 5 lists the items, 6 the overlap with existing TODO items, 7 the choices for Boss.

## 0. Research: packing many options into an admin suite

Method note: this is a design review from established products and published interface guidance, not a fetch of web pages (the environment's network is limited and nothing here depends on a figure that needs a citation). Each pattern is rated against this project's problem: roughly 90 settings, seven Generate forms with about 80 controls, about eight admin pages, one admin at a time, mostly on a desktop but sometimes a phone, and a rule that pages work without script.

### 0.1 Patterns and what they would do here

| Pattern | Seen in | What it gives | Limit | Verdict |
|---|---|---|---|---|
| Command palette (one box reaches every page, task, setting, recent job) | VS Code, GitHub, Linear, Raycast | Zero chrome; scales to hundreds of options; gets out of the way | Not discoverable alone; needs a no-script fallback and proper combobox accessibility | **Adopt as the primary way in**, beside a visible rail |
| Typed commands with arguments that land on a review step | Raycast, Slack slash commands, GitHub Actions "run workflow" | "add sectors radius 12" opens the Review step pre-filled; three keystrokes to run | The parser must stay small and forgiving | **Adopt**: the palette fills the same form the pages use, never runs anything by itself |
| Search-first settings with a "modified only" filter | VS Code, macOS System Settings, Chrome | No scrolling; changed values findable; matches in labels, summaries and keywords | Needs good keywords | **Adopt** for Settings |
| Progressive disclosure (basic, then advanced) | Stripe, GitHub settings, Home Assistant | Fewer controls on screen, same power | An "advanced" bucket gets abused | **Adopt** with a rule: advanced only for options a normal change never needs, and search always reaches them |
| Accordion folds with remembered state | this site's ADM.4 | Hides length | Hides content without reducing the number of choices, which is why Generate still feels bad | **Drop** as the main device |
| Dense table with inline edit (click the value, Enter saves, Esc cancels) | Airtable, Grafana, Kubernetes dashboards | Many rows per screen; value, default and state in one line | Risky for options that need a review (rate limits, restart-only) | **Offer** as a compact density for Settings (answer 7) |
| Master and detail (list left, detail right, no page loads) | cloud consoles, mail clients | Moving between settings groups or jobs costs nothing | Needs a good phone layout (list then detail) | **Adopt** for Settings groups, Jobs and Backups |
| Persistent activity strip | CI dashboards, Vercel, GitHub Actions | The running job is visible from every page, collapses to a thin bar | One more moving part | **Adopt**: one strip across the whole suite, not just Generate |
| Saved presets and "run again" | Jenkins "build with parameters", Blender presets, GitHub saved replies | The common run is one click; long forms stop being re-typed | Needs storage and a way to name and delete | **Adopt** as Run again first, named recipes second (answer 11) |
| Plan, then apply with a diff | Terraform, git | Review is the safety; destructive actions get a typed confirm only when they cannot be undone | Needs a clear summary, not raw JSON | **Keep and reuse**: Review is the same component for settings changes and for Generate runs |
| Undo instead of confirm for reversible changes | Gmail | Fewer dialogs | Not available for restart-only or destructive settings | **Limited**: only for light choices (density, hide toggle) |
| Keyboard shortcuts with a help overlay (`/` to search, `g` then a letter to go) | GitHub, Gmail, Linear | Power use without UI weight | Must not clash with typing; must be listed | **Adopt** a small set (answer 6) |
| URL-addressable state (every setting, task step, job has a link) | GitHub, Grafana | Share and bookmark; Back works | Needs a stable key per item | **Adopt**; settings keep their config key as the anchor |
| Responsive bottom sheets and single-column steps | mobile consoles | Works on a phone | More CSS | **Adopt** for the flow steps and palette |

### 0.2 What "gets out of the way" means here (design rules)

1. **The default screen is status plus the next action.** Nothing else is drawn until asked for. No section of long forms is open on load.
2. **One way to find anything.** The palette, the rail, the hub page and the settings search all read the same registry, so a new page, task or setting appears everywhere with no extra work (this is how `settings.opt` already feeds the settings page, the API and the CLI).
3. **A control that cannot do anything for this admin is not drawn.** The "cannot be changed" items sit behind a switch that is off by default (Boss's rule).
4. **Each choice is made once.** The same shape fields, scatter limits, prevalence and overrides exist once and are reached from whichever task needs them.
5. **Safe things are cheap, dangerous things are far.** Check, add sectors and settings changes sit up front. Reset and New galaxy are last, red and need the typed database name.
6. **The running job follows the admin.** A thin strip on every admin page; the full log opens from it.
7. **Everything works without script;** script adds the palette, live counts and the strip's live updates.
8. **Compact by default,** one density toggle, remembered in this browser.

### 0.3 The suite

```
 top bar:  [ Admin ]  [ Search or run...  Ctrl+K ]                 [job strip: Add sectors 3/6 · 42% · ETA · Cancel]  [account]
 rail:     Overview | Generate | Settings | Jobs | Data | Account
 content:  one column, up to two panes (list left, detail right) at wide widths
```

| Destination | Holds today | Notes |
|---|---|---|
| Overview | admin hub, maintenance mode, API keys | status, shortcuts to the six common tasks, maintenance toggle |
| Generate | the seven forms, one-off system, Recent jobs | section 4: six task cards, five-step flow, History |
| Settings | config.json options | section 3 |
| Jobs | Work Queue (`/admin/queue`, tree, confirmation), job pages | one place for every job including Generate's; ADM.95 builds on it |
| Data | Stats, Backups and restore, Images, Facility types | list and detail panes |
| Account | username, password | unchanged |

The **admin registry** is the new foundation: a small module where each page, task and setting declares a label, one-line summary, keywords, the role that may use it, a danger level (safe, rewrites rows, destroys) and the URL it lives at. Settings keep their model as the source for their entries. The palette searches it on the server (so the no-script path is a plain search page, `/admin/search?q=`), the rail and hub are drawn from it, and a test fails when an admin route is missing from it.

**Palette contents:** pages; tasks (verbs such as "Add sectors", "Check the database", "Make one system"); settings by label, summary and keywords (open the row, and switch or number options can be set from the palette); recent jobs; sectors and systems by name (the existing `/galaxy/locate` search); maintenance mode on and off; the keyboard help. Typed arguments that match a task's fields open that task's Review step pre-filled ("add sectors 12 pc around <sector name>"); nothing runs from the palette without the Review step.

**Keyboard set (small on purpose):** `Ctrl+K` or `/` search; `g` then `g`, `s`, `j`, `d` to go to Generate, Settings, Jobs, Data; `?` help; `Ctrl+Enter` run on a Review step; `Esc` closes. All of it off while typing in a field.

### 0.4 Risks the design must handle

* The palette must be a real accessible combobox (ARIA authoring practices pattern), usable by keyboard and screen reader, with the plain search page as the no-script path.
* A pre-filled command must show exactly what it will do on the Review step; the parser may refuse an ambiguous command rather than guess.
* Settings that need a restart or are overridden by an environment variable must say so in the palette result too, not only on the page.
* The job strip must not poll a heavy endpoint on every page: reuse the existing Server-Sent Events stream and fall back to the light status call only while a job runs.
* Moving other admin pages onto the shell is the large part of the work (ADM.93); the order puts Generate and Settings first and splits the rest into one item per destination.

## Where the redo stands in the plan

Short answer to "do we have a redo somewhere": yes, but late. ADM.55 (audit), ADM.92 (admin shell), ADM.93 (move every page onto it) and ADM.56 (umbrella) are all in 8.6 "edit everything", behind UX.115. Nothing in the plan redesigns Generate or Settings themselves, and ADM.28 ("A simpler Generate page", done, PR #991) was the last try, before the page grew to seven forms. This draft is the redesign itself, with a proposal to pull it forward.

## 1. What the two pages do today

### Generate (`/admin/generate`)

Files: `web/generate_page.py` (1,536 lines), `templates/generate.html` (571), `generate_job.html`, `generate_system.html`, `partials/job_status.html`, `static/generatejobs.js`, `generatefolds.js`, `generatecustomize.js`, `generatestarmix.js`, `generatebuttons.js`, `limitsliders.js`, `web/limit_sliders.py`, `web/system_page.py`.

One URL, seven POST forms told apart by a hidden `action`. Every form becomes a background job (one at a time) that runs the same CLI an admin would type. Sections, top to bottom, all folded with `<details>` and remembered per browser (ADM.4):

| Section | What it does | Controls |
|---|---|---|
| Badges | Database name, planned or not, sector count, bright-star state | none |
| Estimate (interstitial) | After Generate or deep check, shows size, time, disk, then asks again to confirm (PERF.3, ADM.33) | Generate, Generate anyway, Cancel |
| Current job | Step, stages, two progress bars, log terminal, cancel, copy, download | only while a job runs |
| Galaxy shells | Read-only table of shells | none |
| One-off system | Link to a second page that builds one system and never saves it | link |
| New galaxy | Wipes the database, plans, generates a first neighbourhood, scatters | 9 shape fields, 10 prevalence and 3 star-mix fields, 5 overrides, 2 scatter sliders, neutron star and black hole limit, 3 segments (core, halo, dense cells) with 5 more sliders, skip checkbox, radius, minimum density, typed confirmation. Customize is a modal with tabs |
| Generate sectors | Adds sectors, five modes (random start, shell, around a sector, one address, span), the centre named three ways | about 25 fields, the same prevalence and override groups again |
| Plan the galaxy | Rebuilds the density skeleton, then scatters | the same 9 shape fields and the same scatter controls again |
| Redo scatters | Redoes mass, luminosity or phenomena scatters on the current plan | 3 checkboxes, each with its own slider group |
| Add a dimmer step | Extends the bright-star floor downward | 1 field |
| Check the database | Plain or deep check | 1 checkbox |
| Reset | Wipes everything | typed confirmation |
| Recent jobs | Table of the last jobs | pager |

Why it is hard to use, from the code:

1. It is organised by the CLI's commands (plan, galaxy, redo, band, reset), not by what an admin wants to get done.
2. The same controls appear on two or three forms: the shape fields on New galaxy and Plan, the scatter limits on New galaxy, Plan and Redo, prevalence and overrides on New galaxy and Generate sectors. The admin cannot tell which copy matters.
3. Eight things that change the stored plan, the sector list and the whole database sit in one scroll with equal weight; Reset and New galaxy look like the safe forms above them.
4. The running job is far from the forms that started it, and the forms lock without saying which one is blocked or offering a way to prepare the next run.
5. The estimate and the confirmation arrive as a separate interstitial page that replaces the form, so a changed mind means losing the form.
6. Radio fieldsets, `Customize` modals, folds and sliders each solve a different part of the same overload.
7. Most fields have no unit scale or live consequence beside them (the ceiling cost hint is the one good example).

### Settings (`/admin/settings`, ADM.84, ADM.96)

Files: `web/settings_page.py`, `templates/admin_settings.html`, `admin_settings_review.html`, `util/settings.py` (the model), `util/settings_view.py`, `util/settings_writer.py`, `api/settings.py`.

One form over the whole `config.json` model: 89 options in 16 categories, each drawn as a "setting" block: the title made from the path (`tile cache > warm zoom levels`), the unit in brackets, the control, the description, then a meta line with the raw key (`tile_cache.warm_zoom_levels`), `Default: 14`, limits as text, "needs a restart", "overridden by PLANETGEN_X", "set in config.json", a secret tag and a reset checkbox. Read-only options are printed as a `<code>` value plus the reason. Save is two steps (review diff, confirm) with a backup and audit rows, which is good and should stay.

Why it reads like a config file:

1. Titles come from key paths, not words. The key, the default and the source are printed on every row.
2. Control choice is by type only: an integer is a text box, a duration is "seconds" in a text box, a rate limit is the string `600 per minute`, a list is a textarea, a colour is a text box.
3. 31 of the 89 options cannot be changed from the web, yet they are printed in full among the ones that can.
4. Nothing is grouped by purpose: caches are 17 options in one block, search engines 14, rate limits 7 with raw limiter strings.
5. No search, no sticky save, no way to see what you changed before you submit.

## 2. Every setting classified

Counts come from running the model (`settings.iter_fields` with `settings_writer.allowed`) on the repo.

| Class | Count | Meaning |
|---|---|---|
| A. Changeable by an admin from the web | 58 | `x-editable` true, role admin |
| B. Owner only | 7 | `x-editable` true, `x-min_role` owner. Nobody in the web can change them until the Owner account exists (USR.31 is an open 8.1 item), so today they are read-only for everyone |
| C. Server only | 24 | `x-editable` false: the file can only be edited, or `python3 -m planetgen.cli.config set` run, on the server |
| Total | 89 | |

Of the 58 in class A, 13 need a web server restart to apply (the rate limits, the page cache's switches and sizes), and 6 have an environment variable that overrides the file (`ratelimit.default`, `ratelimit.storage_uri`, `cache.max_mb`, `cache.warm_on_change`, `page_cache.enabled`, `page_cache.stored_counts`), in which case a save does nothing to the value in use.

### Class A, by proposed group (what the admin sees)

| Proposed group | Options | Controls proposed |
|---|---|---|
| Site identity | `site_name`, `site.tagline`, `site.language`, `web.title_template` | text, text, language picker, text with a live preview of a page title |
| Look | `web.theme_color_light`, `web.theme_color_dark`, `icons.manifest_enabled`, `icons.favicon_svg`, `favicon_png`, `apple_touch_icon` | colour pickers, a switch; the three icon paths replaced by a link to the Images page (they are already managed there) unless the path is the only way |
| Search engines | `seo.indexing`, `detail_pages`, `description`, `keywords`, `json_ld`, `block_ai_crawlers`, `sitemap_enabled`, `sitemap_max_urls`, `sitemap_cache_hours` | segmented choices, switches, chips for keywords, number with unit |
| Search engines, advanced (collapsed) | `seo.robots_extra`, `robots_override`, `security_txt`, `humans_txt`, `privacy_note` | text areas, with a warning on `robots_override` |
| Link previews | `og.enabled`, `site_name`, `image_alt`, `default_image`, `twitter_card`, `twitter_site`, `generated_cards`, `card_cache_mb` | switch, text, choice, slider with MB |
| Speed and disk | `cache.max_mb`, `cache.warm_on_change`, `cache.warm_settle_seconds`, `tile_cache.warm_levels`, `tile_cache.warm_zoom_levels`, the three `prefetch_*`, `page_cache.enabled`, `page_cache.warm_max_objects`, `page_cache.stored_counts` | slider with the disk size in GB, switches, steppers with the tile size shown ("level 14: 32 pc tiles") |
| Speed and disk, advanced (collapsed) | `page_cache.max_entries`, `max_mb`, `stamp_seconds`, `max_age_seconds`, `shared_seconds` | numbers with unit, "needs restart" mark |
| Abuse protection | `ratelimit.default`, `storage_uri`, and the five `ratelimit.pages.*` | a rate control: a number, a unit (minute, hour, day) and an "unlimited" switch per row; storage URI behind Advanced |
| Housekeeping | `jobs.keep`, `settings.backups_kept`, `mysql.statement_timeout_seconds` | steppers with plain sentences ("keep the last 50 finished jobs") |

### Class B (Owner only, 7): hidden by default

`base_url`, `login_allowlist`. They appear in the "cannot be changed" list with "Only the Owner can change this" until USR.31 and ADM.86 land, then in an Owner section.

### Class C (server only, 24): hidden by default

| Group | Options |
|---|---|
| Database connection | `mysql.host`, `port`, `user`, `password`, `database`, `database_prefix`, `control_database` |
| Redis | `redis.url` |
| Logging | `debug`, `log_file`, `log_dir`, `log_rotation` |
| Proxy | `proxy_fix.x_for`, `x_proto`, `x_host` |
| Security | `admin_cookie_insecure`, `secret_key` |
| Paths | `api_base_url`, `web.assets_dir`, `og.card_cache_dir`, `tile_cache.dir`, `jobs.dir`, `jobs.python` |
| File | `config_version` |

These move into a single "Server files" disclosure: a section that explains in one paragraph how they are changed (`python3 -m planetgen.cli.config set key value`, then restart), then the list with values (secrets as "set" or "empty").

## 3. The proposed Settings page

Principles: words not keys, the control that fits the value, one line that says what it does and what happens at both ends, nothing listed that cannot be changed unless asked for.

Layout:

1. Header line: "Settings", a search box, and a switch "Show settings that cannot be changed (31)", **off by default**. When on, the 31 appear in their groups, greyed, with the reason in plain words (server only, Owner only, set by an environment variable) and no control. The choice is remembered in the browser and also in the URL (`?show=fixed`) so a link can carry it.
2. Left rail (top bar on a phone): the groups from section 2 with a count of changed-but-unsaved options beside each. A group with nothing visible is not listed.
3. Each setting is a row: **plain label**, a one-line description (the first sentence of the model description), the control with its unit written in the control ("1,024 MB"), and a small "Details" disclosure holding the long help, the default, the allowed range, the config key, and where the value in use comes from. "Reset to default" is a link in the row, shown only when the value differs from the default.
4. Controls by kind: switch for yes or no; segmented choice for 2 to 4 choices, select for more; slider with a number box for sizes and counts (this is ADM.97, using the existing `limit_slider`); duration controls that say "15 s" or "6 h"; the rate control above; chips for lists; colour picker for colours; text area only for the long text options. Secrets are never in class A.
5. A sticky save bar at the bottom: "3 changes. Review and save." It stays disabled until something changes. The two-step review stays: a plain-language diff ("Cache size: 1,024 MB to 2,048 MB"), the list of changes that need a restart, a note for any that an environment variable overrides, then confirm. The backup and audit behaviour stays untouched.
6. Environment-overridden options (6) are shown with a warning badge "Set by PLANETGEN_X on the server, so a change here has no effect" (answer 9).
7. If `config.json` cannot be written (`view.problem`), the page opens with that explanation at the top and every control read-only, instead of hiding the problem.

Wording example, today versus proposed:

| Today | Proposed |
|---|---|
| `tile cache > warm zoom levels (levels)`, then "Finest tile level of the zoom run `planetgen warm` builds from the opening view in to the galaxy's centre, on top of the levels above (level 14 is tiles 32 pc across and the finest the map's cut uses; at most 15; 0 builds no run)." then `tile_cache.warm_zoom_levels  Default: 14` | **How close the map is prepared in advance** "Zoom in to this level when pre-building the map" [stepper 14] "Level 14 shows tiles 32 pc wide, the closest the map goes. 0 skips the zoom run." Details: default 14, up to 15, key `tile_cache.warm_zoom_levels` |
| `ratelimit.pages.galaxy_tiles`, text box "600 per minute" | **Map loading limit** [600] per [minute] per visitor, "Stops one visitor from flooding the server; applies after a restart." |

## 4. The proposed Generate page

Principles: organise by what the admin wants to do, one running job always visible, one place for each control, safe actions apart from dangerous ones, and the estimate is a step of the flow instead of a page that replaces it.

### Frame (every Generate screen)

1. **Status strip**: database name, "Planned, edge at shell N" or "Not planned yet", sector count, scatter state. It replaces the badges and the Galaxy shells fold (the shell table moves behind a "Galaxy at a glance" link).
2. **Current job panel**, pinned under the strip on every Generate page while a job runs: title, "step 3 of 6", one overall progress bar with the estimate to finish, Cancel, and the log as a collapsed terminal (Copy and Download stay). When idle, it shows the last job's result in one line with "Show on Galaxy Map", "Full output" and "Run again with the same settings".
3. Forms stay usable while a job runs. Only the final Run button is disabled, with the reason beside it ("Generate sectors is running; Cancel it or wait"). This replaces "the forms unlock when the job finishes" (answer 10: queue one next job, item K9).

### Start screen: six task cards

| Card | Replaces | Danger |
|---|---|---|
| Add sectors to the galaxy | Generate sectors | none |
| Start a new galaxy | New galaxy (reset, plan, first neighbourhood) | destructive |
| Change the stars and phenomena scattered across the galaxy | Plan, Redo scatters, Add a dimmer step | rewrites stored rows; sectors already generated are left alone |
| Make one system (not saved) | One-off system | none |
| Check the database | Check the database | none |
| Wipe the galaxy | Reset | destructive, shown last, red, apart from the others |

Recent jobs becomes a **History** tab next to Start with a link to the Work Queue (ADM.95 and the queue pages already cover the same jobs).

### Every task is the same five steps

**What** (only where there is a choice) → **Where** → **How big** → **Review** → **Run**. A progress rail on the left shows the steps; the draft is carried in the URL query and hidden fields, so Back works, a step can be bookmarked, and the page works without script.

* **Add sectors**: What: "Around a place" or "By shells, bands and wedges". Where for "Around a place": a random populated spot, a named sector (the existing finder), a sector address, or a position in pc, one picker with four tabs (this is the five modes and the three centre choices collapsed into two plus a picker). Where for the grid choice: shell, band and wedge or a range of each, which also covers "one address" and "a whole shell". How big: radius with a slider and a live line "about N sectors, M systems"; or the limit. Customize (collapsed, "Advanced"): prevalence, star mix (with the 100% total checker) and overrides, each with "Usual: X%" kept.
* **Start a new galaxy**: Where: the shape. Every shape field shows "The seed decides" with a "Choose it myself" switch, so the pinned values are obvious and nothing is preset (nothing fixed from any real galaxy). How big: first neighbourhood radius and minimum start density. Advanced: scatter limits (below), prevalence, overrides. Review ends with the typed database name.
* **Change what is scattered**: What: "Rebuild the plan and scatter" (today's Plan), "Redo only some scatters" (today's Redo, with the three kinds as three rows with switches), "Add a dimmer step" (one field). The scatter limits (the two galaxy-wide sliders, the neutron star and black hole limit, core, halo, dense cells) live here as an Advanced block with one "Reset all to defaults" and the "Stored with this plan" values shown in place, instead of repeated on three forms.
* **Make one system**: the same page, on the same frame; the ten tri-state options become segmented "Decide for me / Yes / No" rows.
* **Check**: plain or deep as a two-option choice with the time the deep check will take.
* **Wipe**: What it deletes, what it keeps (admin logins and API keys), the typed database name, Run.

**Review** is the old estimate page, brought into the flow: a plain-language summary ("Add about 340 sectors around the sector named X, within 12 pc, with the usual prevalence"), sectors, systems, stars, database growth, time with N workers, disk check (the refusal and "Generate anyway" stay), the exact command line in a copyable box, then Run. Changing anything is one Back click, with the form intact.

### What is removed

`generatefolds.js` and the fold macro, the Customize modal (`generatecustomize.js`), the five-radio mode fieldset, the repeated scatter and shape controls, the per-segment reset buttons, the separate estimate interstitial, and the "forms below unlock" lock. `generatejobs.js`, `limitsliders.js`, `generatestarmix.js`, the argv builders in `generate_page.py` and the estimate cache stay, since they are the logic and the live job panel; only the templates and their wiring change. No backward compatibility shim for the old URLs or form field names: the old `action` posts go, and `old_urls` is not extended.

## 5. Items needed

Numbers were assigned by the TODO thread (6:25 PM EDT): S1 to S6 are ADM.101 to ADM.106, G1 to G11 are ADM.107 to ADM.117, K1 is ADM.118, K3 to K6 are ADM.119 to ADM.122, K9 is ADM.123; K2 is ADM.92, K7 is ADM.93 and K8 is ADM.55. Boss moved everything to 8.1 phase 1 (answer 1). Phase: "8.1 P1" means 8.1 Phase 1 (the Admin and API housekeeping and Generate page threads), which is where the settings page (ADM.84) already lives. The argv builders and their tests (`test_web_generate.py`, `test_generatebuttons.py`, `test_web_settings.py`) are the safety net and must keep passing unchanged where they test commands rather than markup.

### Suite (TODO filed ADM.98 research and cards, ADM.99 Generate suite, ADM.100 Settings; ADM.55, ADM.92, ADM.93, ADM.56 moved here from 8.6 as 8.2 Phase 1 with this thread as owner, 6:10 PM EDT)

The S and G items below are the children of ADM.100 and ADM.99. The K items are new and have no number yet; they are the suite pieces the research adds.

| ID | Item | Done when | Needs | Phase |
|---|---|---|---|---|
| ADM.118 | The admin registry: pages, tasks and settings declare label, summary, keywords, role, danger level and URL; the rail, hub and search read it | a test fails for any admin route missing from the registry; the settings entries come from the settings model with no second list | ADM.101 | 8.1 P1 |
| ADM.92 | The shell (this is ADM.92): top bar, rail with six destinations, content width, density toggle remembered per browser, phone layout | the worked example is the Overview page; both themes; no console error; keyboard-only run | ADM.118 | 8.1 P1 |
| ADM.119 | Job strip on every admin page (replaces ADM.107's current-job panel as the shared part) | strip shows a running job from every page, Cancel works, the log opens from it, one SSE stream per page and none when idle | ADM.92 | 8.1 P1 |
| ADM.120 | Command palette and `/admin/search`: server search over the registry, accessible combobox, Ctrl+K and `/`, results for pages, tasks, settings, jobs, sectors | screen-reader pattern test; works with script off through the search page; a restart-only or overridden setting says so in its result | ADM.118, ADM.92 | 8.1 P1 |
| ADM.121 | Typed commands that open a task's Review step pre-filled | a table of command strings to argv, byte for byte against the old forms; ambiguous input is refused with a message | ADM.120, ADM.108 | 8.1 P1 |
| ADM.122 | Keyboard set and help overlay | shortcuts off inside fields; the overlay lists every shortcut; tests for each | ADM.92 | 8.1 P1 |
| ADM.93 | Move the other pages onto the shell (this is ADM.93), one change per destination: Overview, Jobs, Data, Account | each page opens in the shell, same table, form, confirmation and job-log parts, TEST.169 browser tests pass | ADM.92 | 8.1 P1 (Boss answers 1 and 2) |
| ADM.55 | Audit of the pages other than Generate and Settings (this is ADM.55 narrowed): ranked defects | report in `/mnt/project-files/reports/`, items filed for each defect above the line | none | 8.1 P1 |

### Settings

| ID | Item | Done when | Needs | Phase |
|---|---|---|---|---|
| ADM.101 | Words for every option: label, one-line summary, group, control, advanced flag (`x-label`, `x-group`, `x-control`, `x-advanced`) in `settings.opt`; the long help becomes Details; the path-made titles go; the `x-category` list is replaced by groups | a test fails for any option without a label, group and summary; no label equals its path text; `config.md`, `config.json.example` and `config.schema.json` regenerate; wording reviewed in a short Boss-readable list | none | 8.1 P1 |
| ADM.102 | View classes: `settings_view.describe` returns `class` (editable, owner, server, overridden) and `fixed_reason`, plus counts | unit tests for each class on a fixture model; the API returns them; the count of fixed options matches `iter_fields` | ADM.101 | 8.1 P1 |
| ADM.103 | The new Settings page: grouped rows, left rail, search, "Show settings that cannot be changed" switch off by default (URL and browser memory), sticky save bar, per-row reset link, Details disclosure | browser test: loads with 31 fixed options hidden, switch shows them with reasons, search filters rows and groups, a change enables the save bar; works with script off; both themes; keyboard-only run | ADM.101, ADM.102, UX.124 | 8.1 P1 |
| ADM.104 | Controls: switch, slider with number box (absorbs ADM.97), duration, rate-limit control, chips, colour picker | each control round-trips its value through `collect` and the writer; a rate string stays valid for the limiter; invalid input shows its error at the row | ADM.103, existing `limit_slider` | 8.1 P1 |
| ADM.105 | Review page in plain language: "Cache size: 1,024 MB to 2,048 MB", restart list, overridden list | review tests updated; the audit rows and backup are byte-identical to today's | ADM.103 | 8.1 P1 |
| ADM.106 | Server files disclosure and the overridden-option badge | text names the exact command; secrets show only set or empty | ADM.102 | 8.1 P1 |

### Generate

| ID | Item | Done when | Needs | Phase |
|---|---|---|---|---|
| ADM.107 | The frame: status strip, pinned current job panel (idle shows the last result), History tab, Work Queue link | browser test opens a running job from every Generate screen; Cancel works; idle panel shows Run again | none | 8.1 P1 |
| ADM.108 | The task flow engine: Start screen with six cards, step rail, draft carried in the URL, Review step with the estimate, summary and command line, Run | any task reaches Review and Back keeps every value; works without script; estimate refusal and Generate anyway still work (ADM.33) | ADM.107, UX.125 | 8.1 P1 |
| ADM.109 | Add sectors task: two choices, the four-tab place picker, how-big with live count | argv for every old mode and centre choice is reproduced byte for byte by a test table; the finder still lists filled sectors | ADM.108 | 8.1 P1 |
| ADM.110 | Start a new galaxy task: seed-decides shape fields, first neighbourhood, typed confirmation | pinned and unpinned fields produce the same argv as today's blank and filled boxes; confirmation refuses a wrong name | ADM.108 | 8.1 P1 |
| ADM.111 | Scatter task: rebuild, redo some, add a dimmer step; one Advanced block of scatter limits with one reset and the stored values shown | one definition of each limit renders once; Reset all restores the defaults; the three old actions build the same jobs as before | ADM.108, UX.116 (done) | 8.1 P1 for the layout; the controls GEN.301 adds are placed in this block by their own items |
| ADM.112 | Prevalence, star mix and overrides as one shared Advanced block (rendered once) | the 100% total checker works; "Usual" values shown; one definition used by Add sectors and Start a new galaxy | ADM.108 | 8.1 P1 |
| ADM.113 | One-off system on the frame, tri-state rows as segmented controls | `test_web_generate` tri-state field test still matches the CLI list | ADM.107, UX.125 | 8.1 P1 |
| ADM.114 | Check and Wipe tasks; History tab | Wipe is last and red and needs the typed name; History lists the same rows as Recent jobs | ADM.107 | 8.1 P1 |
| ADM.115 | Run again from History: re-run a finished job with its stored settings | needs the job record to hold the form; a failed job's Run again goes back to its Review step | ADM.108, ADM.114 | 8.1 P1 or later (answer 11: Run again plus named recipes) |
| ADM.116 | Remove the old Generate page pieces and update `docs/reference/html-interface.md` and the architecture table | the removed files are gone; no unused CSS; docs match | ADM.109 to ADM.114 | 8.1 P1 |
| ADM.117 | Browser tests for every task end to end, in both themes, without script, and on a phone width | a failing step fails the test with a screenshot | ADM.109 to ADM.114 | 8.1 P1; joins TEST.168 and TEST.169 |

## 6. Overlap with existing TODO items

| Item | State | How it relates |
|---|---|---|
| ADM.55 Admin usability audit | moved by TODO to 8.2 Phase 1, owner this thread (6:10 PM EDT) | The Generate and Settings audit is this report; ADM.55 narrows to the other admin pages (ADM.55) |
| ADM.92 Admin shell, ADM.93 Move every page onto it, ADM.56 umbrella | moved by TODO to 8.2 Phase 1, owner this thread | The shell is ADM.92 and the page moves are K7. ADM.107 is not a second shell: its job panel becomes the shared job strip ADM.119 |
| Goal issue | existing #1688 "Admin tools: editing, configuration and incomplete uploads" lists ADM.13, 32, 35, 43, 44, 50, 55, 56, 92, 93 and TEST.169 | TODO is filing a new goal issue for ADM.98 to ADM.100 and the suite; #1688 stays for the editing items |
| ADM.97 Settings slider | open, 8.2 | Absorbed by S4. Move it to the same phase as ADM.104 |
| ADM.86 Owner-tier settings, USR.31 Owner stub | open, 8.1 | Class B (7 options) stays hidden until these land; ADM.102 reads the role |
| ADM.28 Simpler Generate page | done | Its Customize modal and layer specs are replaced by the Advanced blocks |
| ADM.4 folds, ADM.16 prevalence, ADM.45 star mix | done | Their logic and argv builders are kept, their layout is replaced |
| UX.116 limit sliders | done | The sliders are reused in ADM.104 and ADM.111 |
| UX.124 field look, UX.125 Generate on the field macro, UX.126 admin pages on the field macro, UX.49 Shoelace fields | open, 8.1 | ADM.108, ADM.113 and ADM.103 should be built on these; UX.125 becomes part of ADM.108 and UX.126 is split so settings are covered by ADM.103 |
| UX.134 Generate shows the drawn galaxy type, UX.137 halo and satellite switches, GEN.301 scatter controls | open, 9.4, 8.5, 8.1 | These add Generate controls: they land in the new "Start a new galaxy" and Advanced blocks, so each should list ADM.110 or ADM.111 as a prerequisite |
| ADM.95 A remote run is a job; ADM.13 incomplete uploads; ADM.35 control from every screen | open | Use the same frame and job panel from ADM.107 |
| UX.123 Generate page log on the log layer | open | Applies to the new job panel's terminal |
| TEST.168 control sweep, TEST.169 admin browser tests | open | ADM.117 is its Generate share |
| DOC.72 README external files; GEN.590 data check | open | The Run step shows the data check result; it keeps the "refuse unless told to disregard missing files" rule from the 5:53 PM EDT answer |

## 7. Boss's answers (all 13 cards, tapped 6:05 PM to 6:09 PM EDT)

| # | Choice | Answer |
|---|---|---|
| 1 | When | 8.1 Phase 1, Settings first, then shell and palette, then Generate (TODO had filed 8.2 Phase 1; it moves to 8.1 Phase 1) |
| 2 | Shell scope | Every admin page, once (ADM.92 and ADM.93 in this effort) |
| 3 | Palette | Ctrl+K palette plus visible rail, with the plain search page for no-script use |
| 4 | Typed commands | Pre-fill the Review step, never run by themselves |
| 5 | Job strip | Thin strip on every admin page |
| 6 | Keyboard | The small set |
| 7 | Settings saving | Grouped rows with the two-step review, compact table density available |
| 8 | Unchangeable settings | Hidden by default, remembered in this browser and the URL |
| 9 | Environment-overridden | Shown with a warning badge |
| 10 | While a job runs | **Queue one next job** (not the recommended option): forms stay usable and Run queues one job behind the running one |
| 11 | Run again and recipes | Run again from History plus named recipes saved per admin in the control database |
| 12 | Add sectors modes | Two choices with a four-tab place picker |
| 13 | Density | Compact by default with a toggle, remembered in this browser |

Decided without a card: Owner-only options stay hidden with the unchangeable ones until the Owner account exists; Shoelace where it fits (UX.49) and native HTML elsewhere; scatter limits, prevalence, overrides, robots and page-cache internals are Advanced; Reset and New galaxy always need the typed database name.

Consequences of answer 10 (new work): the job runner today allows exactly one active job (`web/jobs.py`, the `active` lock file). Queueing one next job needs a single "next" slot beside the active lock, started by the runner when the active job ends, shown in the job strip with a Cancel, and refused when the slot is full. It must stay one job running at a time. New item ADM.123 below; ADM.108 and ADM.115 and the job strip ADM.119 depend on it. Questions to settle in the build, not cards: the queued job keeps its own estimate check (re-run at start, refuse if the disk no longer fits); a queued Wipe or New galaxy is allowed only with its typed confirmation already given; a failed job cancels the queued one unless the admin chose otherwise (default: cancel it and say so).

**ADM.123 Queue one next job:** a single next slot with start on finish, Cancel, shown in the job strip and History. Done when: a second submit while a job runs queues, a third is refused; finishing starts it; a failure or Cancel clears it with a visible reason; a destructive job queued without its confirmation is refused. Needs K3. Phase 8.1 P1.

## 8. Notes for whoever builds it

* Argv builders and their tests are the contract with the CLI: keep `plan_argv`, `galaxy_argv`, `scatter_argv`, `redo_scatters_argv` and friends, and add tests that pin their output across the change before any template work.
* `settings_view.describe` is shared by the page, the API and the command line, so the new fields must come from the model, not the template.
* The two pages share one rule: a control that cannot do anything for this admin is not drawn.
* No fixed content: shape defaults stay "drawn from the seed" and nothing is prefilled from the Milky Way or the Solar System. Every time shown to the admin is the server's clock formatted as such; this draft's times are EDT.

## 9. What is built

- **Settings words and groups (ADM.101)** are in `util/settings.py` (`opt()` takes `label`, `summary`, `group`, and optionally `control`, `advanced`, `choice_labels` and `slider`); `util/settings_view.py` turns them into the groups and controls the page draws.
- **Settings page (ADM.102 to ADM.106)**: `web/settings_page.py`, `templates/admin_settings.html`, `static/settingspage.js`.
- **Registry (ADM.118)**: `web/adminregistry.py`. Pages and exempt routes are declared there (a test fails for an `/admin` page route in neither); settings entries are built from the settings model; tasks are added when the Generate cards exist.
- **Shell (ADM.92)**: `templates/base.html` draws the rail, the search box and the Compact switch above every page whose `section` is in `adminregistry.SECTION_DESTINATION`. The rail is a bar across the top of the page (it wraps on a phone), not a side column. The Account page wears it too (ADM.93), except while the seeded credentials still have to be changed.
- **Job strip (ADM.119)**: `templates/partials/admin_jobstrip.html`, `static/adminjobstrip.js`, fed by the `strip` field of `/admin/generate/status`.
- **Search and palette (ADM.120), keys (ADM.122)**: `web/admin_search_page.py` (`/admin/search`, `/admin/palette.json`) and `static/adminpalette.js`.
- **Queue one next job (ADM.123)**: `web/jobs.py` (`queue_next`, `queued_job`, `promote_queued`, `discard_queued`; the spec is the jobs directory's `next.json`, mode 0600 like `job.json`), started by `job_runner.run` when a job ends and by any strip poll. It starts after a succeeded or cancelled job and is held after a failed or interrupted one (Start now or Discard); a failed start also holds it.
