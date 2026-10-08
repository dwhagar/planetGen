// html/static/pagecrumbs.js
//
// The breadcrumb of a server-rendered page (base.html's nav.breadcrumbs):
// the page already holds the full list; this folds it onto one line with
// the shared component (breadcrumb.js, NAV.14). Without the script the list
// simply wraps.

const VERSION_QUERY = new URL(import.meta.url).search;
const { createBreadcrumb, stepsOfList } = await import(`./breadcrumb.js${VERSION_QUERY}`);

const nav = document.querySelector("nav.breadcrumbs");
if (nav) createBreadcrumb(nav).set(stepsOfList(nav));
