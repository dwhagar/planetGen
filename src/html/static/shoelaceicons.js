// html/static/shoelaceicons.js
//
// Points Shoelace's icons at files under vendor/shoelace/assets/: the
// "default" library (icons a page names, such as the gear) at icons/, and
// the "system" library (the check in a menu item, the caret on a button) at
// system/. Shoelace's own "system" icons are data: URIs, which the site's
// Content-Security-Policy (`default-src 'self'`) refuses to fetch. components.js imports this first, before
// any component, so no icon is requested through the default resolver.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { registerIconLibrary } = await import(`./vendor/shoelace/utilities/icon-library.js${VERSION_QUERY}`);

const folder = (path) => (name) => new URL(`./vendor/shoelace/assets/${path}/${name}.svg`, import.meta.url).href;

registerIconLibrary("default", { resolver: folder("icons") });
registerIconLibrary("system", { resolver: folder("system") });
