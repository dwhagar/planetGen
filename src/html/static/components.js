// html/static/components.js
//
// The shared Shoelace component set (UX.40): buttons, menus and, as pages
// adopt them, dialogs and form fields. base.html loads it as a module on
// every page. Shoelace is vendored under vendor/shoelace/ (see
// vendor/THIRD_PARTY_NOTICES.txt) and registers its custom elements when
// each file is imported. A page that needs another component adds it to
// COMPONENTS in scripts/vendor_shoelace.py, re-runs that script and
// imports it here.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;

// The icon libraries first, so no icon is requested through the default resolver.
await import(`./shoelaceicons.js${VERSION_QUERY}`);
await Promise.all([
  "button/button",
  "icon-button/icon-button",
  "icon/icon",
  "dropdown/dropdown",
  "menu/menu",
  "menu-item/menu-item",
  "divider/divider",
].map((name) => import(`./vendor/shoelace/components/${name}.js${VERSION_QUERY}`)));
