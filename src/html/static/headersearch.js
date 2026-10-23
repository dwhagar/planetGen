// html/static/headersearch.js
//
// The header's search button (UX.146) drops the search box down; this puts the cursor in it, so a visitor can type
// at once. The dropdown is the <sl-dropdown class="site-search-menu"> of templates/base.html.

document.addEventListener("sl-after-show", (event) => {
  const menu = event.target;
  if (!menu || !menu.matches || !menu.matches("sl-dropdown.site-search-menu")) return;
  const input = menu.querySelector("input[type=search]");
  if (input) input.focus();
});
