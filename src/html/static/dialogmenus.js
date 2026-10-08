// html/static/dialogmenus.js
//
// What picking a menu item does (UX.26, UX.31, UX.27):
//   <sl-menu-item data-dialog="ID">  opens the <sl-dialog id="ID">;
//   <sl-menu-item data-href="URL">   goes to that page (a menu of links);
//   <button data-dialog-open="ID">   opens the <sl-dialog id="ID"> too;
// and any [data-dialog-close] button inside a dialog closes it. The Admin
// menu of templates/partials/admin_menu.html and the system page's
// Navigate menu use them.

document.addEventListener("sl-select", (event) => {
  const item = event.detail && event.detail.item;
  if (item && item.dataset.href) {
    window.location.assign(item.dataset.href);
    return;
  }
  const dialog = item && item.dataset.dialog && document.getElementById(item.dataset.dialog);
  if (dialog && typeof dialog.show === "function") {
    dialog.show();
  }
});

document.addEventListener("click", (event) => {
  // A button with data-dialog-open="ID" opens that dialog (a map's Map help).
  const opener = event.target.closest && event.target.closest("[data-dialog-open]");
  const target = opener && document.getElementById(opener.dataset.dialogOpen);
  if (target && typeof target.show === "function") {
    target.show();
    return;
  }
  const button = event.target.closest && event.target.closest("[data-dialog-close]");
  const dialog = button && button.closest("sl-dialog");
  if (dialog) {
    dialog.hide();
  }
});
