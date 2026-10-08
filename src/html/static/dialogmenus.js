// html/static/dialogmenus.js
//
// Menu items that open a dialog (UX.26, UX.31): an <sl-menu-item
// data-dialog="ID"> opens the <sl-dialog id="ID"> when picked, and any
// [data-dialog-close] button inside a dialog closes it. The Edit menus of
// templates/partials/edit_controls.html are the first users.

document.addEventListener("sl-select", (event) => {
  const item = event.detail && event.detail.item;
  const dialog = item && item.dataset.dialog && document.getElementById(item.dataset.dialog);
  if (dialog && typeof dialog.show === "function") {
    dialog.show();
  }
});

document.addEventListener("click", (event) => {
  const button = event.target.closest && event.target.closest("[data-dialog-close]");
  const dialog = button && button.closest("sl-dialog");
  if (dialog) {
    dialog.hide();
  }
});
