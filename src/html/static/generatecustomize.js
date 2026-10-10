// static/generatecustomize.js
//
// The admin Generate page's "Customize" window (TODO ADM.28): a form marked
// data-customize holds its less common settings in fieldsets marked
// data-customize-tab="Name". Without JavaScript they simply stack in the
// form. With it they move into a window opened by a Customize button, one
// tab per fieldset, so the page keeps only the common actions. The inputs
// stay inside the form (the window is a child of it), so closing the window
// changes nothing about what is submitted.

export const TAB_ATTR = "data-customize-tab";

function make(tag, className, text) {
  const el = document.createElement(tag);
  if (className) el.className = className;
  if (text) el.textContent = text;
  return el;
}

function open(dialog) {
  if (typeof dialog.showModal === "function") dialog.showModal();
  else dialog.setAttribute("open", "");
}

function close(dialog) {
  if (typeof dialog.close === "function") dialog.close();
  else dialog.removeAttribute("open");
}

// Builds the window for `form` and returns {dialog, button, select(i)}, or
// null when the form has no settings groups.
export function enhance(form) {
  const groups = Array.from(form.querySelectorAll("[" + TAB_ATTR + "]"));
  if (!groups.length) return null;
  const dialog = make("dialog", "customize-dialog");
  dialog.setAttribute("aria-label", "Customize");
  const tablist = make("div", "customize-tabs");
  tablist.setAttribute("role", "tablist");
  const panels = make("div", "customize-panels");
  const tabs = groups.map((group, index) => {
    const tab = make("button", "customize-tab", group.getAttribute(TAB_ATTR));
    tab.type = "button";
    tab.setAttribute("role", "tab");
    tab.addEventListener("click", () => select(index));
    tablist.appendChild(tab);
    group.classList.add("customize-group");
    panels.appendChild(group);
    return tab;
  });

  function select(chosen) {
    groups.forEach((group, index) => {
      group.hidden = index !== chosen;
      tabs[index].setAttribute("aria-selected", index === chosen ? "true" : "false");
      tabs[index].classList.toggle("selected", index === chosen);
    });
  }

  const done = make("button", "btn", "Done");
  done.type = "button";
  done.addEventListener("click", () => close(dialog));
  const footer = make("div", "search-actions");
  footer.appendChild(done);
  dialog.appendChild(make("h3", "customize-title", "Customize"));
  dialog.appendChild(tablist);
  dialog.appendChild(panels);
  dialog.appendChild(footer);
  select(0);

  const button = make("button", "btn", "Customize");
  button.type = "button";
  button.addEventListener("click", () => open(dialog));
  form.appendChild(dialog);
  // A value the browser rejects (a share over 100, say) must not hide in a
  // closed window: open it on that tab so the message can show.
  form.addEventListener("invalid", (event) => {
    const index = groups.findIndex((group) => group.contains(event.target));
    if (index < 0) return;
    select(index);
    if (!dialog.open && !dialog.hasAttribute("open")) open(dialog);
  }, true);
  const actions = Array.from(form.children).find((child) => child.classList.contains("search-actions"));
  if (actions) actions.insertBefore(button, actions.childNodes[0] || null);
  else form.appendChild(button);
  return { dialog, button, select };
}

export function init() {
  return Array.from(document.querySelectorAll("form[data-customize]")).map(enhance).filter(Boolean);
}

init();
