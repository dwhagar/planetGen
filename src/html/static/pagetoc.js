// UX.196: the contents list of a page, built from its headings. Every h2 and h3 of <main> becomes a link in a
// "On this page" navigation after the page title (and its sub-head), nested by level. A heading with no id gets one.
// A link to a heading inside a collapsed section opens that section (opensection.js). Pages with fewer than three
// headings, and headings in dialogs or hidden regions, get no list.
(function () {
  const main = document.getElementById("main");
  if (!main) return;
  const heads = [...main.querySelectorAll("h2, h3")].filter(
    h => !h.closest("sl-dialog, dialog, [hidden], .page-toc") && (!h.checkVisibility || h.checkVisibility()));
  if (heads.length < 3) return;

  const used = new Set([...document.querySelectorAll("[id]")].map(el => el.id));
  function idOf(h) {
    if (h.id) return h.id;
    const base = "toc-" + (h.textContent.trim().toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "section");
    let id = base, n = 2;
    while (used.has(id)) id = base + "-" + n++;
    used.add(id);
    h.id = id;
    return id;
  }

  const nav = document.createElement("nav");
  nav.className = "page-toc";
  nav.setAttribute("aria-label", "On this page");
  const title = document.createElement("p");
  title.className = "page-toc-title";
  title.textContent = "On this page";
  nav.append(title);
  const list = document.createElement("ul");
  let sub = null;
  for (const h of heads) {
    const item = document.createElement("li");
    const link = document.createElement("a");
    link.href = "#" + idOf(h);
    link.textContent = h.textContent.trim();
    item.append(link);
    if (h.tagName === "H3" && list.lastElementChild) {
      if (!sub || sub.parentElement !== list.lastElementChild) {
        sub = document.createElement("ul");
        list.lastElementChild.append(sub);
      }
      sub.append(item);
    } else {
      list.append(item);
    }
  }
  nav.append(list);

  const anchor = main.querySelector(".page-subhead") || main.querySelector("h1");
  if (anchor) anchor.after(nav);
})();
