// UX.196: the contents list of a page, built from its headings. Every h2 and h3 of <main> becomes a link in a
// "On this page" navigation after the page title (and its sub-head), nested by level. A heading with no id gets one.
// A link to a heading inside a collapsed section opens that section (opensection.js). Pages with fewer than three
// headings, and headings in dialogs or hidden regions, get no list. A page that already has a contents list of its own gets no
// second one (UX.235): a nav in <main>, or any element flagged `data-own-toc`, with three or more links to headings or sections on
// the same page (the Settings page's group list is one); a menu of links to other pages is not. The Galaxy Map's page (the wide
// main column, `app-main-wide`) gets none: its map takes the page's whole width, and the menu's column would narrow it (TEST.442).
(function () {
  const main = document.getElementById("main");
  if (!main || main.classList.contains("app-main-wide")) return;
  const heads = [...main.querySelectorAll("h2, h3")].filter(
    h => !h.closest("sl-dialog, dialog, [hidden], .page-toc") && (!h.checkVisibility || h.checkVisibility()));
  if (heads.length < 3) return;
  const hasOwnContents = [...main.querySelectorAll("nav, [data-own-toc]")].some(el => {
    if (el.closest(".page-toc")) return false;
    const inPage = [...el.querySelectorAll("a[href^='#']")].filter(a => a.getAttribute("href").length > 1
      && document.getElementById(a.getAttribute("href").slice(1)));
    return inPage.length >= 3;
  });
  if (hasOwnContents) return;

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
  main.classList.add("has-page-toc");

  // If every link fits in one row across the top, show the row; if any would need scrolling, use the grouped menu.
  function fit() {
    main.classList.add("toc-strip");
    if (nav.scrollWidth > nav.clientWidth + 1) main.classList.remove("toc-strip");
  }
  fit();
  window.addEventListener("resize", fit);
  window.addEventListener("load", fit);
})();
