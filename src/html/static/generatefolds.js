// static/generatefolds.js
//
// The admin Generate page's folding sections and its "Around a sector"
// finder (web/templates/generate.html, TODO ADM.4).
//
// Folds: every section is a <details>. Current job always starts open;
// the others open or close as this browser last left them (localStorage),
// except a section the server marked data-fold-keep (a form shown again
// with its error or estimate), which stays open. A link to a section
// (#plan, or the redirect's #current-job) opens it.
//
// Finder: searches sector and star system names through the Galaxy Map's
// /galaxy/locate, or pages through every filled sector from
// /admin/generate/sectors; picking one fills "Sector ID" and selects
// "Around a sector" by a filled sector. Without JavaScript the sections
// still fold (the server opens Current job and a form shown again), and
// the Sector ID is typed by hand.

const STORE_KEY = "planetgen.generate.open";
const SEARCH_DELAY_MS = 300;

function readStored() {
  try {
    const parsed = JSON.parse(window.localStorage.getItem(STORE_KEY) || "null");
    return Array.isArray(parsed) ? parsed : null;
  } catch {
    return null;
  }
}

function writeStored(ids) {
  try {
    window.localStorage.setItem(STORE_KEY, JSON.stringify(ids));
  } catch {
    // Private mode or a full store: the page just forgets.
  }
}

const folds = Array.from(document.querySelectorAll("details[data-fold]"));

function remembered() {
  return folds.filter((fold) => fold.dataset.fold !== "current-job" && fold.open).map((fold) => fold.dataset.fold);
}

function openForHash() {
  const id = decodeURIComponent(window.location.hash.slice(1));
  if (!id) return;
  const target = document.getElementById(id);
  const fold = target && (target.querySelector("details[data-fold]") || target.closest("details[data-fold]"));
  if (fold && !fold.open) fold.open = true;
}

const stored = readStored();
for (const fold of folds) {
  if (fold.dataset.fold === "current-job" || "foldKeep" in fold.dataset) continue;
  fold.open = stored ? stored.includes(fold.dataset.fold) : false;
}
openForHash();
for (const fold of folds) {
  fold.addEventListener("toggle", () => writeStored(remembered()));
}
window.addEventListener("hashchange", openForHash);

// --- "Around a sector" finder -------------------------------------------------

const finder = document.querySelector("[data-finder]");
if (finder) setUpFinder(finder);

function setUpFinder(box) {
  const form = box.closest("form");
  const query = box.querySelector("[data-finder-query]");
  const status = box.querySelector("[data-finder-status]");
  const results = box.querySelector("[data-finder-results]");
  const pager = box.querySelector("[data-finder-pager]");
  const prev = box.querySelector("[data-finder-prev]");
  const next = box.querySelector("[data-finder-next]");
  const idInput = form.querySelector('input[name="center_sector"]');
  const chosen = form.querySelector("[data-finder-chosen]");
  let page = 1;
  let pages = 1;
  let timer = null;
  let request = 0;

  box.hidden = false;

  function address(item) {
    return item.ring === null || item.ring === undefined
      ? "no address"
      : `ring ${item.ring}, layer ${item.layer}, slot ${item.slot}`;
  }

  function pick(id, text) {
    idInput.value = id;
    const byRadio = form.querySelector('input[name="center_by"][value="sector"]');
    const modeRadio = form.querySelector('input[name="mode"][value="center"]');
    if (byRadio) byRadio.checked = true;
    if (modeRadio) modeRadio.checked = true;
    if (chosen) chosen.textContent = `Center: ${text}.`;
  }

  function show(items, describe) {
    results.replaceChildren(
      ...items.map((item) => {
        const li = document.createElement("li");
        const button = document.createElement("button");
        button.type = "button";
        button.className = "sector-finder-pick";
        const text = describe(item);
        button.textContent = text;
        button.addEventListener("click", () => pick(item.sector_id ?? item.id, text));
        li.append(button);
        return li;
      }),
    );
  }

  async function getJson(url) {
    const response = await fetch(url, { headers: { Accept: "application/json" }, credentials: "same-origin" });
    const body = await response.json().catch(() => ({}));
    if (!response.ok) throw new Error(body.error || `HTTP ${response.status}`);
    return body;
  }

  async function search(term) {
    const mine = ++request;
    pager.hidden = true;
    if (!term) {
      results.replaceChildren();
      status.textContent = "";
      return;
    }
    status.textContent = "Searching…";
    try {
      const body = await getJson(`${box.dataset.locateUrl}?q=${encodeURIComponent(term)}`);
      if (mine !== request) return;
      const matches = body.matches || [];
      status.textContent = matches.length
        ? `${matches.length} match${matches.length === 1 ? "" : "es"}.`
        : "No filled sector or star system has that name.";
      show(matches, (m) =>
        m.kind === "sector"
          ? `Sector ${m.name} (ID ${m.sector_id}, ${address(m)})`
          : `${m.name}, in sector ${m.sector_name} (ID ${m.sector_id}, ${address(m)})`,
      );
    } catch (err) {
      if (mine === request) status.textContent = `The search failed: ${err.message}`;
    }
  }

  async function list(wanted) {
    const mine = ++request;
    status.textContent = "Loading…";
    try {
      const body = await getJson(`${box.dataset.listUrl}?page=${wanted}`);
      if (mine !== request) return;
      page = body.page;
      pages = body.pages;
      status.textContent = body.total
        ? `${body.total.toLocaleString()} filled sector${body.total === 1 ? "" : "s"}, nearest the core first. Page ${page} of ${pages}.`
        : "No sector has been generated yet.";
      show(body.items, (s) => `Sector ${s.name} (ID ${s.id}, ${address(s)}, ${s.systems ?? 0} systems)`);
      pager.hidden = pages <= 1;
      prev.disabled = page <= 1;
      next.disabled = page >= pages;
    } catch (err) {
      if (mine === request) status.textContent = `The list could not be loaded: ${err.message}`;
    }
  }

  query.addEventListener("input", () => {
    clearTimeout(timer);
    timer = setTimeout(() => search(query.value.trim()), SEARCH_DELAY_MS);
  });
  query.addEventListener("keydown", (event) => {
    // Enter searches now instead of submitting the Generate form.
    if (event.key === "Enter") {
      event.preventDefault();
      clearTimeout(timer);
      search(query.value.trim());
    }
  });
  box.querySelector("[data-finder-list]").addEventListener("click", () => {
    query.value = "";
    list(1);
  });
  prev.addEventListener("click", () => list(Math.max(1, page - 1)));
  next.addEventListener("click", () => list(Math.min(pages, page + 1)));
}
