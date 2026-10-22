// static/generatefinder.js
//
// The admin Generate page's "Around a sector" finder (web/templates/generate_task.html, TODO ADM.4, ADM.109).
//
// Searches sector and star system names through the Galaxy Map's /galaxy/locate, or pages through every filled
// sector from /admin/generate/sectors; picking one fills "Sector ID" and selects "Around a sector" by a filled
// sector. Without JavaScript the Sector ID is typed by hand.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { formatNumber } = await import(`./numberformat.js${VERSION_QUERY}`);

const SEARCH_DELAY_MS = 300;

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
    return item.shell === null || item.shell === undefined
      ? "no address"
      : `shell ${item.shell}, band ${item.band}, wedge ${item.wedge}`;
  }

  function pick(id, text) {
    idInput.value = id;
    for (const [name, value] of [["what", "place"], ["place", "sector"]]) {
      const radio = form.querySelector(`input[name="${name}"][value="${value}"]`);
      if (radio) {
        radio.checked = true;
        radio.dispatchEvent(new Event("change", { bubbles: true }));
      }
    }
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
        ? `${formatNumber(body.total)} filled sector${body.total === 1 ? "" : "s"}, nearest the core first. Page ${page} of ${pages}.`
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
