// html/static/galaxyregionlist.js
//
// MAP.316: the Menu's list of charted regions, the always-available way to reach known space. Each region is a button
// ("Fardoo Beez Region, 14 percent charted"); the list sorts by share charted or by name (population joins once it is
// stored, MAP.318), and is empty, hidden, while no region is charted. The map says what a button does (`onOpen(blob)`).

const VERSION_QUERY = new URL(import.meta.url).search;
const { sortRegions, chartedText } = await import(`./galaxyregions.js${VERSION_QUERY}`);

// The list in `root` (the `<details>` of the Menu, holding `[data-regions-count]`, `[data-regions-sort]` and
// `[data-regions-list]`): `set(blobs)` replaces the regions, the sort select re-orders them.
export function createRegionList(root, onOpen) {
  const list = root.querySelector("[data-regions-list]");
  const select = root.querySelector("[data-regions-sort]");
  const count = root.querySelector("[data-regions-count]");
  let blobs = [];

  function render() {
    const by = select ? select.value : "share";
    list.textContent = "";
    sortRegions(blobs, by).forEach((blob) => {
      const item = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.className = "galaxy-region-button";
      button.dataset.regionKey = blob.key;
      const title = blob.name || "Region " + blob.key;
      button.textContent = title + ", " + chartedText(blob.share);
      button.title = blob.sectorCount + " of " + blob.capacity + " sectors generated" + (blob.parent ? " (in " + blob.parent + ")" : "");
      button.addEventListener("click", () => onOpen(blob));
      item.appendChild(button);
      list.appendChild(item);
    });
    if (count) count.textContent = blobs.length ? "(" + blobs.length + ")" : "";
    root.hidden = blobs.length === 0;
  }

  if (select) select.addEventListener("change", render);
  return {
    set(next) { blobs = next.slice(); render(); },
    // For the browser tests: the names in the order shown.
    order() { return Array.from(list.querySelectorAll("button")).map((b) => b.textContent); },
  };
}
