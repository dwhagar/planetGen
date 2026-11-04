// html/static/boxfit.js
//
// UX.142: a box of fixed shape on a page (the nav chart, the system diagram) stays inside the browser window: it is never taller than what
// the window leaves under the held header (`mapfit.fitWidth`), however wide its column is, and re-fits when the window is resized. A box asks
// for it with `data-fit-aspect="W/H"` (its width over its height, the `aspect-ratio` its CSS gives it).

const VERSION_QUERY = new URL(import.meta.url).search;
const { watchWindowBox } = await import(`./mapfit.js${VERSION_QUERY}`);

for (const box of document.querySelectorAll("[data-fit-aspect]")) {
  const [width, height] = box.dataset.fitAspect.split("/").map(Number);
  if (width > 0 && height > 0) watchWindowBox(box, width / height);
}
