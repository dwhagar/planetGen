// html/static/mapfit.js
//
// The size of a map's window on its page (MAP.291): a map's viewport takes the whole width of its column, and its height
// is what the browser window has left under the things above it (the header, the search) and above the things under it
// (the toolbar), never less than `min`. So the shape of the viewport follows the shape of the window (wide in a wide
// window, tall in a tall one, with no fixed aspect ratio) and the renderer follows the viewport (mapcore.js
// watchResize). `fitHeight` is the pure part; `watchWindowFit` keeps one viewport fitted. Any page that gives a map a
// window of its own can use it (UX.142: every viewer scales to its window). Tests: tests/js/mapfit.test.mjs.

export const MIN_HEIGHT_PX = 320;
export const GAP_PX = 16;

// The height in pixels a viewport whose top edge lies `top` pixels down the page (not the window) may take in a window
// `windowHeight` tall when `below` pixels of toolbar follow it: the rest of the window, kept at or above `min`.
export function fitHeight(windowHeight, top, below, options) {
  const min = options && options.min !== undefined ? options.min : MIN_HEIGHT_PX;
  const gap = options && options.gap !== undefined ? options.gap : GAP_PX;
  const room = windowHeight - top - below - gap;
  return Math.max(min, Math.round(Number.isFinite(room) ? room : min));
}

// Sets `view.style.height` to the fitting height now, when the window is resized and when the page above the view
// changes height (a notice, wrapped crumbs). `below` is the element after the view (or null). Returns {refit, stop}.
export function watchWindowFit(view, below, options) {
  let busy = false;
  function refit() {
    if (busy) return;
    busy = true;
    try {
      const rect = view.getBoundingClientRect();
      const top = rect.top + (window.scrollY || 0);
      const under = below ? below.getBoundingClientRect().height : 0;
      const wanted = fitHeight(window.innerHeight, top, under, options) + "px";
      if (view.style.height !== wanted) view.style.height = wanted;
    } finally {
      busy = false;
    }
  }
  let observer = null;
  const page = view.closest ? view.closest("main, .app-main") : null;
  if (typeof ResizeObserver !== "undefined") {
    observer = new ResizeObserver(refit);
    if (page) observer.observe(page);
    if (below) observer.observe(below);
  }
  window.addEventListener("resize", refit);
  refit();
  return {
    refit: refit,
    stop: function () {
      if (observer) observer.disconnect();
      window.removeEventListener("resize", refit);
    },
  };
}
