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
// MAP.350: a map that fills the window takes this share of the window's height (what a held header leaves of it); the rest
// is split evenly above and below it. Boss: always 95% of the window, however the window is sized.
export const FILL_SHARE = 0.95;
// How far (pixels) the page may be from the aligned place and still count as there, so a rounding does not scroll it again and again.
const ALIGN_TOLERANCE_PX = 2;

// The height in pixels a viewport whose top edge lies `top` pixels down the page (not the window) may take in a window
// `windowHeight` tall when `below` pixels of toolbar follow it: the rest of the window, kept at or above `min`.
export function fitHeight(windowHeight, top, below, options) {
  const min = options && options.min !== undefined ? options.min : MIN_HEIGHT_PX;
  const gap = options && options.gap !== undefined ? options.gap : GAP_PX;
  const room = windowHeight - top - below - gap;
  return Math.max(min, Math.round(Number.isFinite(room) ? room : min));
}

// MAP.350: the height of a map that fills its browser window: `FILL_SHARE` (95%) of the window's height less what a header held at the top
// of the window covers (`options.covered`, pixels), kept at or above `min`.
export function fillHeight(windowHeight, options) {
  const min = options && options.min !== undefined ? options.min : MIN_HEIGHT_PX;
  const share = options && options.share !== undefined ? options.share : FILL_SHARE;
  const covered = options && Number.isFinite(options.covered) ? options.covered : 0;
  const room = (windowHeight - covered) * share;
  return Math.max(min, Math.round(Number.isFinite(room) ? room : min));
}

// The margin above (and below) a map of `height` pixels in a window `windowHeight` tall under `covered` pixels of header: half of what the map leaves.
export function fillBuffer(windowHeight, height, covered) {
  const room = (windowHeight - (Number.isFinite(covered) ? covered : 0) - height) / 2;
  return Number.isFinite(room) ? Math.max(0, Math.round(room)) : 0;
}

// MAP.350: where to scroll the page (pixels down it) so a map whose top edge is `top` pixels down the page has `buffer` pixels of window
// above it under what a header held at the top of the window covers, or null when it is there already (or the page cannot be scrolled that far; `maxScroll` is the farthest it goes).
export function alignScroll(top, scrollY, maxScroll, options) {
  const buffer = options && options.buffer !== undefined ? options.buffer : 0;
  if (![top, scrollY, maxScroll].every(Number.isFinite)) return null;
  const covered = options && Number.isFinite(options.covered) ? options.covered : 0;
  const wanted = Math.max(0, Math.min(maxScroll, Math.round(top - covered - buffer)));
  return Math.abs(wanted - scrollY) <= ALIGN_TOLERANCE_PX ? null : wanted;
}

// Sets `view.style.height` to the fitting height now, when the window is resized and when the page above the view
// changes height (a notice, wrapped crumbs). `below` is the element after the view (or null). Returns {refit, stop}.
// MAP.350: with `options.fill` the view is 95% of the window's height (less what a held header covers), with the rest split above and below it, and the page is scrolled
// so it sits there, once it has loaded and as what is above it settles, until the reader scrolls the page themselves.
export function watchWindowFit(view, below, options) {
  let busy = false;
  const fill = !!(options && options.fill);
  let readerScrolled = false;
  let ourScroll = false;
  // The height of the header the page holds at the top of the window (position fixed or sticky), which the map is not to go under.
  function covered() {
    const header = document.querySelector("header");
    if (!header) return 0;
    const position = window.getComputedStyle(header).position;
    return position === "fixed" || position === "sticky" ? Math.round(header.getBoundingClientRect().bottom) : 0;
  }
  function align() {
    if (!fill || readerScrolled) return;
    const doc = document.documentElement;
    const rect = view.getBoundingClientRect();
    const hidden = covered();
    const buffer = fillBuffer(window.innerHeight, rect.height, hidden);
    const to = alignScroll(rect.top + (window.scrollY || 0), window.scrollY || 0, Math.max(0, doc.scrollHeight - window.innerHeight), Object.assign({}, options, { covered: hidden, buffer: buffer }));
    if (to === null) return;
    ourScroll = true;
    window.scrollTo(0, to);
    window.requestAnimationFrame(function () { ourScroll = false; });
  }
  function refit() {
    if (busy) return;
    busy = true;
    try {
      const rect = view.getBoundingClientRect();
      const top = rect.top + (window.scrollY || 0);
      const under = below ? below.getBoundingClientRect().height : 0;
      const height = fill ? fillHeight(window.innerHeight, Object.assign({}, options, { covered: covered() })) : fitHeight(window.innerHeight, top, under, options);
      const wanted = height + "px";
      if (view.style.height !== wanted) view.style.height = wanted;
      align();
    } finally {
      busy = false;
    }
  }
  function onScroll() {
    if (!ourScroll) readerScrolled = true;
  }
  if (fill) window.addEventListener("scroll", onScroll, { passive: true });
  if (fill) window.addEventListener("load", refit);
  let observer = null;
  const page = view.closest ? view.closest("main, .app-main") : null;
  if (typeof ResizeObserver !== "undefined") {
    observer = new ResizeObserver(refit);
    if (page) observer.observe(page);
    if (below) observer.observe(below);
    if (fill && document.querySelector("header")) observer.observe(document.querySelector("header"));
  }
  window.addEventListener("resize", refit);
  refit();
  return {
    refit: refit,
    stop: function () {
      if (observer) observer.disconnect();
      window.removeEventListener("resize", refit);
      window.removeEventListener("scroll", onScroll);
      window.removeEventListener("load", refit);
    },
  };
}
