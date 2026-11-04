// html/static/levelprogress.js
//
// MAP.364: the progress bars the Galaxy Map shows while the nodes of a map level load. A request names the nodes it asks for
// (`m.shell.band.wedge` keys); a node's level is its `m` (3 shells deep, a sector at 1). Each level with nodes still to come gets
// one bar: nodes received over nodes requested for that level. A bar shows only if the load is still going 400 ms after it began, so a
// level served from this browser's store or the server's cache never flashes one. A load that fails, or a node that arrives marked
// incomplete, turns its level's bar into a line in plain words that says it will try again. Every bar is a `role="progressbar"` named
// for its level; the colors are the map viewport's own tokens, so both themes read.

export const SHOW_DELAY_MS = 400;
export const LINGER_MS = 250;

// The level of a node key: its first number.
export function levelOf(key) {
  const m = Number(String(key).split(".")[0]);
  return Number.isFinite(m) && m > 0 ? m : null;
}

// What the map calls a level, in plain words.
export function levelName(m) {
  return m === 1 ? "sectors" : "blocks " + m + " shells deep";
}

export function createLevelProgress(host, options) {
  const doc = host.ownerDocument || document;
  const timers = options && options.timers ? options.timers : { set: (fn, ms) => setTimeout(fn, ms), clear: (handle) => clearTimeout(handle) };
  let box = null;
  let levels = new Map();       // m -> { total, done: Set, stalled: text or null, row }
  let showTimer = null;
  let removeTimer = null;
  let active = false;

  function ensureBox() {
    if (box) return;
    box = doc.createElement("div");
    box.className = "map-level-progress";
    host.appendChild(box);
  }

  function paint() {
    if (!box) return;
    levels.forEach((level, m) => {
      if (level.hidden) return;
      if (!level.row) {
        level.row = doc.createElement("div");
        level.row.className = "map-level-row";
        box.appendChild(level.row);
      }
      const row = level.row;
      if (level.stalled) {
        row.className = "map-level-row map-level-stalled";
        row.setAttribute("role", "status");
        row.removeAttribute("aria-valuenow");
        row.removeAttribute("aria-valuemax");
        row.removeAttribute("aria-valuemin");
        row.removeAttribute("aria-label");
        row.textContent = level.stalled;
        return;
      }
      if (row.getAttribute("role") !== "progressbar") {
        row.textContent = "";
        row.className = "map-level-row";
        row.setAttribute("role", "progressbar");
        row.setAttribute("aria-label", "Loading " + levelName(m));
        row.setAttribute("aria-valuemin", "0");
        const label = doc.createElement("span");
        label.className = "map-level-label";
        const track = doc.createElement("span");
        track.className = "map-level-track";
        const fill = doc.createElement("span");
        fill.className = "map-level-fill";
        track.appendChild(fill);
        row.appendChild(label);
        row.appendChild(track);
      }
      const have = level.done.size;
      row.setAttribute("aria-valuemax", String(level.total));
      row.setAttribute("aria-valuenow", String(have));
      row.setAttribute("aria-valuetext", have + " of " + level.total + " " + levelName(m));
      row.children[0].textContent = levelName(m).charAt(0).toUpperCase() + levelName(m).slice(1) + ": " + have + " of " + level.total;
      row.children[1].children[0].style.width = (level.total ? Math.round((100 * have) / level.total) : 0) + "%";
    });
  }

  function clearBars() {
    if (box) {
      box.remove();
      box = null;
    }
    levels.forEach((level) => { level.row = null; });
  }

  function show() {
    showTimer = null;
    if (!active) return;
    ensureBox();
    paint();
  }

  function cancelRemoval() {
    if (removeTimer !== null) {
      timers.clear(removeTimer);
      removeTimer = null;
    }
  }

  return {
    // A request for these node keys begins. A request that follows one still showing keeps its bars, so a map that goes again does not blink.
    start(keys) {
      cancelRemoval();
      const wasShowing = !!box;
      levels = new Map();
      keys.forEach((key) => {
        const m = levelOf(key);
        if (m === null) return;
        if (!levels.has(m)) levels.set(m, { total: 0, done: new Set(), stalled: null, row: null, keys: new Set() });
        const level = levels.get(m);
        if (!level.keys.has(key)) {
          level.keys.add(key);
          level.total += 1;
        }
      });
      if (box) {
        clearBars();
      }
      active = levels.size > 0;
      if (showTimer !== null) {
        timers.clear(showTimer);
        showTimer = null;
      }
      if (active) {
        if (wasShowing) show();
        else showTimer = timers.set(show, SHOW_DELAY_MS);
      }
    },
    // These keys arrived whole.
    arrived(keys) {
      keys.forEach((key) => {
        const level = levels.get(levelOf(key));
        if (level && level.keys.has(key)) level.done.add(key);
      });
      if (box) paint();
    },
    // The load did not finish: every level with nodes still to come says so. `retrySeconds` is the pause before the next try.
    stalled(retrySeconds) {
      if (showTimer !== null) {
        timers.clear(showTimer);
        showTimer = null;
      }
      let any = false;
      levels.forEach((level, m) => {
        if (level.done.size >= level.total) {
          if (level.row) level.row.remove();      // a complete level has nothing to say
          level.row = null;
          level.hidden = true;
          return;
        }
        any = true;
        level.stalled = "Could not load all of the " + levelName(m) + (retrySeconds ? ". Trying again in " + retrySeconds + (retrySeconds === 1 ? " second." : " seconds.") : ". Trying again.");
      });
      if (any) {
        active = true;
        ensureBox();
        paint();
      } else {
        this.finish();
      }
    },
    // Nodes that came back marked incomplete: their levels say the rest is still being worked out.
    incomplete(keys) {
      keys.forEach((key) => {
        const level = levels.get(levelOf(key));
        if (level) level.stalled = "Some of the " + levelName(levelOf(key)) + " are still being worked out. Trying again soon.";
      });
      if (keys.length) {
        ensureBox();
        paint();
      }
    },
    // The request is over and every level is complete: the bars go after a moment.
    finish() {
      if (showTimer !== null) {
        timers.clear(showTimer);
        showTimer = null;
      }
      active = false;
      cancelRemoval();
      if (!box) {
        levels = new Map();
        return;
      }
      // A level that still says it is stalled (an incomplete node) keeps its line until the next request begins.
      if (Array.from(levels.values()).some((level) => level.stalled)) {
        return;
      }
      levels.forEach((level) => { level.done = new Set(level.keys); });
      paint();
      removeTimer = timers.set(() => {
        removeTimer = null;
        clearBars();
        levels = new Map();
      }, LINGER_MS);
    },
    // A request left behind (aborted for a newer one): nothing to show for it.
    abandon() {
      if (showTimer !== null) {
        timers.clear(showTimer);
        showTimer = null;
      }
      active = false;
    },
    // For tests and the page's own hooks.
    state() {
      return {
        shown: !!box,
        levels: Array.from(levels.entries()).map(([m, level]) => ({ m, total: level.total, done: level.done.size, stalled: level.stalled })),
      };
    },
  };
}
