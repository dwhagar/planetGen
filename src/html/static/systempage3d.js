// html/static/systempage3d.js
//
// The system page's switch between the Diagram (systemmap.js, the flat SVG)
// and the 3D view (systemview3d.js) (MAP.74). The page always loads the
// diagram; the 3D view's script and the system's scene are fetched the first
// time 3D is chosen, so a visitor who never switches pays nothing for it. The
// last view used is remembered in the browser; `?view=3d` and
// `?object=<ref>` (a star, planet, moon or comet, e.g. `planet:12`) open the
// 3D view and select a body, and the address follows the view (MAP.67) so
// Back, reload and bookmarks return to it. Without WebGL the switch says so
// and the diagram stays. A list of the bodies stands in for the canvas for
// screen readers.

const VERSION_QUERY = new URL(import.meta.url).search;
const STORAGE_KEY = "planetgen.systemView";

const root = document.getElementById("sysview3d");
const diagram = document.getElementById("sysmap-diagram");
const diagramButton = document.getElementById("sysmap-view-diagram");
const button3d = document.getElementById("sysmap-view-3d");

function remembered() {
  try {
    return window.localStorage.getItem(STORAGE_KEY);
  } catch (err) {
    return null;
  }
}

function remember(view) {
  try {
    window.localStorage.setItem(STORAGE_KEY, view);
  } catch (err) {
    // A private window: the choice just isn't kept.
  }
}

function hasWebGl() {
  try {
    const canvas = document.createElement("canvas");
    return !!(window.WebGLRenderingContext && (canvas.getContext("webgl2") || canvas.getContext("webgl")));
  } catch (err) {
    return false;
  }
}

function setUrl(view, ref) {
  const url = new URL(window.location.href);
  if (view === "3d") url.searchParams.set("view", "3d");
  else url.searchParams.delete("view");
  if (view === "3d" && ref) url.searchParams.set("object", ref);
  else url.searchParams.delete("object");
  window.history.replaceState(window.history.state, "", url.pathname + url.search + url.hash);
}

if (root && diagram && diagramButton && button3d) {
  let viewModule = null;
  let clockModule = null;
  let view = null;
  let clock = null;
  let loading = null;
  let current = "diagram";
  const els = {
    canvas: document.getElementById("sysview3d-canvas"),
    labels: document.getElementById("sysview3d-labels"),
    tip: document.getElementById("sysview3d-tip"),
    info: document.getElementById("sysmap-info"),
    play: document.getElementById("sysview3d-play"),
    slower: document.getElementById("sysview3d-slower"),
    faster: document.getElementById("sysview3d-faster"),
    now: document.getElementById("sysview3d-now"),
    rate: document.getElementById("sysview3d-rate"),
    scale: document.getElementById("sysview3d-scale"),
    scaleNote: document.getElementById("sysview3d-scale-note"),
    reset: document.getElementById("sysview3d-reset"),
    list: document.getElementById("sysview3d-list"),
  };

  function showRate() {
    els.play.textContent = clock.playing() ? "Pause" : "Play";
    els.play.setAttribute("aria-pressed", clock.playing() ? "true" : "false");
    els.rate.textContent = clock.playing() ? clock.rate().label : "Paused";
  }

  function fillList(list) {
    els.list.textContent = "";
    list.forEach((item) => {
      const li = document.createElement("li");
      const open = document.createElement("button");
      open.type = "button";
      open.className = "btn btn-small btn-secondary";
      open.textContent = item.name + " (" + item.kind + ")";
      open.addEventListener("click", () => view.select(item.ref));
      li.appendChild(open);
      els.list.appendChild(li);
    });
  }

  async function start() {
    if (view) return view;
    if (loading) return loading;
    loading = (async () => {
      const response = await fetch(root.dataset.sceneUrl, { headers: { Accept: "application/json" } });
      if (!response.ok) throw new Error("scene " + response.status);
      const scene = await response.json();
      viewModule = await import(`./systemview3d.js${VERSION_QUERY}`);
      clockModule = await import(`./orbitclock.js${VERSION_QUERY}`);
      clock = clockModule.createOrbitClock({ epochUnix: scene.epoch_unix, nowMs: Date.now() });
      view = viewModule.createSystemView({
        canvasEl: els.canvas, labelsEl: els.labels, tooltipEl: els.tip, infoEl: els.info, scene: scene, clock: clock,
        onSelect: (entry) => { if (current === "3d") setUrl("3d", entry ? entry.ref : null); },
      });
      els.scaleNote.textContent = view.note();
      fillList(view.list());
      showRate();
      root.dataset.ready = "true";
      return view;
    })();
    return loading;
  }

  function show(which) {
    current = which;
    const is3d = which === "3d";
    diagram.hidden = is3d;
    root.hidden = !is3d;
    diagramButton.setAttribute("aria-pressed", is3d ? "false" : "true");
    button3d.setAttribute("aria-pressed", is3d ? "true" : "false");
    remember(which);
    setUrl(which, is3d && view && view.selected() ? view.selected().ref : null);
  }

  // Without WebGL the switch says so and the diagram stays.
  function noWebGl() {
    show("diagram");
    button3d.disabled = true;
    button3d.title = "3D needs WebGL, which this browser does not have.";
    if (!document.getElementById("sysview3d-nogl")) {
      const note = document.createElement("p");
      note.id = "sysview3d-nogl";
      note.className = "hint";
      note.textContent = "3D needs WebGL, which this browser does not have. The diagram shows the same system.";
      button3d.parentElement.after(note);
    }
  }

  async function choose3d(ref) {
    if (!hasWebGl()) {
      noWebGl();
      return;
    }
    show("3d");
    try {
      await start();
    } catch (err) {
      show("diagram");
      button3d.disabled = true;
      button3d.title = "The 3D view could not be loaded.";
      return;
    }
    if (ref) {
      view.select(ref);
      view.flyTo(ref);
    }
  }

  diagramButton.addEventListener("click", () => show("diagram"));
  button3d.addEventListener("click", () => choose3d(null));
  els.play.addEventListener("click", () => { if (clock) { clock.toggle(); showRate(); } });
  els.slower.addEventListener("click", () => { if (clock) { clock.slower(); showRate(); } });
  els.faster.addEventListener("click", () => { if (clock) { clock.faster(); showRate(); } });
  els.now.addEventListener("click", () => { if (clock) { clock.now(Date.now()); showRate(); } });
  els.reset.addEventListener("click", () => { if (view) view.reset(); });
  els.scale.addEventListener("change", () => {
    if (!view) return;
    view.setMode(els.scale.value);
    els.scaleNote.textContent = view.note();
  });

  const params = new URLSearchParams(window.location.search);
  if (params.get("view") === "3d" || params.get("object") || remembered() === "3d") {
    choose3d(params.get("object"));
  }
}
