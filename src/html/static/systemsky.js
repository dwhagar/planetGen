// html/static/systemsky.js
//
// VIEW.35 (S17 of the System Sky plan, reports/system-sky-from-the-system-2026-11-04.md): the System View's background becomes the sky its system sees.
// This is the page's controller for it. It asks the sky route (VIEW.31, docs/design/interface/system-sky-wire-format.md) for the packed stars and the
// nebula and remnant sprites, asks the band route (VIEW.13) for the glow image the sky names, builds the three layers (skyband.js VIEW.34, skysprites.js
// VIEW.36, skysphere.js VIEW.33), turns them from the galaxy's axes into the scene's and hands them to the view (`view.attachSky`), which updates them
// every frame from the camera and the orbit clock.
//
//   * A preference switches the sky on and off (the browser remembers it; on when none is) and a status line says what is shown.
//   * A sky whose tiers are still being filled in (`meta.report.filling`) is drawn as far as it is built, says so, and is asked for again every
//     POLL_MS (at most MAX_POLLS times) until it is complete; the new sky replaces the old.
//   * A sky that cannot be had (no place in the galaxy, a failed request) leaves the old backdrop, the view's own background, and says why. A browser
//     without WebGL never gets this far: the 3D view itself is not offered (systempage3d.js).
//
// The layers test depth against the system's bodies (the sky is far past them) so a planet in front of a star hides it, which the layers' own tests, with
// no bodies, do not need. A system's orbital plane is not tied to the galactic frame in the stored model, so the galactic plane is laid in the scene's
// reference plane: the galaxy's north pole (z) is the scene's up (y), the galactic x stays x. A stand-in, replaced when a system gets an orientation.
// `THREE` is passed in, so plain node can test the rest.

const VERSION_QUERY = new URL(import.meta.url).search;
const { unpackSky } = await import(`./skypack.js${VERSION_QUERY}`);
const { formatNumber } = await import(`./numberformat.js${VERSION_QUERY}`);
const { createSkySphere, skyYears } = await import(`./skysphere.js${VERSION_QUERY}`);
const { createSkySprites } = await import(`./skysprites.js${VERSION_QUERY}`);
const { createSkyBand, imageFromHalfFloats } = await import(`./skyband.js${VERSION_QUERY}`);
const { skyLimit } = await import(`./skylight.js${VERSION_QUERY}`);
const { pickSkyStar, describeSkyStar, skyStarsOnScreen } = await import(`./skypick.js${VERSION_QUERY}`);

export const STORAGE_KEY = "planetgen.systemSky";
// A sky still filling in is asked for again this often (ms), this many times.
export const POLL_MS = 30000;
export const MAX_POLLS = 20;
// A window this narrow or narrower is a phone: a lower limiting magnitude and fewer stars.
export const PHONE_MAX_WIDTH_PX = 700;

// The rotation that turns the galaxy's axes into the scene's: a quarter turn about x, (x, y, z) -> (x, z, -y), as the quaternion [x, y, z, w].
export const GALAXY_TO_SCENE = [-Math.SQRT1_2, 0, 0, Math.SQRT1_2];

// The preference: `true` (on) unless the browser remembers "off". `storage` is a Storage or null; every access is guarded.
export function readPreference(storage) {
  try {
    return !storage || storage.getItem(STORAGE_KEY) !== "off";
  } catch (err) {
    return true;
  }
}

export function writePreference(storage, on) {
  try {
    if (storage) storage.setItem(STORAGE_KEY, on ? "on" : "off");
  } catch (err) {
    // A private window: the choice just isn't kept.
  }
}

export function isPhone(win) {
  try {
    return !!win && (win.innerWidth || Infinity) <= PHONE_MAX_WIDTH_PX;
  } catch (err) {
    return false;
  }
}

// The sky route's URL for a limit: `base` is the page's `data-sky-url`, which may already carry a query.
export function skyUrlFor(base, limit) {
  return base + (base.indexOf("?") >= 0 ? "&" : "?") + "limit=" + encodeURIComponent(String(limit));
}

// The band route's URL: `meta.band.href` already holds `?limit=`; the key asks for the stored copy by name.
export function bandUrlFor(band) {
  return band.href + (band.href.indexOf("?") >= 0 ? "&" : "?") + "key=" + encodeURIComponent(band.key);
}

function parsecs(pc) {
  return pc >= 1000 ? (pc / 1000).toFixed(1) + " kpc" : Math.round(pc) + " pc";
}

// The status line for a state: `{ kind: "off" | "loading" | "ready" | "filling" | "unavailable", report, count, status }`.
export function noteFor(state) {
  switch (state.kind) {
    case "off": return "Sky background off.";
    case "loading": return "Loading the sky this system sees.";
    case "unavailable":
      return state.status === 404 ? "This system has no place in the galaxy yet, so it has no sky to show." : "The sky could not be loaded; the plain background is shown.";
    case "filling": {
      const report = state.report || {};
      const built = Number.isFinite(report.built_to_pc) ? ", built to " + parsecs(report.built_to_pc) : "";
      return "Sky of " + formatNumber(state.count) + " stars" + built + ", filling in.";
    }
    case "ready": return "Sky of " + formatNumber(state.count) + " stars, with the glow of the fainter ones.";
    default: return "";
  }
}

// Every material of the layer's objects tests depth against the bodies.
function testDepth(objects) {
  objects.forEach((object) => { object.material.depthTest = true; });
}

// Builds the layers of an unpacked sky and joins them for the view: `{ objects, update(camera, state), dispose(), count }`. `state` is `{ years,
// pixelRatio }`, `years` the orbit clock's years after the scene's epoch.
export function buildLayer(THREE, sky, band, options) {
  const o = options || {};
  const sphere = createSkySphere(THREE, sky, { phone: o.phone, limit: o.limit, pixelRatio: o.pixelRatio || 1 });
  const sprites = createSkySprites(THREE, sky.meta.nebulae || [], {});
  const glow = band ? createSkyBand(THREE, band, { seed: o.seed }) : null;
  const objects = [glow && glow.object, ...sprites.objects, sphere.object].filter(Boolean);
  testDepth(objects);
  const turn = new THREE.Quaternion(...GALAXY_TO_SCENE);
  objects.forEach((object) => object.quaternion.copy(turn));
  let lastYears = 0;
  // The camera's view of the sky for skypick.js: the galaxy's axes into the camera's (3 x 3, row-major) and the projection's scales.
  function viewOf(camera, size) {
    camera.updateMatrixWorld();
    const rotation = new THREE.Matrix4().extractRotation(camera.matrixWorldInverse)
      .multiply(new THREE.Matrix4().makeRotationFromQuaternion(turn));
    const e = rotation.elements;
    return {
      matrix: [e[0], e[4], e[8], e[1], e[5], e[9], e[2], e[6], e[10]],
      p00: camera.projectionMatrix.elements[0], p11: camera.projectionMatrix.elements[5], width: size.width, height: size.height,
    };
  }
  return {
    objects, count: sphere.count, sphere, sprites, band: glow,
    // VIEW.37: the star under a point of the canvas (`pointer` `{ x, y }` CSS pixels from its top left), as `{ title, lines, href }`, or null.
    // `size` is `{ width, height }` of the canvas.
    screenStars(camera, size) {
      return skyStarsOnScreen(sphere.kept, lastYears, viewOf(camera, size), { limit: sphere.limit() });
    },
    pick(camera, pointer, size) {
      const view = viewOf(camera, size);
      const found = pickSkyStar(sphere.kept, lastYears, view, pointer, { limit: sphere.limit() });
      return found ? Object.assign({ pick: found }, describeSkyStar(found, sky.meta && sky.meta.charted, o.starUrl)) : null;
    },
    update(camera, state) {
      const years = skyYears(sky, o.sceneEpochUnix, state.years);
      lastYears = years;
      sphere.update(camera, { years, pixelRatio: state.pixelRatio });
      sprites.update(camera);
      if (glow) glow.update(camera);
    },
    dispose() {
      sphere.dispose();
      sprites.dispose();
      if (glow) glow.dispose();
    },
  };
}

// The controller. `options`: `skyUrl` (the sky route), `view` (a system view with `attachSky(layer)` and `detachSky()`), `sceneEpochUnix`, `storage`,
// `phone`, `pixelRatio`, `starUrl` (a star's page with `{id}` for its printed ID), `onState(state, note)`, `fetch` (default the global one), `setTimeout`/`clearTimeout` (default the globals), `pollMs`.
// Returns `{ enable(), disable(), enabled(), state(), dispose(), reload() }`; `enable` and `reload` resolve when the sky is drawn or has failed.
export function createSystemSky(THREE, options) {
  const o = options;
  const fetcher = o.fetch || ((...args) => fetch(...args));
  const later = o.setTimeout || ((fn, ms) => setTimeout(fn, ms));
  const cancel = o.clearTimeout || ((id) => clearTimeout(id));
  const pollMs = Number.isFinite(o.pollMs) ? o.pollMs : POLL_MS;
  const limit = skyLimit(!!o.phone);
  let on = readPreference(o.storage);
  let current = { kind: on ? "loading" : "off" };
  let layer = null;
  let timer = null;
  let polls = 0;
  let serial = 0;

  function set(next) {
    current = next;
    if (o.onState) o.onState(current, noteFor(current));
  }

  function drop() {
    if (timer !== null) cancel(timer);
    timer = null;
    if (layer) {
      o.view.detachSky();
      layer.dispose();
      layer = null;
    }
  }

  async function fetchBytes(url) {
    const response = await fetcher(url, { headers: { Accept: "application/octet-stream" } });
    if (!response.ok) {
      const error = new Error("sky " + response.status);
      error.status = response.status;
      throw error;
    }
    return new Uint8Array(await response.arrayBuffer());
  }

  async function load() {
    const mine = ++serial;
    let sky;
    try {
      sky = unpackSky(await fetchBytes(skyUrlFor(o.skyUrl, limit)));
    } catch (err) {
      if (mine === serial && on) {
        drop();
        set({ kind: "unavailable", status: err.status || 0 });
      }
      return current;
    }
    let band = null;
    const wanted = sky.meta && sky.meta.band;
    if (wanted && wanted.href) {
      try {
        band = imageFromHalfFloats(THREE, await fetchBytes(bandUrlFor(wanted)), wanted.width, wanted.height);
      } catch (err) {
        band = null;           // the stars and sprites still draw; only the glow is missing
      }
    }
    if (mine !== serial || !on) return current;
    const next = buildLayer(THREE, sky, band, {
      phone: o.phone, limit, pixelRatio: o.pixelRatio, sceneEpochUnix: o.sceneEpochUnix, starUrl: o.starUrl,
    });
    drop();
    layer = next;
    o.view.attachSky(layer);
    const report = (sky.meta && sky.meta.report) || {};
    set({ kind: report.filling ? "filling" : "ready", report, count: layer.count, band: !!band });
    if (report.filling && polls < MAX_POLLS) {
      timer = later(() => { timer = null; polls += 1; load(); }, pollMs);
    }
    return current;
  }

  return {
    enabled: () => on,
    state: () => current,
    async enable() {
      on = true;
      writePreference(o.storage, true);
      polls = 0;
      set({ kind: "loading" });
      return load();
    },
    disable() {
      on = false;
      serial += 1;
      writePreference(o.storage, false);
      drop();
      set({ kind: "off" });
    },
    reload: () => load(),
    dispose() {
      on = false;
      serial += 1;
      drop();
    },
  };
}
