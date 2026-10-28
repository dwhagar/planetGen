// html/static/galaxymarkerlayer.js
//
// The Galaxy Map's marker layer, the drawing half (MAP.219; the choosing is galaxymarkers.js): a pool of screen-size
// sprites, one per marker the map is showing, each `MARKER_PX` pixels across whatever its distance, drawn over the
// stars (no depth test), with the glyph of its kind in the kind's color and a dark outline so it reads on any
// background. `set` takes what `selectMarkers` chose, `update` keeps the sizes right as the canvas changes (call it
// every frame) and `entries` lists what is drawn for the picker (mappick.js, a `points` layer).

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { GLYPHS, MARKER_PX, tokenOfItem } = await import(`./galaxymarkers.js${VERSION_QUERY}`);
const { shapesOf } = await import(`./poiicons.js${VERSION_QUERY}`);

const TEXTURE_PX = 64;

// A crowd's mark (MAP.305) is drawn this much larger than a lone one, with its count in a badge at the corner.
const CROWD_SCALE = 1.4;
const COUNT_CAP = 99;

// A glyph's texture, drawn once per shape and color.
function makeGlyphTexture(glyph, color, outline) {
  const canvas = document.createElement("canvas");
  canvas.width = canvas.height = TEXTURE_PX;
  GLYPHS[glyph](canvas.getContext("2d"), TEXTURE_PX, color, outline);
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

// A point-of-interest symbol's shapes (MAP.320), read once from its file under static/poi/ (`shapesOf`): a promise.
const shapeFiles = new Map();
function shapesOfIcon(id) {
  if (!shapeFiles.has(id)) {
    const url = new URL(`poi/${id}.svg${VERSION_QUERY}`, import.meta.url);
    shapeFiles.set(id, fetch(url).then(function (response) {
      if (!response.ok) throw new Error("no icon file " + url.pathname);
      return response.text();
    }).then(function (text) {
      const shapes = shapesOf(text);
      shapes.paths = shapes.paths.map(function (path) { return { path: new Path2D(path.d), evenodd: path.evenodd }; });
      return shapes;
    }));
  }
  return shapeFiles.get(id);
}

// The texture of a symbol in `color` with the dark or light edge `outline` (the paint order of the pages' symbols, edge
// under the fill). It is blank until the file has been read; the map draws every frame, so it shows from the next one.
function makeIconTexture(icon, color, outline, count) {
  const canvas = document.createElement("canvas");
  canvas.width = canvas.height = TEXTURE_PX;
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  shapesOfIcon(icon).then(function (shapes) {
    const ctx = canvas.getContext("2d");
    ctx.scale(TEXTURE_PX / shapes.size, TEXTURE_PX / shapes.size);
    ctx.lineJoin = "round";
    ctx.lineWidth = 1.5;
    ctx.strokeStyle = outline;
    ctx.fillStyle = color;
    shapes.paths.forEach(function (shape) {
      ctx.stroke(shape.path);
      ctx.fill(shape.path, shape.evenodd ? "evenodd" : "nonzero");
    });
    if (count > 1) drawCount(ctx, shapes.size, count, color, outline);
    texture.needsUpdate = true;
  }).catch(function (error) { console.error(error); });
  return texture;
}

// The count of a crowd in a round badge at the symbol's top right corner, in a canvas whose drawing is `size` units square.
function drawCount(ctx, size, count, color, outline) {
  const r = size * 0.27;
  const x = size - r;
  const y = r;
  ctx.beginPath();
  ctx.arc(x, y, r, 0, Math.PI * 2);
  ctx.fillStyle = outline;
  ctx.fill();
  ctx.lineWidth = size * 0.04;
  ctx.strokeStyle = color;
  ctx.stroke();
  ctx.fillStyle = color;
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  const text = count > COUNT_CAP ? COUNT_CAP + "+" : String(count);
  ctx.font = "bold " + Math.round(r * (text.length > 2 ? 0.95 : 1.25)) + "px sans-serif";
  ctx.fillText(text, x, y + r * 0.06);
}

// The color a CSS custom property of the page has now (style.css's tokens; they differ by theme).
function tokenColor(token) {
  return getComputedStyle(document.documentElement).getPropertyValue(token).trim() || "#888888";
}

// A marker layer in `scene`, for `camera` on `canvasEl`. `outline` is the glyph's outline color (CSS). `options`:
// `listIn` (an element) gets a visually hidden list with a button per marker drawn, for the keyboard and screen readers
// (a canvas has no focusable parts, as the Sector Map's list); `onChoose(entry)` is called when one is activated.
export function createMarkerLayer(scene, camera, canvasEl, outline, options) {
  const o = options || {};
  let list = null;
  if (o.listIn) {
    list = document.createElement("ul");
    list.className = "starmap-sr-list starmap-sr-markers sr-only";
    list.setAttribute("aria-label", "Markers on the map");
    o.listIn.appendChild(list);
  }
  const group = new THREE.Group();
  group.renderOrder = 4;
  scene.add(group);
  const textures = new Map();
  const sprites = [];
  let drawn = [];

  function textureOf(kind, item) {
    if (kind.icon) {
      const color = tokenColor(tokenOfItem(kind, item));
      const count = item && item.count > 1 ? item.count : 0;
      const iconKey = kind.icon + "|" + color + "|" + outline + "|" + count;
      if (!textures.has(iconKey)) textures.set(iconKey, makeIconTexture(kind.icon, color, outline, count));
      return textures.get(iconKey);
    }
    const key = kind.glyph + "|" + kind.color + "|" + outline;
    if (!textures.has(key)) textures.set(key, makeGlyphTexture(kind.glyph, kind.color, outline));
    return textures.get(key);
  }

  function spriteAt(index) {
    if (!sprites[index]) {
      const sprite = new THREE.Sprite(new THREE.SpriteMaterial({
        transparent: true, depthWrite: false, depthTest: false, sizeAttenuation: false,
      }));
      sprite.renderOrder = 4;
      group.add(sprite);
      sprites[index] = sprite;
    }
    return sprites[index];
  }

  // Shows `chosen` (`[{kind, item}]`); the sprites left over are hidden.
  function set(chosen) {
    drawn = chosen.map(function (entry) {
      return { x: entry.item.x, y: entry.item.y, z: entry.item.z, kind: entry.kind, item: entry.item };
    });
    drawn.forEach(function (entry, index) {
      const sprite = spriteAt(index);
      sprite.material.map = textureOf(entry.kind, entry.item);
      sprite.material.needsUpdate = true;
      sprite.position.set(entry.x, entry.y, entry.z);
      sprite.visible = true;
    });
    for (let index = drawn.length; index < sprites.length; index++) sprites[index].visible = false;
    if (list) {
      list.textContent = "";
      drawn.forEach(function (entry) {
        const item = document.createElement("li");
        const button = document.createElement("button");
        button.type = "button";
        button.textContent = entry.kind.describe(entry.item);
        button.addEventListener("click", function () { if (o.onChoose) o.onChoose(entry); });
        item.appendChild(button);
        list.appendChild(item);
      });
    }
    update();
  }

  // A sprite without size attenuation is scaled in units of the view's height at distance 1 (mappick.createRing).
  function update() {
    if (!drawn.length) return;
    // A marker on a body that moves (MAP.321) names the entry that holds where the body is now (`item.follow`: x, y, z).
    for (let index = 0; index < drawn.length; index++) {
      const follow = drawn[index].item.follow;
      if (!follow) continue;
      drawn[index].x = follow.x;
      drawn[index].y = follow.y;
      drawn[index].z = follow.z;
      sprites[index].position.set(follow.x, follow.y, follow.z);
    }
    const heightPx = canvasEl.clientHeight || 1;
    const scale = (MARKER_PX / heightPx) * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    for (let index = 0; index < drawn.length; index++) {
      const size = drawn[index].item.merged ? scale * CROWD_SCALE : scale;
      sprites[index].scale.set(size, size, 1);
    }
  }

  return {
    set: set,
    update: update,
    entries: function () { return drawn; },
    group: group,
    // Hides or shows the whole layer (a view with no sectors to be seen through).
    setVisible: function (on) { group.visible = !!on; },
  };
}
