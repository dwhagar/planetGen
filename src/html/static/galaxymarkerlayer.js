// html/static/galaxymarkerlayer.js
//
// The Galaxy Map's marker layer, the drawing half (MAP.219; the choosing is galaxymarkers.js): a pool of screen-size
// sprites, one per marker the map is showing, each `MARKER_PX` pixels across whatever its distance, drawn over the
// stars (no depth test), with the glyph of its kind in the kind's color and a dark outline so it reads on any
// background. `set` takes what `selectMarkers` chose, `update` keeps the sizes right as the canvas changes (call it
// every frame) and `entries` lists what is drawn for the picker (mappick.js, a `points` layer).

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { GLYPHS, MARKER_PX } = await import(`./galaxymarkers.js${VERSION_QUERY}`);

const TEXTURE_PX = 64;

// A glyph's texture, drawn once per shape and color.
function makeGlyphTexture(glyph, color, outline) {
  const canvas = document.createElement("canvas");
  canvas.width = canvas.height = TEXTURE_PX;
  GLYPHS[glyph](canvas.getContext("2d"), TEXTURE_PX, color, outline);
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
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

  function textureOf(kind) {
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
      sprite.material.map = textureOf(entry.kind);
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
    const heightPx = canvasEl.clientHeight || 1;
    const scale = (MARKER_PX / heightPx) * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    for (let index = 0; index < drawn.length; index++) sprites[index].scale.set(scale, scale, 1);
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
