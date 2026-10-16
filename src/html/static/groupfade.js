// html/static/groupfade.js
//
// Draws a whole scene group at a share of its strength (MAP.152), for a sector's or a system's scene
// fading in as the camera comes up to it: every material's own strength is scaled from the value it was
// built with (a shader's `fade`, `maxAlpha` or `glowStrength` uniform; a plain material's `opacity`), and at 0
// the group is not drawn at all. Calling it again with another share rescales from the same built values,
// so the shares do not compound.

const UNIFORMS = ["fade", "maxAlpha", "glowStrength"];

// Draws `group` at `share` (0 to 1, held to it); returns the share used.
export function fadeGroup(group, share) {
  const faded = Math.min(1, Math.max(0, share));
  group.visible = faded > 0;
  group.traverse(function (object) {
    const material = object.material;
    if (!material) return;
    const base = material.userData.fadeBase || (material.userData.fadeBase = {});
    if (material.isShaderMaterial) {
      UNIFORMS.forEach(function (name) {
        const uniform = material.uniforms[name];
        if (!uniform) return;
        if (base[name] === undefined) base[name] = uniform.value;
        uniform.value = base[name] * faded;
      });
      return;
    }
    if (base.opacity === undefined) {
      base.opacity = material.opacity;
      base.transparent = material.transparent;
    }
    material.opacity = base.opacity * faded;
    material.transparent = base.transparent || faded < 1;
  });
  return faded;
}
