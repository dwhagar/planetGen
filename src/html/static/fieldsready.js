// html/static/fieldsready.js
//
// The form fields are Shoelace elements (UX.125, templates/partials/ui_fields.html), defined by components.js a moment
// after the page's own scripts first run. Until then an `sl-input` has no `value`, no `checked` and no events, so a
// script that reads or listens to the fields waits for `fieldsReady()` first. It resolves when the three field elements
// are defined, or after WAIT_MS if they never are (the page then shows the plain twins in its <noscript> and the
// script finds no field to drive), or at once where there is no custom element registry.

export const WAIT_MS = 3000;
export const FIELD_ELEMENTS = ["sl-input", "sl-select", "sl-checkbox"];

export function fieldsReady() {
  if (typeof customElements === "undefined") return Promise.resolve();
  const defined = Promise.all(FIELD_ELEMENTS.map((name) => customElements.whenDefined(name)));
  const gaveUp = new Promise((resolve) => setTimeout(resolve, WAIT_MS));
  return Promise.race([defined, gaveUp]);
}
