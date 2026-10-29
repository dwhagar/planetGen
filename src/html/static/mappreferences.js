// html/static/mappreferences.js
//
// MAP.339: the per-browser settings the Galaxy Map reads (UX.236 plans the page that sets them and the account storage they move to;
// this is the reader the map calls, so a key that is missing, unreadable or out of range reads as its default). A value is kept in
// the browser's `localStorage` under `planetgen.pref.` plus the key, as the number's text. `storage` is `localStorage` or anything with
// `getItem`, which can be missing or throw (private windows), so every read is guarded. No DOM and no three.js:
// tests/js/mappreferences.test.mjs holds it.

export const PREFIX = "planetgen.pref.";

// Each number setting: its key, default and the range a stored value must lie in.
export const NUMBER_SETTINGS = {
  // The side in pixels of the square one black hole mark may stand on (a smaller number shows more marks); one mark per 120 by
  // 120 pixels by default (Boss, 6:22 PM EDT 2026-10-28).
  "map.blackHoleDensity": { fallback: 120, min: 30, max: 1000 },
};

// The setting `key` as a number: the stored value when it is a finite number inside its range, else the default.
export function readNumber(storage, key) {
  const setting = NUMBER_SETTINGS[key];
  if (!setting) throw new Error("unknown number setting: " + key);
  try {
    const text = storage && storage.getItem(PREFIX + key);
    if (text === null || text === undefined || String(text).trim() === "") return setting.fallback;
    const value = Number(text);
    return Number.isFinite(value) && value >= setting.min && value <= setting.max ? value : setting.fallback;
  } catch (error) {
    return setting.fallback;
  }
}

// The black hole density setting: the side in pixels of the square one mark may stand on.
export function blackHoleDensityPx(storage) {
  return readNumber(storage, "map.blackHoleDensity");
}
