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
  // The clear margin in pixels the label engine keeps between one name and the next (MAP.299).
  "map.labelBuffer": { fallback: 4, min: 0, max: 40 },
};

// Each setting with a fixed set of words: its key, default and the words a stored value must be one of.
export const CHOICE_SETTINGS = {
  // MAP.336: where Go to a region lands: "frame" (the region's charted part fills the view, the default; Boss, 6:19 PM EDT 2026-10-28) or
  // "densest" (the cluster of its most densely filled sectors).
  "map.regionLanding": { fallback: "frame", choices: ["frame", "densest"] },
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

// The label buffer setting: the clear margin in pixels between two names.
export function labelBufferPx(storage) {
  return readNumber(storage, "map.labelBuffer");
}

// The setting `key` as one of its words: the stored value when it is one, else the default.
export function readChoice(storage, key) {
  const setting = CHOICE_SETTINGS[key];
  if (!setting) throw new Error("unknown choice setting: " + key);
  try {
    const text = storage && storage.getItem(PREFIX + key);
    const word = text === null || text === undefined ? "" : String(text).trim();
    return setting.choices.includes(word) ? word : setting.fallback;
  } catch (error) {
    return setting.fallback;
  }
}

// MAP.336: where Go to a region lands, "frame" or "densest".
export function regionLanding(storage) {
  return readChoice(storage, "map.regionLanding");
}

// --- Writing the settings, and the words the preferences page shows (UX.236) ----------------------------------------------------------

// What the preferences page says about each setting: its name, a plain-words description, and for a number the unit its value is in.
// Every key of NUMBER_SETTINGS and CHOICE_SETTINGS has an entry (tests/js/mappreferences.test.mjs holds that), listed in this order.
export const PREFERENCE_TEXT = {
  "map.regionLanding": {
    group: "Galaxy Map",
    label: "Where Go to a region lands",
    description: "When you pick a region on the Galaxy Map, the camera flies there and stops.",
    choices: {
      frame: "Frame the region: the charted part of the region fills the view.",
      densest: "Land on the densest part: the cluster of the region's most heavily filled sectors.",
    },
  },
  "map.blackHoleDensity": {
    group: "Galaxy Map",
    label: "Black hole marks",
    description: "How crowded the Galaxy Map's black hole marks may get: the side of the square of screen that one mark may stand on. "
      + "A smaller number shows more marks.",
    unit: "pixels",
  },
  "map.labelBuffer": {
    group: "Galaxy Map",
    label: "Space between names",
    description: "The clear margin the Galaxy Map keeps between one name and the next. A bigger number shows fewer names.",
    unit: "pixels",
  },
};

// Every setting in the order the page lists them: {key, kind: "choice" | "number", fallback, ...its range or words, ...its text}.
export function preferenceList() {
  const entries = [];
  Object.keys(PREFERENCE_TEXT).forEach(function (key) {
    const text = PREFERENCE_TEXT[key];
    if (CHOICE_SETTINGS[key]) entries.push(Object.assign({ key: key, kind: "choice" }, CHOICE_SETTINGS[key], text));
    else if (NUMBER_SETTINGS[key]) entries.push(Object.assign({ key: key, kind: "number" }, NUMBER_SETTINGS[key], text));
  });
  return entries;
}

// The setting `key` as it stands: the word or number the map would read (the default when nothing valid is stored).
export function readSetting(storage, key) {
  return CHOICE_SETTINGS[key] ? readChoice(storage, key) : readNumber(storage, key);
}

// Stores `value` for the setting `key`; false (and nothing stored) for a value the reader would not accept, or a store that refuses.
export function writeSetting(storage, key, value) {
  const choice = CHOICE_SETTINGS[key];
  const number = NUMBER_SETTINGS[key];
  if (!choice && !number) throw new Error("unknown setting: " + key);
  let text = String(value).trim();
  if (choice) {
    if (!choice.choices.includes(text)) return false;
  } else {
    const parsed = Number(text);
    if (text === "" || !Number.isFinite(parsed) || parsed < number.min || parsed > number.max) return false;
    text = String(parsed);
  }
  try {
    storage.setItem(PREFIX + key, text);
    return true;
  } catch (error) {
    return false;
  }
}

// Forgets the stored value of one setting, so it reads as its default; false when the store refuses.
export function clearSetting(storage, key) {
  try {
    storage.removeItem(PREFIX + key);
    return true;
  } catch (error) {
    return false;
  }
}

// Forgets every stored setting, so each reads as its default; false when the store refuses.
export function resetSettings(storage) {
  return Object.keys(PREFERENCE_TEXT).every(function (key) { return clearSetting(storage, key); });
}
