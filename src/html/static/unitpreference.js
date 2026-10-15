// html/static/unitpreference.js
//
// The reader's unit preference (UX.78): Automatic (the unit ladders as they
// are), Metric only or Customary. Kept in localStorage, which can be missing
// or throw (private windows, blocked site data), so every access is guarded
// and the page works with Automatic without it.

export const SYSTEMS = ["auto", "metric", "customary"];
export const LABELS = { auto: "Automatic", metric: "Metric only", customary: "Customary" };

const KEY = "planetgen.units";

// The saved preference, "auto" when there is none or storage is unavailable.
export function unitSystem() {
  try {
    const value = globalThis.localStorage.getItem(KEY);
    return SYSTEMS.includes(value) ? value : "auto";
  } catch (e) {
    return "auto";
  }
}

// Saves `system` (not remembered, but still returned, when storage is unavailable).
export function setUnitSystem(system) {
  const chosen = SYSTEMS.includes(system) ? system : "auto";
  try {
    if (chosen === "auto") {
      globalThis.localStorage.removeItem(KEY);
    } else {
      globalThis.localStorage.setItem(KEY, chosen);
    }
  } catch (e) {
    // Applied to this page only.
  }
  return chosen;
}
