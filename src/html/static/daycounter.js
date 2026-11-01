// html/static/daycounter.js
//
// The day-and-night counter of an object page (UX.243): a clock that counts the object's own day at its real rate, with the day number, the time
// into the day and whether it is day or night and for how long. It is apart from spin.js, which shows every body turning far faster than it does
// (stand-in SP1); the page says so beside the clock.
//
// The page puts what is stored in data attributes on `[data-day-counter]`: `data-day-hours` (the solar day), `data-night-hours` (the night in it
// at the equator, absent when unknown) and `data-epoch-ms` (where every object's days are counted from, so two pages agree on the day number).
// The day starts at sunrise, daylight runs for the day less the night, then the night.
//
// No DOM at import; `tests/js/daycounter.test.mjs` holds the numbers.

export const MS_PER_HOUR = 3600000;

// Where the clock stands `nowMs` milliseconds since 1970, for a day of `dayHours` counted from `epochMs` and a night of `nightHours` (null when not known):
// `{day, hoursIn, isDay, hoursLeft}` -- the solar day number (0 is the first), the hours since sunrise, whether it is daylight (null when the night is
// not known), and the hours left of the daylight or the night now going.
export function clockAt(nowMs, epochMs, dayHours, nightHours) {
  if (!(dayHours > 0) || !Number.isFinite(dayHours)) return null;
  const elapsed = (nowMs - epochMs) / MS_PER_HOUR;
  const day = Math.floor(elapsed / dayHours);
  const hoursIn = elapsed - day * dayHours;
  if (nightHours === null || nightHours === undefined || !(nightHours >= 0) || nightHours >= dayHours) {
    return { day, hoursIn, isDay: nightHours >= dayHours ? false : null, hoursLeft: null };
  }
  const daylight = dayHours - nightHours;
  const isDay = hoursIn < daylight;
  return { day, hoursIn, isDay, hoursLeft: isDay ? daylight - hoursIn : dayHours - hoursIn };
}

// A span of hours as "7 h 12 min 33 s" (hours only when the span is whole hours; seconds always shown for a span under a day, so the clock moves).
export function spanText(hours) {
  if (!Number.isFinite(hours) || hours < 0) return "–";
  const total = Math.floor(hours * 3600);
  const h = Math.floor(total / 3600);
  const m = Math.floor((total % 3600) / 60);
  const s = total % 60;
  return `${h.toLocaleString("en-US")} h ${String(m).padStart(2, "0")} min ${String(s).padStart(2, "0")} s`;
}

// The sentences of the clock for one reading, or null for a body with no usable day.
export function sentences(reading, dayHours, nightHours) {
  if (!reading) return null;
  const parts = {
    day: `Day ${(reading.day + 1).toLocaleString("en-US")} since the count began`,
    inDay: `${spanText(reading.hoursIn)} into a day of ${spanText(dayHours)}`,
  };
  if (nightHours === 0) {
    parts.phase = "No night: it is daylight all of the day.";
  } else if (reading.isDay === null) {
    parts.phase = "";
  } else if (reading.hoursLeft === null) {
    parts.phase = "It is night all of the day.";
  } else {
    parts.phase = `It is ${reading.isDay ? "day" : "night"}: ${spanText(reading.hoursLeft)} of ${reading.isDay ? "daylight" : "night"} left.`;
  }
  return parts;
}

function number(text) {
  if (text === undefined || text === null || text === "") return null;
  const value = Number(text);
  return Number.isFinite(value) ? value : null;
}

// Starts the clock in `element` (a `[data-day-counter]`), updating each second; returns the stop function.
export function start(element, now = () => Date.now(), setTimer = setInterval, clearTimer = clearInterval) {
  const dayHours = number(element.dataset.dayHours);
  const nightHours = number(element.dataset.nightHours);
  const epoch = number(element.dataset.epochMs);
  const target = (name) => element.querySelector(`[data-day-counter-${name}]`);
  const show = () => {
    const parts = sentences(clockAt(now(), epoch, dayHours, nightHours), dayHours, nightHours);
    if (!parts) return;
    for (const [name, text] of Object.entries(parts)) {
      const node = target(name === "inDay" ? "in-day" : name);
      if (node) node.textContent = text;
    }
  };
  show();
  const timer = setTimer(show, 1000);
  return () => clearTimer(timer);
}

export function startAll(root = document) {
  for (const element of root.querySelectorAll("[data-day-counter]")) {
    if (element.dataset.dayHours) start(element);
  }
}

if (typeof document !== "undefined") startAll();
