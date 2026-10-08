// html/static/orbitclock.js
//
// The time control of the 3D system view (MAP.70): a clock that gives the
// years past the scene's epoch for the positions of orbitpositions.js. It
// starts at "now" (the real time since the epoch), plays at a rate (real
// time, then faster), pauses, and jumps back to now. Plain
// state, no DOM: the view calls tick(nowMs) each frame and draws
// positionsAt(scene, clock.years()).

export const SECONDS_PER_YEAR = 365.25 * 86400;

// Years of scene time per real second at each speed step; the first is real time.
export const RATES = [
  { label: "1×", yearsPerSecond: 1 / SECONDS_PER_YEAR },
  { label: "1 day/s", yearsPerSecond: 86400 / SECONDS_PER_YEAR },
  { label: "1 month/s", yearsPerSecond: 30.44 * 86400 / SECONDS_PER_YEAR },
  { label: "1 year/s", yearsPerSecond: 1 },
  { label: "10 years/s", yearsPerSecond: 10 },
  { label: "100 years/s", yearsPerSecond: 100 },
];

// options: epochUnix (seconds; the scene's epoch, or null), nowMs (Date.now()).
// The clock starts on the real present, running in real time.
export function createOrbitClock(options) {
  const epochMs = options.epochUnix == null ? options.nowMs : options.epochUnix * 1000;
  const realYears = (ms) => (ms - epochMs) / 1000 / SECONDS_PER_YEAR;
  let sceneYears = realYears(options.nowMs);
  let playing = true;
  let rateIndex = 0;
  let lastMs = options.nowMs;

  const api = {
    // Scene time as years after the epoch.
    years: function () { return sceneYears; },
    playing: function () { return playing; },
    rate: function () { return RATES[rateIndex]; },
    rateIndex: function () { return rateIndex; },
    // Advances to wall-clock `nowMs`; only a playing clock moves, at its rate.
    tick: function (nowMs) {
      const dtSeconds = Math.max(0, (nowMs - lastMs) / 1000);
      lastMs = nowMs;
      if (playing) sceneYears += dtSeconds * RATES[rateIndex].yearsPerSecond;
      return sceneYears;
    },
    play: function () { playing = true; },
    pause: function () { playing = false; },
    toggle: function () { playing = !playing; },
    faster: function () { rateIndex = Math.min(RATES.length - 1, rateIndex + 1); playing = true; },
    slower: function () { rateIndex = Math.max(0, rateIndex - 1); },
    // Back to the real present, playing in real time.
    now: function (nowMs) {
      sceneYears = realYears(nowMs);
      rateIndex = 0;
      playing = true;
      lastMs = nowMs;
    },
  };
  return api;
}
