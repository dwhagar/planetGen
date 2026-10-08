// html/static/navpick.js
//
// Picking a NAV course on the maps (NAV.29, NAV.33). The Galaxy Map and the
// Sector Map hold one pick state: which end is being chosen ("from", the
// start, or "to", the destination) and the end already chosen. The state
// starts from the NAV page's "Pick on Galaxy Map" (`?pick=to&from=...`) or
// from nothing, and a click on "Start Here" or "End Here" in an object's
// info panel moves it on without leaving the map: the first end is kept,
// the view stays where it is, and the user zooms out, in and across to find
// the other end the same way. Only the second end leaves for the NAV page,
// with both ends in its URL.
//
// The page-level URLs (`?pick=...`) are the state's `query()`; the panel's
// buttons are `actionsFor(endpoint, name)`. Plain data in, plain data out:
// no DOM here but the banner's.

const START_LABEL = "Start Here";
const END_LABEL = "End Here";

function labelOf(pick) {
  return pick === "from" ? START_LABEL : END_LABEL;
}

function encodeEndpoint(endpoint) {
  return encodeURIComponent(endpoint).replace(/%3A/gi, ":");
}

// options: navUrl (the NAV page's URL), pick ("from", "to" or empty), other
// (the end already chosen, an endpoint string such as "system:12", or empty),
// otherName (its name, for the banner), cancelUrl (where Cancel goes for a
// pick the NAV page began), onChange(state) (after every change).
export function createNavPick(options) {
  const navUrl = String(options.navUrl || "/nav");
  let pick = options.pick === "from" || options.pick === "to" ? options.pick : null;
  let other = pick && options.other ? String(options.other) : null;
  let otherName = other ? options.otherName || other : null;
  const cancelUrl = options.cancelUrl || null;
  let began = false; // the pick was begun on the map, so Cancel clears it in place

  function changed() {
    if (options.onChange) options.onChange(api);
  }

  // The NAV page's URL with both ends.
  function courseUrl(endpoint) {
    const start = pick === "from" ? endpoint : other;
    const end = pick === "from" ? other : endpoint;
    const join = navUrl.indexOf("?") < 0 ? "?" : "&";
    return navUrl + join + "from=" + encodeEndpoint(start) + "&to=" + encodeEndpoint(end);
  }

  // Begins picking the other end with `endpoint` as this one: "Start Here"
  // chooses the destination next, "End Here" the start.
  function begin(side, endpoint, name) {
    began = true;
    pick = side === "from" ? "to" : "from";
    other = endpoint;
    otherName = name || endpoint;
    changed();
  }

  function clear() {
    began = false;
    pick = null;
    other = null;
    otherName = null;
    changed();
  }

  const api = {
    active: function () { return pick !== null; },
    pick: function () { return pick; },
    other: function () { return other; },
    otherName: function () { return otherName; },
    // The other end's field and value, as the bookmarks menu keeps them.
    keepName: function () { return pick === "from" ? "to" : "from"; },
    // The query every URL of the map carries while a pick is on.
    query: function () {
      if (!pick) return "";
      return "?pick=" + pick + (other ? "&" + api.keepName() + "=" + encodeEndpoint(other) : "");
    },
    // The panel's NAV buttons for an object: with no pick, "Start Here" and
    // "End Here"; while picking, only the end being chosen (NAV.30: nothing
    // that would leave the course). Each is {label, href} or {label, onClick}.
    actionsFor: function (endpoint, name) {
      if (!endpoint) return [];
      if (!pick) {
        return [
          { label: START_LABEL, onClick: function () { begin("from", endpoint, name); } },
          { label: END_LABEL, onClick: function () { begin("to", endpoint, name); } },
        ];
      }
      if (other) return [{ label: labelOf(pick), href: courseUrl(endpoint), primary: true }];
      return [{ label: labelOf(pick), primary: true, onClick: function () { begin(pick, endpoint, name); } }];
    },
    // The banner's words: [lead, rest], e.g. ["Start: Sol. Choosing a destination:", " press ..."].
    bannerParts: function (what) {
      if (!pick) return null;
      const choosing = pick === "from" ? "Choosing a start" : "Choosing a destination";
      const lead = other && began
        ? (pick === "from" ? "Destination: " : "Start: ") + otherName + ". " + choosing
        : choosing;
      return [lead, " press “" + labelOf(pick) + "” to take this " + (what || "object") + ", or use a bookmark. "];
    },
    // Where Cancel goes: the NAV page's own for a pick it began; none for
    // one begun here (Cancel clears that in place).
    cancelUrl: function () { return began ? null : cancelUrl; },
    clear: clear,
    begin: begin,
  };
  return api;
}

export const NAV_LABELS = { start: START_LABEL, end: END_LABEL };
