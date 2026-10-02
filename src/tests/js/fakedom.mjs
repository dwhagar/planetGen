// tests/js/fakedom.mjs
//
// A small stand-in for the browser's DOM, enough to run the site's page
// scripts (html/static/*.js) under node without a browser: elements with
// attributes, `dataset`, `classList`, text, children, simple CSS
// selectors (tag, #id, .class, [attr], [attr="value"], descendant and
// comma lists), events that capture and bubble, a `<select>`'s options,
// and a window with fake timers, location, history (pushState, back,
// forward and popstate), localStorage and fetch.
//
// installDom() puts `window`, `document`, `location`, `history` and the
// rest on globalThis, so a page module imported after it sees them as
// globals, just as in a page. It is not a browser: layout is whatever a
// test sets (el.rect, el.clientWidth), and nothing is drawn.

import { readFileSync } from "node:fs";

export class FakeEvent {
  constructor(type, init) {
    Object.assign(this, init || {});
    this.type = type;
    this.bubbles = !init || init.bubbles !== false;
    this.defaultPrevented = false;
    this.propagationStopped = false;
  }

  preventDefault() {
    this.defaultPrevented = true;
  }

  stopPropagation() {
    this.propagationStopped = true;
  }

  stopImmediatePropagation() {
    this.propagationStopped = true;
  }
}

class EventTargetish {
  constructor() {
    this.listeners = new Map();
  }

  addEventListener(type, fn, options) {
    if (!fn) return;
    const capture = options === true || !!(options && options.capture);
    const once = !!(options && options.once);
    if (!this.listeners.has(type)) this.listeners.set(type, []);
    this.listeners.get(type).push({ fn, capture, once });
  }

  removeEventListener(type, fn, options) {
    const capture = options === true || !!(options && options.capture);
    const list = this.listeners.get(type) || [];
    this.listeners.set(type, list.filter((l) => l.fn !== fn || l.capture !== capture));
  }

  fire(event, capturePhase) {
    const list = (this.listeners.get(event.type) || []).slice();
    for (const listener of list) {
      if (capturePhase !== null && listener.capture !== capturePhase) continue;
      if (listener.once) this.removeEventListener(event.type, listener.fn, { capture: listener.capture });
      event.currentTarget = this;
      listener.fn.call(this, event);
      if (event.propagationStopped) return;
    }
  }

  // The path from the outermost ancestor down to this target.
  eventPath() {
    return [this];
  }

  dispatchEvent(event) {
    if (!event.target) event.target = this;
    const path = this.eventPath();
    for (let n = 0; n < path.length - 1 && !event.propagationStopped; n++) path[n].fire(event, true);
    if (!event.propagationStopped) this.fire(event, null);
    for (let n = path.length - 2; n >= 0 && event.bubbles && !event.propagationStopped; n--) path[n].fire(event, false);
    return !event.defaultPrevented;
  }
}

class TextNode {
  constructor(data) {
    this.nodeType = 3;
    this.data = String(data);
    this.parentNode = null;
  }

  get textContent() {
    return this.data;
  }

  remove() {
    if (this.parentNode) this.parentNode.removeChild(this);
  }
}

function camelToData(name) {
  return "data-" + name.replace(/[A-Z]/g, (c) => "-" + c.toLowerCase());
}

function dataToCamel(name) {
  return name.slice(5).replace(/-([a-z])/g, (_, c) => c.toUpperCase());
}

// --- Selectors -------------------------------------------------------------

function splitOutside(text, separator) {
  const parts = [];
  let depth = 0;
  let quote = null;
  let current = "";
  for (const c of text) {
    if (quote) {
      if (c === quote) quote = null;
    } else if (c === '"' || c === "'") {
      quote = c;
    } else if (c === "[" || c === "(") {
      depth++;
    } else if (c === "]" || c === ")") {
      depth--;
    } else if (depth === 0 && (separator === " " ? /\s/.test(c) : c === separator)) {
      if (current.trim()) parts.push(current.trim());
      current = "";
      continue;
    }
    current += c;
  }
  if (current.trim()) parts.push(current.trim());
  return parts;
}

function parseCompound(text) {
  const out = { tag: null, id: null, classes: [], attrs: [], not: [] };
  const re = /^([a-zA-Z][\w-]*|\*)|#([\w-]+)|\.([\w-]+)|\[([\w-]+)(?:([~^$*|]?=)(?:"([^"]*)"|'([^']*)'|([^\]]*)))?\]|:not\(([^)]*)\)/g;
  let consumed = 0;
  let match;
  while ((match = re.exec(text)) && match[0]) {
    if (match.index !== consumed) break;
    consumed += match[0].length;
    if (match[1]) out.tag = match[1] === "*" ? null : match[1].toLowerCase();
    else if (match[2]) out.id = match[2];
    else if (match[3]) out.classes.push(match[3]);
    else if (match[4]) out.attrs.push({ name: match[4], op: match[5] || null, value: match[6] ?? match[7] ?? match[8] ?? null });
    else if (match[9] !== undefined) out.not.push(parseCompound(match[9].trim()));
  }
  if (consumed !== text.length) throw new Error("fakedom: unsupported selector " + JSON.stringify(text));
  return out;
}

function matchesCompound(el, c) {
  if (c.tag && el.localName !== c.tag) return false;
  if (c.id && el.getAttribute("id") !== c.id) return false;
  for (const name of c.classes) if (!el.classList.contains(name)) return false;
  for (const a of c.attrs) {
    const value = el.getAttribute(a.name);
    if (value === null) return false;
    if (a.op === "=" && value !== a.value) return false;
    if (a.op === "~=" && !value.split(/\s+/).includes(a.value)) return false;
    if (a.op === "^=" && !value.startsWith(a.value)) return false;
  }
  for (const n of c.not) if (matchesCompound(el, n)) return false;
  return true;
}

function matchesSelector(el, selector) {
  return splitOutside(selector, ",").some(function (one) {
    const parts = splitOutside(one, " ").map(parseCompound);
    if (!matchesCompound(el, parts[parts.length - 1])) return false;
    let up = el.parentElement;
    for (let n = parts.length - 2; n >= 0; n--) {
      while (up && !matchesCompound(up, parts[n])) up = up.parentElement;
      if (!up) return false;
      up = up.parentElement;
    }
    return true;
  });
}

// --- Elements --------------------------------------------------------------

const BOOLEAN_PROPS = ["hidden", "disabled", "checked", "required"];
const NUMBER_PROPS = ["width", "height"];
const STRING_PROPS = { id: "id", className: "class", type: "type", href: "href", name: "name", min: "min", max: "max", step: "step", title: "title" };

export class FakeElement extends EventTargetish {
  constructor(tagName, ownerDocument) {
    super();
    this.nodeType = 1;
    this.localName = String(tagName).toLowerCase();
    this.tagName = this.localName.toUpperCase();
    this.ownerDocument = ownerDocument || null;
    this.attributes = new Map();
    this.childNodes = [];
    this.parentNode = null;
    this.style = { setProperty(name, value) { this[name] = value; }, removeProperty(name) { delete this[name]; } };
    this.rect = null;
    this.clientWidth = 0;
    this.clientHeight = 0;
    this.offsetWidth = 0;
    this.offsetHeight = 0;
    this.scrollTop = 0;
    this.scrollHeight = 0;
    this._value = null;
    this._selected = null;
    const el = this;
    this.dataset = new Proxy({}, {
      get(_, key) {
        return typeof key === "string" ? (el.getAttribute(camelToData(key)) ?? undefined) : undefined;
      },
      set(_, key, value) {
        el.setAttribute(camelToData(key), value);
        return true;
      },
      has(_, key) {
        return el.hasAttribute(camelToData(key));
      },
      deleteProperty(_, key) {
        el.removeAttribute(camelToData(key));
        return true;
      },
      ownKeys() {
        return Array.from(el.attributes.keys()).filter((n) => n.startsWith("data-")).map(dataToCamel);
      },
      getOwnPropertyDescriptor(_, key) {
        return el.hasAttribute(camelToData(key)) ? { enumerable: true, configurable: true, value: el.getAttribute(camelToData(key)) } : undefined;
      },
    });
    this.classList = {
      contains: (name) => (el.getAttribute("class") || "").split(/\s+/).includes(name),
      add: (...names) => names.forEach((name) => { if (!el.classList.contains(name)) el.setAttribute("class", ((el.getAttribute("class") || "") + " " + name).trim()); }),
      remove: (...names) => names.forEach((name) => el.setAttribute("class", (el.getAttribute("class") || "").split(/\s+/).filter((c) => c && c !== name).join(" "))),
      toggle: (name, force) => {
        const on = force === undefined ? !el.classList.contains(name) : !!force;
        if (on) el.classList.add(name);
        else el.classList.remove(name);
        return on;
      },
    };
  }

  eventPath() {
    const path = [];
    for (let node = this; node; node = node.parentNode) path.unshift(node);
    const doc = path[0] && path[0].nodeType === 9 ? path[0] : null;
    if (doc && doc.defaultView) path.unshift(doc.defaultView);
    return path;
  }

  getAttribute(name) {
    return this.attributes.has(name) ? this.attributes.get(name) : null;
  }

  setAttribute(name, value) {
    this.attributes.set(name, String(value));
  }

  removeAttribute(name) {
    this.attributes.delete(name);
  }

  hasAttribute(name) {
    return this.attributes.has(name);
  }

  get children() {
    return this.childNodes.filter((n) => n.nodeType === 1);
  }

  get firstElementChild() {
    return this.children[0] || null;
  }

  get parentElement() {
    return this.parentNode && this.parentNode.nodeType === 1 ? this.parentNode : null;
  }

  get isConnected() {
    let node = this;
    while (node.parentNode) node = node.parentNode;
    return node.nodeType === 9;
  }

  appendChild(child) {
    if (typeof child === "string") child = new TextNode(child);
    if (child.parentNode) child.parentNode.removeChild(child);
    child.parentNode = this;
    this.childNodes.push(child);
    return child;
  }

  append(...nodes) {
    nodes.forEach((n) => this.appendChild(n));
  }

  prepend(...nodes) {
    nodes.reverse().forEach((n) => this.insertBefore(typeof n === "string" ? new TextNode(n) : n, this.childNodes[0] || null));
  }

  insertBefore(child, before) {
    if (child.parentNode) child.parentNode.removeChild(child);
    child.parentNode = this;
    const at = before ? this.childNodes.indexOf(before) : -1;
    if (at < 0) this.childNodes.push(child);
    else this.childNodes.splice(at, 0, child);
    return child;
  }

  removeChild(child) {
    const at = this.childNodes.indexOf(child);
    if (at >= 0) this.childNodes.splice(at, 1);
    child.parentNode = null;
    return child;
  }

  remove() {
    if (this.parentNode) this.parentNode.removeChild(this);
  }

  replaceChildren(...nodes) {
    this.childNodes.slice().forEach((n) => this.removeChild(n));
    nodes.forEach((n) => this.appendChild(n));
  }

  get textContent() {
    return this.childNodes.map((n) => n.textContent).join("");
  }

  set textContent(value) {
    this.childNodes.slice().forEach((n) => this.removeChild(n));
    if (value !== "" && value !== null && value !== undefined) this.appendChild(new TextNode(value));
  }

  get innerHTML() {
    return this.textContent;
  }

  set innerHTML(value) {
    if (value) throw new Error("fakedom: innerHTML parsing is not supported");
    this.textContent = "";
  }

  descendants() {
    const out = [];
    const walk = (node) => node.children.forEach((child) => { out.push(child); walk(child); });
    walk(this);
    return out;
  }

  querySelectorAll(selector) {
    return this.descendants().filter((el) => matchesSelector(el, selector));
  }

  querySelector(selector) {
    return this.descendants().find((el) => matchesSelector(el, selector)) || null;
  }

  matches(selector) {
    return matchesSelector(this, selector);
  }

  closest(selector) {
    for (let el = this; el && el.nodeType === 1; el = el.parentNode) if (el.matches(selector)) return el;
    return null;
  }

  getElementById(id) {
    return this.descendants().find((el) => el.getAttribute("id") === id) || null;
  }

  getBoundingClientRect() {
    const r = this.rect || { left: 0, top: 0, width: 0, height: 0 };
    return { left: r.left, top: r.top, width: r.width, height: r.height, x: r.left, y: r.top, right: r.left + r.width, bottom: r.top + r.height };
  }

  setPointerCapture() {}

  scrollIntoView() {}

  // A 2D context that draws nothing (textures made for three.js), and no
  // WebGL.
  getContext(kind) {
    if (kind !== "2d") return null;
    const gradient = { addColorStop() {} };
    return new Proxy({ canvas: this }, {
      get(target, prop) {
        if (prop in target) return target[prop];
        if (prop === "measureText") return (text) => ({ width: 8 * String(text).length });
        if (prop === "createRadialGradient" || prop === "createLinearGradient") return () => gradient;
        if (prop === "getImageData" || prop === "createImageData") return (x, y, w, hgt) => ({ data: new Uint8ClampedArray(4 * (w || 1) * (hgt || 1)), width: w, height: hgt });
        return () => {};
      },
      set(target, prop, value) {
        target[prop] = value;
        return true;
      },
    });
  }

  releasePointerCapture() {}

  contains(node) {
    for (let up = node; up; up = up.parentNode) if (up === this) return true;
    return false;
  }

  focus() {
    if (this.ownerDocument) this.ownerDocument.activeElement = this;
  }

  blur() {
    if (this.ownerDocument && this.ownerDocument.activeElement === this) this.ownerDocument.activeElement = null;
  }

  click() {
    if (this.disabled) return;
    this.dispatchEvent(new FakeEvent("click", { button: 0 }));
  }

  // <input>, <select> and <option>.
  get value() {
    if (this.localName === "progress") return this.hasAttribute("value") ? Number(this.getAttribute("value")) : 0;
    if (this.localName === "select") {
      const chosen = this.selectedOptions[0];
      return chosen ? chosen.value : "";
    }
    if (this._value !== null) return this._value;
    if (this.localName === "option") return this.getAttribute("value") ?? this.textContent;
    return this.getAttribute("value") ?? "";
  }

  set value(v) {
    if (this.localName === "progress") {
      this.setAttribute("value", v);
      return;
    }
    if (this.localName === "select") {
      this.options.forEach((o) => { o._selected = o.value === String(v); });
      return;
    }
    this._value = String(v);
  }

  get options() {
    return this.querySelectorAll("option");
  }

  get selectedOptions() {
    const options = this.options;
    const picked = options.filter((o) => o.selected);
    if (picked.length) return [picked[0]];
    return options.length && this.localName === "select" ? [options.find((o) => !o.disabled) || options[0]] : [];
  }

  get selected() {
    return this._selected === null ? this.hasAttribute("selected") : this._selected;
  }

  set selected(on) {
    const select = this.closest("select");
    if (on && select) select.options.forEach((o) => { o._selected = false; });
    this._selected = !!on;
  }
}

for (const prop of BOOLEAN_PROPS) {
  Object.defineProperty(FakeElement.prototype, prop, {
    get() { return this.hasAttribute(prop); },
    set(on) { if (on) this.setAttribute(prop, ""); else this.removeAttribute(prop); },
  });
}

for (const prop of NUMBER_PROPS) {
  Object.defineProperty(FakeElement.prototype, prop, {
    get() { return Number(this.getAttribute(prop) || 0); },
    set(v) { this.setAttribute(prop, v); },
  });
}

for (const [prop, attr] of Object.entries(STRING_PROPS)) {
  Object.defineProperty(FakeElement.prototype, prop, {
    get() { return this.getAttribute(attr) ?? ""; },
    set(v) { this.setAttribute(attr, v); },
  });
}

export class FakeDocument extends FakeElement {
  constructor() {
    super("#document", null);
    this.nodeType = 9;
    this.ownerDocument = this;
    this.readyState = "complete";
    this.activeElement = null;
    this.defaultView = null;
    this.documentElement = this.createElement("html");
    this.appendChild(this.documentElement);
    this.head = this.createElement("head");
    this.body = this.createElement("body");
    this.documentElement.append(this.head, this.body);
  }

  createElement(tag) {
    return new FakeElement(tag, this);
  }

  createElementNS(_, tag) {
    return new FakeElement(tag, this);
  }

  createTextNode(text) {
    return new TextNode(text);
  }
}

// Builds an element tree from a plain description:
// h("div", {id: "x", "data-job": "7"}, [h("span", {}, "text")]).
export function h(tag, attrs, children) {
  const el = globalThis.document.createElement(tag);
  for (const [name, value] of Object.entries(attrs || {})) {
    if (value === false || value === null || value === undefined) continue;
    el.setAttribute(name, value === true ? "" : value);
  }
  if (typeof children === "string") el.textContent = children;
  else (children || []).forEach((child) => el.appendChild(typeof child === "string" ? new TextNode(child) : child));
  return el;
}

// --- Window ----------------------------------------------------------------

class FakeTimers {
  constructor() {
    this.now = 0;
    this.nextId = 1;
    this.queue = [];
  }

  set(fn, ms) {
    const id = this.nextId++;
    this.queue.push({ id, at: this.now + Math.max(0, Number(ms) || 0), fn });
    return id;
  }

  clear(id) {
    this.queue = this.queue.filter((t) => t.id !== id);
  }

  // Runs every timer due within `ms` from now, in order (timers they set
  // that also fall due run too).
  advance(ms) {
    const until = this.now + ms;
    for (;;) {
      this.queue.sort((p, q) => p.at - q.at || p.id - q.id);
      const next = this.queue[0];
      if (!next || next.at > until) break;
      this.queue.shift();
      this.now = next.at;
      next.fn();
    }
    this.now = until;
  }

  pending() {
    return this.queue.map((t) => t.at - this.now).sort((p, q) => p - q);
  }
}

class FakeLocation {
  constructor(href) {
    this.assigned = [];
    this.reloads = 0;
    this.set(href);
  }

  set(href) {
    const url = new URL(href, this.url ? this.url.href : "http://localhost/");
    this.url = url;
  }

  get href() { return this.url.href; }
  get origin() { return this.url.origin; }
  get pathname() { return this.url.pathname; }
  get search() { return this.url.search; }
  get hash() { return this.url.hash; }

  assign(href) {
    this.assigned.push(href);
  }

  replace(href) {
    this.assigned.push(href);
  }

  reload() {
    this.reloads += 1;
  }

  toString() {
    return this.href;
  }
}

class FakeHistory {
  constructor(win) {
    this.win = win;
    this.entries = [{ state: null, href: win.location.href }];
    this.index = 0;
  }

  get state() {
    return this.entries[this.index].state;
  }

  get length() {
    return this.entries.length;
  }

  pushState(state, _title, url) {
    if (url !== undefined && url !== null) this.win.location.set(url);
    this.entries = this.entries.slice(0, this.index + 1);
    this.entries.push({ state: structuredClone(state), href: this.win.location.href });
    this.index += 1;
  }

  replaceState(state, _title, url) {
    if (url !== undefined && url !== null) this.win.location.set(url);
    this.entries[this.index] = { state: structuredClone(state), href: this.win.location.href };
  }

  go(delta) {
    const to = this.index + delta;
    if (to < 0 || to >= this.entries.length || !delta) return;
    this.index = to;
    this.win.location.set(this.entries[to].href);
    const state = this.entries[to].state;
    // A browser fires popstate after the current task.
    queueMicrotask(() => this.win.dispatchEvent(new FakeEvent("popstate", { state: structuredClone(state), bubbles: false })));
  }

  back() {
    this.go(-1);
  }

  forward() {
    this.go(1);
  }
}

class FakeStorage {
  constructor() {
    this.map = new Map();
  }

  get length() { return this.map.size; }
  key(n) { return Array.from(this.map.keys())[n] ?? null; }
  getItem(k) { return this.map.has(k) ? this.map.get(k) : null; }
  setItem(k, v) { this.map.set(k, String(v)); }
  removeItem(k) { this.map.delete(k); }
  clear() { this.map.clear(); }
}

// A Response-like answer for the fake fetch.
export function jsonResponse(body, status) {
  status = status || 200;
  return { ok: status >= 200 && status < 300, status, json: async () => body, text: async () => JSON.stringify(body) };
}

class FakeWindow extends EventTargetish {
  constructor(href) {
    super();
    this.timers = new FakeTimers();
    this.location = new FakeLocation(href || "http://localhost/");
    this.history = new FakeHistory(this);
    this.localStorage = new FakeStorage();
    this.devicePixelRatio = 1;
    this.innerWidth = 1280;
    this.innerHeight = 800;
    // fetch(url, options) -> a Promise of a Response; tests replace
    // `fetchHandler` (url, options) -> body | Response.
    this.fetchCalls = [];
    this.fetchHandler = () => { throw new Error("fakedom: no fetch handler set"); };
    this.alerts = [];
    this.confirmAnswer = true;
  }

  setTimeout(fn, ms) { return this.timers.set(fn, ms); }
  clearTimeout(id) { this.timers.clear(id); }
  setInterval(fn, ms) {
    const timers = this.timers;
    let id;
    const tick = () => { id = timers.set(tick, ms); fn(); };
    id = timers.set(tick, ms);
    return id;
  }
  clearInterval(id) { this.timers.clear(id); }
  requestAnimationFrame(fn) { return this.timers.set(() => fn(this.timers.now), 16); }
  cancelAnimationFrame(id) { this.timers.clear(id); }
  matchMedia(query) { return { matches: false, media: query, addEventListener() {}, removeEventListener() {} }; }
  getComputedStyle(el) { return { getPropertyValue: () => "", color: (el && el.style && el.style.color) || "rgb(0, 0, 0)", overflowX: "visible" }; }
  alert(text) { this.alerts.push(text); }
  confirm() { return this.confirmAnswer; }

  fetch(url, options) {
    this.fetchCalls.push({ url: String(url), options: options || {} });
    try {
      const answer = this.fetchHandler(String(url), options || {});
      return Promise.resolve(answer).then((a) => (a && typeof a.json === "function" ? a : jsonResponse(a)));
    } catch (err) {
      return Promise.reject(err);
    }
  }
}

// Puts a fresh window and document on globalThis and returns the window
// (its `document` is window.document).
export function installDom(href) {
  const win = new FakeWindow(href);
  const doc = new FakeDocument();
  doc.defaultView = win;
  win.document = doc;
  win.CustomEvent = class extends FakeEvent {
    constructor(type, init) {
      super(type, init);
      this.detail = init ? init.detail : undefined;
    }
  };
  win.Event = FakeEvent;
  const globals = {
    window: win, document: doc, location: win.location, history: win.history, localStorage: win.localStorage,
    CustomEvent: win.CustomEvent, Event: FakeEvent,
    fetch: (url, options) => win.fetch(url, options),
    getComputedStyle: (el) => win.getComputedStyle(el),
    requestAnimationFrame: (fn) => win.requestAnimationFrame(fn),
    cancelAnimationFrame: (id) => win.cancelAnimationFrame(id),
    matchMedia: (q) => win.matchMedia(q),
  };
  for (const [name, value] of Object.entries(globals)) {
    Object.defineProperty(globalThis, name, { value, configurable: true, writable: true });
  }
  return win;
}

// Lets pending promise callbacks (fetch answers, popstate) run.
export async function settle(rounds) {
  for (let n = 0; n < (rounds || 10); n++) await new Promise((resolve) => setImmediate(resolve));
}

// The static/ directory's URL, for importing a page module.
export const STATIC_URL = new URL("../../html/static/", import.meta.url);

let importCount = 0;

// Loads a page script afresh (its top-level code runs again), as a page
// would load it with `?v=...`. A module (one with import/export) is
// imported with a fresh query; a plain script (no import or export,
// which node would load once as CommonJS and cache) is run as a
// function body against the current globals, in strict mode like a
// module.
export async function importPage(name) {
  importCount += 1;
  const url = new URL(name + "?fakedom=" + importCount, STATIC_URL);
  const source = readFileSync(new URL(name, STATIC_URL), "utf8");
  if (/^\s*(import|export)\b|\bimport\.meta\b|\bawait import\(/m.test(source)) return import(url.href);
  new Function('"use strict";\n' + source + "\n//# sourceURL=" + url.href)();
  return {};
}
