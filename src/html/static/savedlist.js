// html/static/savedlist.js
//
// The saved-list store (UX.129): one list of saved entries of one kind (bookmarks now, saved courses for NAV.18),
// behind a backend that may answer later than it is asked. The page scripts that show a list (bookmarks.js) sit on
// top of this and know nothing of where the list is kept.
//
// A backend is {available(), load(), save(entries), subscribe(fn)}:
// - available(): false when nothing can be saved here (storage blocked, no database named); then load() gives [].
// - load(): the stored entries, or a promise of them.
// - save(entries): true when the whole list was kept, or a promise of that.
// - subscribe(fn): calls fn() whenever the stored list may have changed somewhere else (another tab, another copy of
//   the script on the page); the store reads it again.
// `localBackend` is the browser's own localStorage, one key per database, the format the bookmarks have always had
// (a JSON array of entries). An account backend (USR.7, NAV.19) is another object of the same shape.
//
// A store is made from a backend and what makes an entry of this kind:
// - keyOf(entry): the string an entry is known by, or null when it is not an entry of this kind;
// - build(entry): the form that is saved (default: the entry as given);
// - max: the most entries there can be;
// - words: {invalid, full, blocked}, the reasons `add` gives (full may be a function of max);
// - stamp: a function giving what to record with every new entry (the galaxy instance stamp UX.46 needs, as
//   `entry.stamp`), or nothing.
//
// Reads (`list`, `find`) answer at once from what was last loaded; `ready` settles when the first load has. Writes
// (`add`, `rename`, `remove`) read the stored list again first, so a change made elsewhere is never overwritten, then
// save, and answer when the backend has. `onChange` is called after every change, here or elsewhere.

export const CHANGE_EVENT = "planetgen-saved-change";

// The list kept in `localStorage` under `prefix + suffix()` (suffix() is "" when no database is named: then
// nothing is kept). `storage()` gives the Storage object, and throws or gives null where there is none.
export function localBackend({ prefix, suffix, storage }) {
  const getStorage = storage || function () { return window.localStorage; };

  function key() {
    const name = suffix();
    return name ? prefix + name : null;
  }

  return {
    available: function () { return key() !== null; },
    load: function () {
      const k = key();
      if (!k) return [];
      try {
        const parsed = JSON.parse(getStorage().getItem(k) || "[]");
        return Array.isArray(parsed) ? parsed : [];
      } catch (e) {
        return [];
      }
    },
    save: function (entries) {
      const k = key();
      if (!k) return false;
      try {
        getStorage().setItem(k, JSON.stringify(entries));
      } catch (e) {
        return false;
      }
      window.dispatchEvent(new CustomEvent(CHANGE_EVENT, { detail: { key: k } }));
      return true;
    },
    subscribe: function (fn) {
      window.addEventListener(CHANGE_EVENT, function (event) {
        if (!event.detail || event.detail.key === key()) fn();
      });
      window.addEventListener("storage", function (event) {
        if (event.key === null || event.key === key()) fn();
      });
    },
  };
}

export function createStore({ backend, keyOf, build, max, words, stamp }) {
  const reasons = Object.assign({
    invalid: "That can't be saved.",
    full: function (most) { return "The list is full (" + most + " at most); delete one first."; },
    blocked: "This browser isn't letting the site save.",
  }, words || {});
  const listeners = [];
  let entries = [];
  let writes = Promise.resolve();

  function clean(raw) {
    return (Array.isArray(raw) ? raw : []).filter(function (entry) { return keyOf(entry) !== null; }).slice(0, max);
  }

  function notify() {
    listeners.slice().forEach(function (fn) { fn(); });
  }

  // Reads the stored list again. True when it differs from what was held.
  async function reload() {
    const next = clean(await backend.load());
    const changed = JSON.stringify(next) !== JSON.stringify(entries);
    entries = next;
    return changed;
  }

  // One write at a time, in the order asked.
  function enqueue(task) {
    const run = writes.then(task, task);
    writes = run.catch(function () {});
    return run;
  }

  function indexOf(list, key) {
    return list.findIndex(function (entry) { return keyOf(entry) === key; });
  }

  async function keep(next) {
    if (!(await backend.save(next))) return false;
    entries = next;
    notify();
    return true;
  }

  const ready = reload();
  backend.subscribe(function () {
    reload().then(function (changed) { if (changed) notify(); });
  });

  return {
    ready: ready,
    available: function () { return backend.available(); },
    // Every entry, oldest first.
    list: function () { return entries.slice(); },
    find: function (key) {
      const index = indexOf(entries, key);
      return index < 0 ? null : entries[index];
    },
    // Saves `entry`: {ok, entry} or {ok: false, reason}. Already saved is ok, as it was.
    add: function (entry) {
      return enqueue(async function () {
        const key = keyOf(entry);
        if (key === null) return { ok: false, reason: reasons.invalid };
        await reload();
        const index = indexOf(entries, key);
        if (index >= 0) return { ok: true, entry: entries[index] };
        if (entries.length >= max) {
          return { ok: false, reason: typeof reasons.full === "function" ? reasons.full(max) : reasons.full };
        }
        const saved = build ? build(entry) : Object.assign({}, entry);
        const mark = stamp ? stamp() : null;
        if (mark != null) saved.stamp = mark;
        if (!(await keep(entries.concat([saved])))) return { ok: false, reason: reasons.blocked };
        return { ok: true, entry: saved };
      });
    },
    // True when the entry known by `key` now has `name`.
    rename: function (key, name) {
      return enqueue(async function () {
        const trimmed = String(name).trim();
        await reload();
        const index = indexOf(entries, key);
        if (index < 0 || !trimmed) return false;
        const next = entries.slice();
        next[index] = Object.assign({}, next[index], { name: trimmed });
        return keep(next);
      });
    },
    // True when the entry known by `key` is gone.
    remove: function (key) {
      return enqueue(async function () {
        await reload();
        const index = indexOf(entries, key);
        if (index < 0) return false;
        return keep(entries.filter(function (_entry, n) { return n !== index; }));
      });
    },
    // Calls `fn` after every change, here or elsewhere.
    onChange: function (fn) {
      listeners.push(fn);
    },
  };
}
