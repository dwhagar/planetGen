// static/apiversion.js -- API.31: a page's script cannot set `User-Agent`, so every request it makes to this site's
// own `/api/` carries the API version in the `X-PlanetGen-API-Version` header instead, the header the server already
// sends back (API.22) and logs (API.15). The page names both in
// `<meta name="planetgen-api-version" content="N" data-header="X-PlanetGen-API-Version">`.

/** Whether `address` (a string, a URL or a Request's url) is this site's own `/api/`, from a page at `base`. */
export function isApiAddress(address, base) {
  try {
    const url = new URL(address, base);
    return url.origin === new URL(base).origin && url.pathname.startsWith("/api/");
  } catch (error) {
    return false;
  }
}

/** Wraps `target.fetch` and, when present, `target.XMLHttpRequest` so a request to the site's `/api/` sets `header`. */
export function install(target, base, header, version) {
  if (!header || !version) return;
  if (typeof target.fetch === "function") {
    const original = target.fetch.bind(target);
    target.fetch = (input, init) => {
      const address = typeof input === "string" || input instanceof URL ? String(input) : input.url;
      if (!isApiAddress(address, base)) return original(input, init);
      const headers = new Headers((init && init.headers) || (typeof input === "object" && input.headers) || undefined);
      headers.set(header, version);
      return original(input, { ...init, headers });
    };
  }
  const Request = target.XMLHttpRequest;
  if (Request && Request.prototype && typeof Request.prototype.open === "function") {
    const open = Request.prototype.open;
    Request.prototype.open = function (method, address, ...rest) {
      const result = open.call(this, method, address, ...rest);
      if (isApiAddress(String(address), base)) this.setRequestHeader(header, version);
      return result;
    };
  }
}

if (typeof document !== "undefined") {
  const meta = document.querySelector('meta[name="planetgen-api-version"]');
  if (meta) install(window, window.location.href, meta.dataset.header, meta.getAttribute("content"));
}
