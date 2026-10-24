// html/static/quantity.js
//
// The browser's Quantity (UX.128): the mirror of planetgen/physics/quantity.py (GEN.569). A number of one
// kind in one unit of the registry (unitregistry.js), kept as given and answered in every other unit of its
// kind. A script that holds a value next to its unit holds it in one of these and asks for the unit it
// needs, `quantity.as("km")`, instead of multiplying by a factor of its own.
//
// The truth is never converted: every other form is the truth scaled by the two units' factors
// (`convert`), kept after the first read, so no drift builds up. A value must be finite and inside its
// kind's bounds (`check`), and a unit of another kind is refused with a message naming both, as in Python.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { KINDS, UNITS, check, convert, toBase, unit } = await import(`./unitregistry.js${VERSION_QUERY}`);

const unitOfKind = (kind, key) => {
  const row = unit(key);
  if (row.kind !== kind) throw new Error(`${row.symbol} (${key}) measures ${row.kind}, not ${kind}`);
  return row;
};

export class Quantity {
  // `kind` is a key of KINDS, `unit` a key of UNITS whose kind is `kind`; `value` is finite and inside the
  // kind's bounds. Throws an Error otherwise (Python raises ValueError).
  constructor(kind, value, unitKey) {
    if (!(kind in KINDS)) throw new Error(`unknown kind ${kind}`);
    unitOfKind(kind, unitKey);
    this.kind = kind;
    this.unit = unitKey;
    this.value = check(Number(value), unitKey);
    this.forms = null;
  }

  // This quantity's number in `key` (a unit of this kind): the truth itself for its own unit, otherwise
  // the truth scaled by the two units' factors, kept after the first read.
  as(key) {
    if (key === this.unit) return this.value;
    unitOfKind(this.kind, key);
    this.forms ??= new Map();
    if (!this.forms.has(key)) this.forms.set(key, convert(this.value, this.unit, key));
    return this.forms.get(key);
  }

  // The quantity in its kind's base unit, as a number (offsets applied).
  base() {
    return toBase(this.value, this.unit);
  }

  // Changes the truth to `value` of `key` (this quantity's own unit when omitted); every form follows.
  set(value, key = this.unit) {
    unitOfKind(this.kind, key);
    this.value = check(Number(value), key);
    this.unit = key;
    this.forms = null;
  }

  toJSON() {
    return { kind: this.kind, value: this.value, unit: this.unit };
  }

  static fromJSON({ kind, value, unit: key }) {
    return new Quantity(kind, value, key);
  }
}

// A Quantity of the unit `key`, its kind the unit's own: `held("km", 384400)`.
export const held = (key, value) => new Quantity(unit(key).kind, value, key);

export { UNITS };
