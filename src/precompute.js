/* Precompute the shape catalog at the bench conditions, so the app opens on a real gel instead of a spinner.
 * jsc src/core.js src/precompute.js > data/catalog.json
 */
'use strict';
const G = globalThis.GelArt;

const load = (f, circ, name) => { const r = G.parseSequences(readFile(f))[0]; r.circular = circ; r.name = name; return r; };
const puc = load('data/pUC19_L09137.fa', true, 'pUC19');
const pbr = load('data/pBR322_J01749.fa', true, 'pBR322');

// Gen091226 rig: 36.8 cm electrode gap, 4 V/cm = 147 V, 1.2% TAE, one comb row with 100 mm of running room.
const GEL = { pct: 1.2, vcm: 4, top: 6, len: 106, rows: 20, loadNg: 500, visMin: 2, visFull: 20 };
const LANES = 20;
const TIMES = [45, 60, 90, 120, 150, 180];
const PLASMIDS = [puc, pbr];

const cat = G.catalog({
  plasmids: PLASMIDS, enzymes: G.ENZYMES, gel: GEL, lanes: LANES, times: TIMES, maxEnz: 3,
  variantOpts: { offsetStride: 2, heightStride: 3, widths: [LANES, 12, 9] }
});

cat.then(res => {
  const out = {
    gel: GEL, lanes: LANES, times: TIMES, electrodeGapCm: 36.8,
    plasmids: PLASMIDS.map(p => ({ name: p.name, len: p.seq.length, circular: p.circular })),
    variants: res.variants, combos: res.combos,
    shapes: res.entries.filter(e => e.bestOfShape).map(e => ({
      name: e.name, size: e.size, y0: e.y0, time: e.time,
      quality: +e.quality.toFixed(4), on: e.on,
      target: e.target.map(col => col.join('')),
      lanes: e.lanes.map(l => l.empty ? null : { p: l.options[0].plasmid, e: l.options[0].enzymes, T: l.options[0].temp, err: +l.options[0].err.toFixed(3) })
    }))
  };
  print(JSON.stringify(out));
}).catch(e => print('ERR ' + e + '\n' + e.stack));
