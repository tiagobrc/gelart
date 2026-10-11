/* Catalog over the real PlasmidDB. jsc src/core.js src/precompute_db.js > data/catalog_db.json */
'use strict';
const G = globalThis.GelArt;
const db = JSON.parse(readFile('data/plasmiddb.json'));
const load = (f, circ, name) => { const r = G.parseSequences(readFile(f))[0]; r.circular = circ; r.name = name; return r; };
const PLASMIDS = [load('data/pUC19_L09137.fa', true, 'pUC19'), load('data/pBR322_J01749.fa', true, 'pBR322')]
  .concat(db.map(d => ({ name: d.name, seq: d.seq, circular: d.circular })));
const GEL = { pct: 1.2, vcm: 4, top: 6, len: 106, rows: 20, loadNg: 500, visMin: 2, visFull: 20 };
const LANES = 20, TIMES = [60, 90, 120, 150, 180];
const t0 = Date.now();
G.catalog({
  plasmids: PLASMIDS, enzymes: G.panelEnzymes('lab'), gel: GEL, lanes: LANES, times: TIMES, maxEnz: 2, minAct: 50,
  variantOpts: { offsetStride: 2, heightStride: 3, widths: [LANES, 12, 9] }
}).then(res => {
  print(JSON.stringify({
    gel: GEL, lanes: LANES, times: TIMES, electrodeGapCm: 36.8, maxEnz: 2, panel: 'lab', minAct: 50,
    plasmids: PLASMIDS.map(p => ({ name: p.name, len: p.seq.length, circular: p.circular })),
    variants: res.variants, combos: res.combos, ms: Date.now() - t0,
    shapes: res.entries.filter(e => e.bestOfShape).map(e => ({
      name: e.name, size: e.size, y0: e.y0, time: e.time, quality: +e.quality.toFixed(4), on: e.on,
      target: e.target.map(c => c.join('')),
      lanes: e.lanes.map(l => l.empty ? null : {
        p: l.options[0].plasmid, e: l.options[0].enzymes, T: l.options[0].temp, err: +l.options[0].err.toFixed(3),
        b: l.options[0].buffer, bw: l.options[0].bufWorst, bs: l.options[0].bufStar, pr: l.options[0].products
      })
    }))
  }));
}).catch(e => print('ERR ' + e + '\n' + e.stack));
