// run: jsc src/buffers.js src/core.js test/test_reuse.js
const G = globalThis.GelArt;
const db = JSON.parse(readFile('data/plasmiddb.json'));
const load = (f, c, n) => { const r = G.parseSequences(readFile(f))[0]; r.circular = c; r.name = n; return r; };
const PL = [load('data/pUC19_L09137.fa', true, 'pUC19'), load('data/pBR322_J01749.fa', true, 'pBR322')]
  .concat(db.map(d => ({ name: d.name, seq: d.seq, circular: true })));
const cat = JSON.parse(readFile('data/catalog.json'));
const ENZ = G.panelEnzymes('lab');
let fail = 0;
const ok = (c, l, x) => { if (!c) fail++; print((c ? 'ok   ' : 'FAIL ') + l + (x ? '  [' + x + ']' : '')); };

// Wave was the pathological case: quality went to -82% at tol 2.0 before the bound existed.
const SHAPES = ['Wave', 'Triangle', 'Heart', 'Double helix', 'Smiley'];
const TOLS = [0, 0.25, 0.5, 1, 2, 5, 100];
const MINQ = 0.6;
(async () => {
  for (const nm of SHAPES) {
    const s = cat.shapes.filter(x => x.name === nm)[0];
    const g = Object.assign({}, cat.gel); g.time = s.time; delete g.model;
    const base = await G.solve({
      plasmids: PL, enzymes: ENZ, target: s.target.map(c => c.split('').map(Number)),
      gel: g, times: [s.time], maxEnz: 2, minAct: 50
    });
    const bq = 1 - base.lanes.reduce((a, l) => a + (l.empty ? 0 : l.options[0].err), 0) / s.on;
    const bd = new Set(base.lanes.filter(l => !l.empty).map(l => l.options[0].plasmid + ':' + l.options[0].enzymes.join('+'))).size;
    // At tol 0 every column keeps its own best digest, so this IS the base per-column floor.
    const baseWorstCol = base.reuseAt(0, 0, MINQ).worstColumnQuality;
    let worst = 1, worstTol = null, maxD = 0;
    for (const tol of TOLS) {
      const r = base.reuseAt(0, tol, MINQ);
      if (r.quality < worst) { worst = r.quality; worstTol = tol; }
      maxD = Math.max(maxD, r.digests);
      ok(r.quality >= 0, `${nm} tol ${tol}: quality not negative`, (r.quality * 100).toFixed(1) + '%');
      // The guarantee is per column: reuse never pushes a column below minQ, and never below its own
      // achievable best if that was already under minQ. Comparing against the OVERALL quality was wrong -
      // a low-weight column can sit far below the weighted average in the fidelity-optimal solve too.
      ok(r.worstColumnQuality >= Math.min(MINQ, baseWorstCol) - 1e-9,
         `${nm} tol ${tol}: no column pushed below the floor`,
         'worst col ' + (r.worstColumnQuality * 100).toFixed(0) + '%, base worst ' + (baseWorstCol * 100).toFixed(0) + '%');
      ok(r.digests <= bd, `${nm} tol ${tol}: never more reactions than the base`, r.digests + ' vs ' + bd);
    }
    const z = base.reuseAt(0, 0, MINQ);
    ok(Math.abs(z.quality - bq) < 0.02, `${nm} tol 0 matches the fidelity-optimal solve`,
       (z.quality * 100).toFixed(1) + '% vs ' + (bq * 100).toFixed(1) + '%');
    print(`     ${nm}: base ${bd} digests @ ${(bq * 100).toFixed(1)}%  |  worst over all tolerances ${(worst * 100).toFixed(1)}% (tol ${worstTol})`);
  }
  print(fail ? '\n' + fail + ' FAILED' : '\nreuse bound holds at every tolerance');
})().catch(e => print('ERR ' + e + '\n' + e.stack));
