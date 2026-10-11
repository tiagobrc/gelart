// run: jsc src/core.js test/test_solver.js
const G = globalThis.GelArt;
const load = (f, circ) => { const r = G.parseSequences(readFile(f))[0]; r.circular = circ; return r; };
const lam = load('data/lambda_J02459.fa', false), puc = load('data/pUC19_L09137.fa', true), pbr = load('data/pBR322_J01749.fa', true);
lam.seq = lam.seq.slice(0, 37588) + 'T' + lam.seq.slice(37589);   // ind1
const times = []; for (let t = 30; t <= 120; t += 10) times.push(t);
const show = (target, lanes, gel, pl) => {
  const pred = lanes.map(l => l.empty ? new Array(gel.rows).fill(0) : G.rowIntensity(G.laneBands(G.digest(pl[l.options[0].plasmid], l.options[0].enzymes), gel), gel));
  for (let r = 0; r < gel.rows; r++) print('  ' + target.map(c => c[r] ? '#' : '.').join('') + '   ' + pred.map(c => c[r] > 0.66 ? '#' : c[r] > 0.33 ? '+' : '.').join(''));
};
async function run(label, pl, loadNg, enzymes, maxEnz) {
  const gel = { pct: 1, vcm: 5, top: 6, len: 60, rows: 16, loadNg };
  let t0 = Date.now();
  const target = G.placeBitmap(G.PARAMETRIC['Double helix'](12, 10), 12, 16, 3);
  const res = await G.solve({ plasmids: pl, enzymes, target, gel, times, maxEnz });
  print(`\n== ${label}: helix, ${res.combos} combos, ${Date.now() - t0} ms, best time ${res.time} min`);
  show(target, res.lanes, Object.assign({}, gel, { time: res.time }), pl);
  res.lanes.forEach((l, i) => l.empty || print(`  lane ${i + 1}: ${l.options[0].enzymes.join('+')} err ${l.options[0].err.toFixed(2)}`));
  t0 = Date.now();
  const cat = await G.catalog({ plasmids: pl, enzymes, gel, lanes: 12, times, maxEnz });
  print(`-- catalog: ${cat.variants} variants, ${cat.combos} combos, ${Date.now() - t0} ms`);
  cat.entries.filter(e => e.bestOfShape).forEach(e => print(`  ${(e.quality * 100).toFixed(0).padStart(3)}%  ${e.name} ${e.size} @row ${e.y0}, ${e.time} min`));
  const top = cat.entries[0]; show(top.target, top.lanes, Object.assign({}, gel, { time: top.time }), pl);
}
const common = G.ENZYMES.filter(e => e.group === 'six' || e.group === 'deg');
(async () => {
  await run('pUC19, 6-cutters+degenerate', [puc], 500, common, 3);
  await run('pUC19, all enzymes', [puc], 500, G.ENZYMES, 3);
  await run('pUC19 + pBR322, all enzymes', [puc, pbr], 500, G.ENZYMES, 3);
  await run('lambda, 6-cutters+degenerate', [lam], 1000, common, 3);
})().catch(e => print('ERR ' + e + '\n' + e.stack));
