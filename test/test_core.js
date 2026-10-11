// run: jsc src/core.js test/test_core.js   (from the project root)
const G = globalThis.GelArt;
const load = (f, circ) => { const r = G.parseSequences(readFile(f))[0]; r.circular = circ; return r; };
const lam = load('data/lambda_J02459.fa', false), puc = load('data/pUC19_L09137.fa', true), pbr = load('data/pBR322_J01749.fa', true);
let fail = 0;
const eq = (label, got, exp) => { const ok = JSON.stringify(got) === JSON.stringify(exp); if (!ok) fail++; print((ok ? 'ok   ' : 'FAIL ') + label + (ok ? '' : '\n   got ' + JSON.stringify(got) + '\n   exp ' + JSON.stringify(exp))); };
eq('lambda EcoRI', G.digest(lam, ['EcoRI']), [21226, 7421, 5804, 5643, 4878, 3530]);
eq('lambda BamHI', G.digest(lam, ['BamHI']), [16841, 7233, 6770, 6527, 5626, 5505]);
// J02459 is WILD-TYPE lambda: 6 HindIII sites. NEB's lambda is cI857 Sam7, whose ind1 mutation adds a 7th
// site splitting 6682 into 6557 + 125 — so the NEB datasheet fragment list does not apply to this accession.
eq('lambda HindIII (wt J02459, 6 sites)', G.digest(lam, ['HindIII']), [23130, 9416, 6682, 4361, 2322, 2027, 564]);
eq('pUC19 DraI', G.digest(puc, ['DraI']), [1975, 692, 19]);
eq('pUC19 PvuII', G.digest(puc, ['PvuII']), [2364, 322]);
for (const e of ['EcoRI', 'BamHI', 'HindIII', 'PstI', 'SalI']) eq('pBR322 ' + e + ' single cutter', G.digest(pbr, [e]), [4361]);
const hae = G.digest(pbr, ['HaeIII']); eq('pBR322 HaeIII n=22, top 587/540/502', [hae.length, hae[0], hae[1], hae[2]], [22, 587, 540, 502]);
eq('pBR322 BsaI (type IIS) sums to length', G.digest(pbr, ['BsaI']).reduce((a, b) => a + b, 0), 4361);
eq('circular digest fragments sum to plasmid length', G.digest(puc, ['HaeIII']).reduce((a, b) => a + b, 0), 2686);
eq('enzyme table parsed, no dup names', [G.ENZYMES.length > 150, new Set(G.ENZYMES.map(e => e.name)).size === G.ENZYMES.length], [true, true]);
eq('all sites valid', G.ENZYMES.filter(e => !(e.c1 >= 0) || !/^[ACGTRYSWKMBDHVN]+$/.test(e.site)).map(e => e.name), []);
// back-compat: default model has p=1, d0=0, so distanceMm must be unchanged by the timecourse extension
const m1 = G.modelFor(1);
eq('linear model still linear in time', Math.abs(G.distanceMm(1000, m1, 5, 120) - 2 * G.distanceMm(1000, m1, 5, 60)) < 1e-9, true);
eq('field from electrode gap, not gel length', G.fieldVcm(100, 15).toFixed(2), '6.67');
print(fail ? fail + ' FAILED' : 'all passed');
