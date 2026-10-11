// run: jsc src/core.js test/test_calib.js
const G = globalThis.GelArt;
let fail = 0;
const ok = (label, cond, extra) => { if (!cond) fail++; print((cond ? 'ok   ' : 'FAIL ') + label + (extra ? '  [' + extra + ']' : '')); };
let s = 12345; const rnd = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; };
const gauss = sd => sd * Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
const LADDER = G.LADDERS['NEB 1 kb Plus'].sizes;
const GEL_LEN = 60, VCM = 5;

function synth(truth, times, noise) {
  const pts = [];
  for (const t of times) for (const len of LADDER) {
    const mm = G.distanceMm(len, truth, VCM, t) + gauss(noise);
    if (mm > 1 && mm < GEL_LEN) pts.push({ len, t, mm });     // bands still on the gel
  }
  return pts;
}
// 1. drifting gel: truth has p < 1 and a well offset
const truth = { muS: 2.7, muL: 0.5, gamma: 4200, p: 0.88, d0: 1.1, tref: 60 };
const times = [14, 27, 41, 55, 69];
const pts = synth(truth, times, 0.3);
const fit = G.fitTimecourse(pts, VCM, { pct: 1 });
ok('timecourse fit converges below noise floor', fit.recommended.rmse < 0.45, 'rmse ' + fit.recommended.rmse.toFixed(3) + ' mm, n=' + pts.length);
ok('drift model beats linear on a drifting gel', fit.full.rmse < fit.simple.rmse * 0.9, 'linear ' + fit.simple.rmse.toFixed(2) + ' vs drift ' + fit.full.rmse.toFixed(2) + ' mm');
ok('recovers p', Math.abs(fit.recommended.p - truth.p) < 0.06, 'p ' + fit.recommended.p.toFixed(3) + ' vs ' + truth.p);
// held-out extrapolation: a time never fitted
let worst = 0;
for (const len of LADDER) { const d = G.distanceMm(len, truth, VCM, 48); if (d < GEL_LEN) worst = Math.max(worst, Math.abs(G.distanceMm(len, fit.recommended, VCM, 48) - d)); }
ok('predicts an unmeasured run time', worst < 1.0, 'worst error ' + worst.toFixed(2) + ' mm at 48 min');
// 2. a well-behaved linear gel must NOT be talked into a drift model
const lin = { muS: 2.9, muL: 0.55, gamma: 5000, p: 1, d0: 0, tref: 60 };
const fitLin = G.fitTimecourse(synth(lin, times, 0.3), VCM, { pct: 1 });
ok('linear gel keeps the simple model', fitLin.recommended.p === 1, 'p=' + fitLin.recommended.p.toFixed(3) + ', linear rmse ' + fitLin.simple.rmse.toFixed(2));
ok('drift diagnostic flat when there is no drift', Math.max(...fitLin.drift.map(d => Math.abs(d.meanResid))) < 0.5,
   fitLin.drift.map(d => d.t + 'min:' + d.meanResid.toFixed(2)).join(' '));
ok('drift diagnostic trends when there IS drift', Math.max(...fit.drift.map(d => Math.abs(d.meanResid))) > 0.5,
   fit.drift.map(d => d.t + 'min:' + d.meanResid.toFixed(2)).join(' '));
// 3. one badly mismeasured band must not wreck the fit (Huber)
const dirty = pts.slice(); dirty[10] = { len: dirty[10].len, t: dirty[10].t, mm: dirty[10].mm + 12 };
ok('robust to one outlier band', Math.abs(G.fitTimecourse(dirty, VCM, { pct: 1 }).recommended.gamma - fit.recommended.gamma) / fit.recommended.gamma < 0.15,
   'gamma ' + G.fitTimecourse(dirty, VCM, { pct: 1 }).recommended.gamma.toFixed(0) + ' vs clean ' + fit.recommended.gamma.toFixed(0));
// 4. visibility
const vt = [{ mass: 0.5, vis: 0 }, { mass: 1, vis: 0 }, { mass: 3, vis: 0.1 }, { mass: 8, vis: 0.45 }, { mass: 15, vis: 0.9 }, { mass: 25, vis: 1 }, { mass: 60, vis: 1 }];
const fv = G.fitVisibility(vt);
ok('visibility floor fitted', fv.visMin > 0.5 && fv.visMin < 6 && fv.visFull > fv.visMin, 'visMin ' + fv.visMin.toFixed(2) + ' ng, visFull ' + fv.visFull.toFixed(1) + ' ng');
// 5. how many rows the geometry can honestly carry
for (const r of [0.3, 0.8, 1.5]) { const mr = G.maxRows({ len: GEL_LEN, top: 6, rows: 16 }, r, 1.5); print(`     rmse ${r} mm -> max ${mr.rows} rows (pitch ${mr.pitch.toFixed(2)} mm, usable ${mr.usable} mm; 16-row grid pitch ${mr.currentPitch.toFixed(2)} mm)`); }
// 6. bench timing table from the uncalibrated prior
const m = G.modelFor(1), usable = 55, small = 100;
const tmax = 60 * usable / G.distanceMm(small, m, VCM, 60);
print(`\n  1% TAE, ${VCM} V/cm, ${usable} mm usable: 100 bp reaches the end at t=${tmax.toFixed(0)} min`);
print('  frac  t(min) | ' + [10000, 3000, 1000, 500, 100].map(l => String(l).padStart(6)).join(''));
for (const f of [0.2, 0.4, 0.6, 0.8, 1.0]) {
  const t = tmax * f;
  print(`  ${(f * 100).toFixed(0).padStart(4)}%  ${t.toFixed(0).padStart(6)} | ` + [10000, 3000, 1000, 500, 100].map(l => { const d = G.distanceMm(l, m, VCM, t); return (d > usable ? 'off' : d.toFixed(1)).padStart(6); }).join(''));
}
print(fail ? '\n' + fail + ' FAILED' : '\nall passed');
