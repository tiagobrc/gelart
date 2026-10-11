/* Gel Art calibration — generate a bench sheet, then fit the migration model to the measured bands.
 *
 *   jsc src/core.js calibration/calibrate.js -- template
 *   jsc src/core.js calibration/calibrate.js -- fit calibration/measurements.csv
 *
 * Edit CONFIG to match the rig BEFORE generating the template.
 */
'use strict';
const G = globalThis.GelArt;

const CONFIG = {
  // Rig as recorded in Gen091226: big box, 14.5 in = 36.8 cm electrode to electrode, 25 cm tray.
  // House rule (Genotyping protocol standard): specify V/cm, never volts; stay at or below 5 V/cm.
  electrodeGapCm: 36.8,
  vcmTarget: 4,            // volts are derived from this, not typed in
  agarosePct: 1.2,
  buffer: '1x TAE',
  ladder: 'NEB 1 kb Plus',
  ladderLanes: [1, 10, 20], // left / centre / right of the 20-well comb, to expose smiling
  usableMm: 100,           // running room below the wells; 1.2-1.5% needs >=70 mm, 100 mm is comfortable
  topMm: 6,                // art grid starts this far below the wells
  gridRows: 20,
  // Both optional. With neither, the migration model still fits completely from the ladder alone - but
  // hold out one TIMEPOINT instead (mark that photo's rows role=check), or the only error estimate left
  // is the fit residual, which flatters itself.
  validation: [],   // e.g. [{ lane: 2, plasmid: 'A_5', enzymes: ['SphI'] }]
  // Without bands that FAIL to appear, visMin is unfittable: set this to run a two-fold series past the floor.
  dilution: null,   // e.g. { lane0: 4, plasmid: 'A_5', enzymes: ['SphI'], series: [250,125,62,31,16,8] }
  smallestBp: 100,         // sets t_max: the run stops when this band reaches usableMm
  fractions: [0.2, 0.4, 0.6, 0.8, 1.0],   // nominal only; photographing every 15 min is equivalent
  bandThicknessMm: 2.5     // house rule: "a band is ~2-3 mm wide, so 3.2 mm is about one band width"
};
CONFIG.volts = Math.round(CONFIG.vcmTarget * CONFIG.electrodeGapCm);

const args = (typeof arguments !== 'undefined' ? arguments : []);
const mode = args[0] || 'template';
const vcm = G.fieldVcm(CONFIG.volts, CONFIG.electrodeGapCm);
const prior = G.modelFor(CONFIG.agarosePct);
const pad = (s, n) => String(s).padStart(n);

function timepoints() {
  const tmax = 60 * CONFIG.usableMm / G.distanceMm(CONFIG.smallestBp, prior, vcm, 60);
  return CONFIG.fractions.map(f => ({ frac: f, t: Math.round(tmax * f) }));
}

// ---------------------------------------------------------------- template
function template() {
  const sizes = G.LADDERS[CONFIG.ladder].sizes, tps = timepoints();
  const out = [
    '# Gel Art migration calibration — measurements',
    '# Fill the mm column at the bench/on the screen. Leave a row blank if the band is off the gel or unreadable.',
    '# mm = leading (anode-side) edge of the well -> CENTRE of the band. Be consistent; a constant offset is fitted out as d0.',
    '# t_min = POWER-ON minutes (stopwatch pauses whenever the supply is off).',
    '# pred_mm is the UNCALIBRATED prior, shown only so you know roughly where to look. Do not copy it into mm.',
    '#',
    '# t_min below is NOMINAL. The prior disagrees with Gen091226 by roughly 2x, so do not follow it as a',
    '# schedule - photograph on a fixed interval, stop on the ladder, and overwrite t_min with the minutes',
    '# you actually recorded. Measure 5-6 photos spread across the run; more frames than that add work,',
    '# not information. (Optional: -- replan <this file> after the first photo rescales the nominal times.)',
    '#',
    `# volts=${CONFIG.volts}`,
    `# electrode_gap_cm=${CONFIG.electrodeGapCm}`,
    `# vcm=${vcm.toFixed(3)}`,
    `# agarose_pct=${CONFIG.agarosePct}`,
    `# buffer=${CONFIG.buffer}`,
    `# ladder=${CONFIG.ladder}`,
    `# usable_mm=${CONFIG.usableMm}`,
    `# top_mm=${CONFIG.topMm}`,
    `# grid_rows=${CONFIG.gridRows}`,
    '# gel_temp_start_C=',
    '# gel_temp_end_C=',
    'role,lane,label,len_bp,frac,t_min,mm,pred_mm,vis,mass_ng,notes'
  ];
  for (const tp of tps) for (const lane of CONFIG.ladderLanes) for (const len of sizes) {
    const d = G.distanceMm(len, prior, vcm, tp.t);
    out.push(`fit,${lane},ladder,${len},${tp.frac},${tp.t},,${d > CONFIG.usableMm ? 'offgel' : d.toFixed(1)},,,`);
  }
  // Fragment sizes are computed here from the sequence, not transcribed.
  let db = null;
  try { db = JSON.parse(readFile('data/plasmiddb.json')); } catch (e) { db = null; }
  const plasmid = nm => { const d = db && db.filter(x => x.name === nm)[0]; return d ? { name: d.name, seq: d.seq, circular: true } : null; };

  if (!CONFIG.validation.length && !CONFIG.dilution) {
    out.push('# --- ladder only. No digest is held out, so hold out a TIMEPOINT instead: when the rows are');
    out.push('#     filled, change role from fit to check on ONE photo. Fitting the rest and predicting that');
    out.push('#     photo tests exactly what the model is for - an unmeasured run time.');
    out.push('# --- visMin / visFull are NOT measured by this run and keep their default 2 / 20 ng.');
  }
  if (CONFIG.validation.length) out.push('# --- validation: real digests held OUT of the fit.');
  for (const v of CONFIG.validation) {
    const pl = plasmid(v.plasmid);
    if (!pl) { out.push(`# ${v.plasmid} not found in data/plasmiddb.json`); continue; }
    const frags = G.digest(pl, v.enzymes);
    for (const tp of tps) for (const L of frags) {
      const d = G.distanceMm(L, prior, vcm, tp.t);
      out.push(`check,${v.lane},${v.plasmid}/${v.enzymes.join('+')},${L},${tp.frac},${tp.t},,${d > CONFIG.usableMm ? 'offgel' : d.toFixed(1)},,,`);
    }
  }

  const dp = CONFIG.dilution ? plasmid(CONFIG.dilution.plasmid) : null;
  if (dp) {
    out.push('# --- dilution series: score vis 0 (invisible) / 0.5 (faint) / 1 (bright) on ONE photo only.');
    out.push('# Use photo 1 or 2: EtBr cast in the gel migrates to the cathode, so late frames lose stain at the bottom.');
    const frags = G.digest(dp, CONFIG.dilution.enzymes), total = frags.reduce((a, b) => a + b, 0);
    CONFIG.dilution.series.forEach((ng, i) => {
      for (const L of frags) {
        const mass = (L / total * ng).toFixed(2);
        out.push(`vis,${CONFIG.dilution.lane0 + i},${CONFIG.dilution.plasmid}/${CONFIG.dilution.enzymes.join('+')} ${ng}ng,${L},,,,,${''},${mass},`);
      }
    });
  }
  return out.join('\n') + '\n';
}

// ---------------------------------------------------------------- fit
function parseCsv(text) {
  const meta = {}, rows = [];
  let header = null;
  for (const raw of text.split('\n')) {
    const line = raw.trim();
    if (!line) continue;
    if (line[0] === '#') { const m = line.slice(1).trim().match(/^([a-z_]+)=(.*)$/); if (m) meta[m[1]] = m[2]; continue; }
    const cells = line.split(',').map(c => c.trim());
    if (!header) { header = cells; continue; }
    const r = {}; header.forEach((h, i) => (r[h] = cells[i] === undefined ? '' : cells[i]));
    rows.push(r);
  }
  return { meta, rows };
}

function bar(v, scale) {
  const n = Math.round(Math.abs(v) * scale);
  return v < 0 ? pad('', Math.max(0, 12 - n)) + '#'.repeat(n) + '|' : pad('', 12) + '|' + '#'.repeat(n);
}

function fit(path) {
  const { meta, rows } = parseCsv(readFile(path));
  const V = meta.vcm ? +meta.vcm : vcm;
  const pct = meta.agarose_pct ? +meta.agarose_pct : CONFIG.agarosePct;
  const usable = meta.usable_mm ? +meta.usable_mm : CONFIG.usableMm;
  const top = meta.top_mm ? +meta.top_mm : CONFIG.topMm;
  const gridRows = meta.grid_rows ? +meta.grid_rows : CONFIG.gridRows;

  const num = r => ({ len: +r.len_bp, t: +r.t_min, mm: +r.mm, lane: +r.lane, label: r.label });
  const usableRow = r => r.len_bp && r.t_min && r.mm !== '' && isFinite(+r.mm);
  const fitPts = rows.filter(r => r.role === 'fit' && usableRow(r)).map(num);
  const checkPts = rows.filter(r => r.role === 'check' && usableRow(r)).map(num);
  const visPts = rows.filter(r => r.vis !== '' && r.mass_ng !== '' && isFinite(+r.vis) && isFinite(+r.mass_ng))
    .map(r => ({ mass: +r.mass_ng, vis: +r.vis }));

  print(`GEL      ${pct}% ${meta.buffer || CONFIG.buffer}, ${meta.volts || CONFIG.volts} V over ${meta.electrode_gap_cm || CONFIG.electrodeGapCm} cm = ${V.toFixed(2)} V/cm`);
  print(`POINTS   ${fitPts.length} fit, ${checkPts.length} validation, ${visPts.length} intensity`);
  if (fitPts.length < 12) { print('\nNot enough measured bands to fit (need >= 12, ideally >= 40).'); return; }

  const res = G.fitTimecourse(fitPts, V, { pct });
  const { simple, full, recommended, drift } = res;
  const chosen = recommended === full ? 'drift + offset' : 'linear in time';
  print('');
  print(`STAGE 1  linear in time     rmse ${simple.rmse.toFixed(2)} mm`);
  print(`STAGE 2  + drift + offset   rmse ${full.rmse.toFixed(2)} mm   p ${full.p.toFixed(3)}  d0 ${full.d0 >= 0 ? '+' : ''}${full.d0.toFixed(2)} mm`);
  print(`CHOSEN   ${chosen}${recommended === full ? '  (drift model earned it: >10% better)' : '  (drift model did not earn the extra parameters)'}`);

  print('\nDRIFT DIAGNOSTIC — mean residual of the LINEAR model per timepoint');
  print('  a monotonic trend means migration is NOT linear in time; flat means it is');
  const scale = 4;
  for (const d of drift) print(`  ${pad(d.t, 4)} min ${pad(d.meanResid.toFixed(2), 6)} mm  ${bar(d.meanResid, scale)}`);

  // The fit is Huber-robust, so a mis-measured band surfaces here instead of dragging the model.
  fitPts.forEach(p => { p.e = G.distanceMm(p.len, recommended, V, p.t) - p.mm; });
  const OUT = 3.0;
  const outliers = fitPts.filter(p => Math.abs(p.e) > OUT).sort((a, b) => Math.abs(b.e) - Math.abs(a.e));
  const clean = fitPts.filter(p => Math.abs(p.e) <= OUT);
  const cleanRmse = Math.sqrt(clean.reduce((s, p) => s + p.e * p.e, 0) / Math.max(1, clean.length));
  const line = p => `  ${pad(p.len, 6)} bp @ ${pad(p.t, 3)} min lane ${pad(p.lane, 2)}: ${p.e >= 0 ? '+' : ''}${p.e.toFixed(2)} mm`;

  if (outliers.length) {
    print(`\nOUTLIERS  ${outliers.length} band(s) off by more than ${OUT} mm — re-measure these before trusting the rest`);
    outliers.slice(0, 6).forEach(p => print(line(p)));
  }
  print(`\nRESIDUALS rmse ${recommended.rmse.toFixed(2)} mm all points / ${cleanRmse.toFixed(2)} mm excluding outliers`);
  clean.slice().sort((a, b) => Math.abs(b.e) - Math.abs(a.e)).slice(0, 4).forEach(p => print(line(p)));

  // lane-to-lane spread at matched (len, t) = smiling; outliers excluded so one bad band cannot fake a smile
  const byKey = {};
  clean.forEach(p => { const k = p.len + '@' + p.t; (byKey[k] = byKey[k] || []).push(p.mm); });
  const spreads = Object.values(byKey).filter(v => v.length > 1).map(v => Math.max(...v) - Math.min(...v)).sort((a, b) => a - b);
  if (spreads.length) {
    const med = spreads[Math.floor(spreads.length / 2)], p90 = spreads[Math.min(spreads.length - 1, Math.floor(spreads.length * 0.9))];
    print(`\nSMILING   median lane-to-lane spread ${med.toFixed(2)} mm, 90th pct ${p90.toFixed(2)} mm (${spreads.length} matched sets)`);
    if (p90 > 1.5) print('  > 1.5 mm: drop the voltage, or keep the art in the central lanes only.');
  }

  if (checkPts.length) {
    let s = 0;
    for (const p of checkPts) s += (G.distanceMm(p.len, recommended, V, p.t) - p.mm) ** 2;
    print(`\nVALIDATION (held out) rmse ${Math.sqrt(s / checkPts.length).toFixed(2)} mm over ${checkPts.length} bands  <-- the honest number`);
  } else {
    print('\nVALIDATION  none — fit rmse flatters itself. Add `check` rows from a real digest.');
  }

  const mr = G.maxRows({ len: usable + top, top, rows: gridRows }, cleanRmse, CONFIG.bandThicknessMm);
  print(`\nGRID      usable ${mr.usable.toFixed(0)} mm / rmse ${cleanRmse.toFixed(2)} mm / band ${CONFIG.bandThicknessMm} mm`);
  print(`          max honest rows ${mr.rows} (pitch ${mr.pitch.toFixed(2)} mm)`);
  print(`          ${gridRows}-row grid has pitch ${mr.currentPitch.toFixed(2)} mm — ${mr.rows >= gridRows ? 'OK' : 'TOO FINE, drop to ' + mr.rows + ' rows'}`);

  if (visPts.length >= 4) {
    const fv = G.fitVisibility(visPts);
    print(`\nVISIBILITY ${visPts.length} scored bands: ${fv.span.lo} invisible, ${fv.span.mid} partial, ${fv.span.hi} bright`);
    print(fv.visMinDetermined
      ? `           visMin  ${fv.visMin.toFixed(2)} ng — faintest band the imager picks up`
      : '           visMin  UNDETERMINED — needs at least one invisible and one partial band. Add a 2-fold dilution series down to ~0.25 ng.');
    print(fv.visFullDetermined
      ? `           visFull ${fv.visFull.toFixed(1)} ng — fully bright`
      : '           visFull UNDETERMINED — needs a saturated band and a partial one.');
  } else if (visPts.length) {
    print(`\nVISIBILITY only ${visPts.length} scored bands — need >= 4 spanning invisible to bright.`);
  }

  const m = recommended;
  print('\nMODEL — paste into the app');
  print(`  {"pct":${pct},"vcm":${V.toFixed(3)},"muS":${m.muS.toFixed(4)},"muL":${m.muL.toFixed(4)},"gamma":${m.gamma.toFixed(0)},"p":${m.p.toFixed(4)},"d0":${m.d0.toFixed(3)},"tref":60}`);
}

// After photo 1: correct the remaining timepoints for however wrong the prior turned out to be.
function replan(path) {
  const { meta, rows } = parseCsv(readFile(path));
  const V = meta.vcm ? +meta.vcm : vcm;
  const pct = meta.agarose_pct ? +meta.agarose_pct : CONFIG.agarosePct;
  const usable = meta.usable_mm ? +meta.usable_mm : CONFIG.usableMm;
  const pr = G.modelFor(pct);
  const pts = rows.filter(r => r.role === 'fit' && r.mm !== '' && isFinite(+r.mm) && r.len_bp && r.t_min)
    .map(r => ({ len: +r.len_bp, t: +r.t_min, mm: +r.mm }));
  if (!pts.length) { print('No measured bands yet — fill the mm column for the first photo.'); return; }

  const t0 = Math.min(...pts.map(p => p.t));
  const first = pts.filter(p => p.t === t0);
  const ratios = first.map(p => p.mm / G.distanceMm(p.len, pr, V, p.t)).sort((a, b) => a - b);
  const scale = ratios[Math.floor(ratios.length / 2)];
  const spread = ratios[ratios.length - 1] / ratios[0];

  print(`PHOTO 1   ${first.length} bands at ${t0} min power-on, ${V.toFixed(2)} V/cm, ${pct}%`);
  print(`SCALE     observed / predicted = ${scale.toFixed(2)}x (median) — the gel runs ${scale > 1 ? 'FASTER' : 'SLOWER'} than the prior`);
  if (spread > 1.6) print(`  WARNING  ratio varies ${ratios[0].toFixed(2)}-${ratios[ratios.length - 1].toFixed(2)} across sizes: the prior's SHAPE is wrong,`);
  if (spread > 1.6) print('           not just its speed. A single rescale is rough here — keep all 5 photos.');

  const tmaxOld = 60 * usable / G.distanceMm(CONFIG.smallestBp, pr, V, 60);
  const tmax = tmaxOld / scale;
  print(`\nt_max     ${tmax.toFixed(0)} min for ${CONFIG.smallestBp} bp to reach ${usable} mm (prior said ${tmaxOld.toFixed(0)} min)`);
  print('REVISED TIMEPOINTS — overwrite t_min on the unmeasured rows');
  for (const f of CONFIG.fractions) {
    const t = Math.round(tmax * f);
    print(`  ${pad((f * 100).toFixed(0) + '%', 5)}  ${pad(t, 4)} min${t <= t0 ? '   <- already past' : ''}`);
  }
  print('\nJudge the endpoint by the ladder, not the dye front. Stop before the smallest band you care about runs off.');
}

if (mode === 'template') print(template());
else if (mode === 'fit') fit(args[1] || 'calibration/measurements.csv');
else if (mode === 'replan') replan(args[1] || 'calibration/measurements.csv');
else print('usage: -- template | -- replan <csv> | -- fit <csv>');
