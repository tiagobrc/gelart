/* Gel Art core -- digest engine, migration model, lane solver. No DOM; runs in a browser or any JS shell. */
(function (root) {
'use strict';

// ---------------------------------------------------------------- enzymes
// name|site|incubation  degC|flags   -- palindromic sites mark the top-strand cut with ^;
// non-palindromic sites give (c1/c2): top / bottom strand cut, counted from the site start.
// flags: dam = blocked by Dam methylation, damo = blocked only when the site overlaps GATC.
const ENZYME_TABLE = {
six: `AatII|GACGT^C AclI|AA^CGTT Acc65I|G^GTACC AfeI|AGC^GCT AflII|C^TTAAG AgeI|A^CCGGT ApaI|GGGCC^C|25
ApaLI|G^TGCAC AscI|GG^CGCGCC AseI|AT^TAAT AsiSI|GCGAT^CGC AvrII|C^CTAGG BamHI|G^GATCC BclI|T^GATCA|37|dam
BglII|A^GATCT BmtI|GCTAG^C BsiWI|C^GTACG BspEI|T^CCGGA|37|damo BspHI|T^CATGA|37|damo BsrGI|T^GTACA
BssHII|G^CGCGC|50 BstBI|TT^CGAA|65 BstZ17I|GTA^TAC ClaI|AT^CGAT|37|damo DraI|TTT^AAA EagI|C^GGCCG
EcoRI|G^AATTC EcoRV|GAT^ATC FseI|GGCCGG^CC FspI|TGC^GCA HindIII|A^AGCTT HpaI|GTT^AAC KasI|G^GCGCC
KpnI|GGTAC^C MfeI|C^AATTG MluI|A^CGCGT MscI|TGG^CCA NaeI|GCC^GGC NarI|GG^CGCC NcoI|C^CATGG NdeI|CA^TATG
NgoMIV|G^CCGGC NheI|G^CTAGC NotI|GC^GGCCGC NruI|TCG^CGA|37|damo NsiI|ATGCA^T PacI|TTAAT^TAA PciI|A^CATGT
PmeI|GTTT^AAAC PmlI|CAC^GTG PsiI|TTA^TAA PspOMI|G^GGCCC PstI|CTGCA^G PvuI|CGAT^CG PvuII|CAG^CTG
SacI|GAGCT^C SacII|CCGC^GG SalI|G^TCGAC SbfI|CCTGCA^GG ScaI|AGT^ACT SfoI|GGC^GCC SmaI|CCC^GGG|25
SnaBI|TAC^GTA SpeI|A^CTAGT SphI|GCATG^C SrfI|GCCC^GGGC SspI|AAT^ATT StuI|AGG^CCT SwaI|ATTT^AAAT|25
XbaI|T^CTAGA|37|damo XhoI|C^TCGAG XmaI|C^CCGGG ZraI|GAC^GTC`,
deg: `AccI|GT^MKAC AflIII|A^CRYGT AhdI|GACNNN^NNGTC AleI|CACNN^NNGTG AlwNI|CAGNNN^CTG ApoI|R^AATTY
AvaI|C^YCGRG AvaII|G^GWCC BanI|G^GYRCC BanII|GRGCY^C BglI|GCCNNNN^NGGC BlpI|GC^TNAGC BsaBI|GATNN^NNATC|60|damo
BsaHI|GR^CGYC BsaJI|C^CNNGG|60 BsaWI|W^CCGGW|60 BsiEI|CGRY^CG|60 BsiHKAI|GWGCW^C|65 BslI|CCNNNNN^NNGG|55
Bsp1286I|GDGCH^C BsrFI|R^CCGGY BstAPI|GCANNNN^NTGC|60 BstEII|G^GTNACC BstXI|CCANNNNN^NTGG BstYI|R^GATCY|60
Bsu36I|CC^TNAGG BtgI|C^CRYGG DraIII|CACNNN^GTG EaeI|Y^GGCCR EcoNI|CCTNN^NNNAGG EcoO109I|RG^GNCCY
HaeII|RGCGC^Y HincII|GTY^RAC MslI|CAYNN^NNRTG MspA1I|CMG^CKG MwoI|GCNNNNN^NNGC|60 NlaIV|GGN^NCC NspI|RCATG^Y
PflMI|CCANNNN^NTGG PpuMI|RG^GWCCY PshAI|GACNN^NNGTC RsrII|CG^GWCCG SexAI|A^CCWGGT SfcI|C^TRYAG
SfiI|GGCCNNNN^NGGCC|50 SmlI|C^TYRAG|55 StyI|C^CWWGG Tth111I|GACN^NNGTC|65 XcmI|CCANNNNN^NNNNTGG XmnI|GAANN^NNTTC`,
four: `AluI|AG^CT BfaI|C^TAG BstNI|CC^WGG|60 BstUI|CG^CG|60 CviQI|G^TAC|25 DdeI|C^TNAG DpnII|^GATC|37|dam
Fnu4HI|GC^NGC HaeIII|GG^CC HhaI|GCG^C HinP1I|G^CGC HinfI|G^ANTC HpaII|C^CGG Hpy188I|TCN^GA|37|damo
Hpy188III|TC^NNGA HpyCH4III|ACN^GT HpyCH4IV|A^CGT HpyCH4V|TG^CA MboI|^GATC|37|dam MluCI|^AATT MseI|T^TAA
MspI|C^CGG NciI|CC^SGG NlaIII|CATG^ RsaI|GT^AC Sau3AI|^GATC Sau96I|G^GNCC ScrFI|CC^NGG TaqI|T^CGA|65|damo
TfiI|G^AWTC|65 TseI|G^CWGC|65 Tsp45I|^GTSAC|65`,
iis: `AcuI|CTGAAG(22/20) AlwI|GGATC(9/10)|37|dam BbsI|GAAGAC(8/12) BbvI|GCAGC(13/17) BccI|CCATC(9/10)
BceAI|ACGGC(17/19) BciVI|GTATCC(12/11) BmrI|ACTGGG(11/10) BpmI|CTGGAG(22/20) BpuEI|CTTGAG(22/20)
BsaI|GGTCTC(7/11) BseRI|GAGGAG(16/14) BseYI|CCCAGC(1/5) BsgI|GTGCAG(22/20) BsmAI|GTCTC(6/10)|55
BsmBI|CGTCTC(7/11)|55 BsmFI|GGGAC(15/19)|65 BsmI|GAATGC(7/5)|65 BspMI|ACCTGC(10/14) BsrBI|CCGCTC(3/3)
BsrDI|GCAATG(8/6)|65 BsrI|ACTGG(6/4)|65 BssSI|CACGAG(1/5) BtgZI|GCGATG(16/20)|60 BtsI|GCAGTG(8/6)
EarI|CTCTTC(7/10) EciI|GGCGGA(17/15) Esp3I|CGTCTC(7/11) FauI|CCCGC(9/11)|55 FokI|GGATG(14/18)
HgaI|GACGC(10/15) HphI|GGTGA(13/12)|37|damo HpyAV|CCTTC(11/10) MboII|GAAGA(13/12)|37|damo MlyI|GAGTC(10/10)
MmeI|TCCRAC(26/24) PaqCI|CACCTGC(11/15) PleI|GAGTC(9/10) SapI|GCTCTTC(8/11) SfaNI|GCATC(10/14)`
};
const GROUP_LABELS = { six: '6+ bp cutters', deg: 'Degenerate / interrupted', four: '4\u20135 bp frequent cutters', iis: 'Type IIS / non-palindromic' };

const IUPAC = { A: 'A', C: 'C', G: 'G', T: 'T', R: 'AG', Y: 'CT', S: 'CG', W: 'AT', K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT' };
const COMP = { A: 'T', C: 'G', G: 'C', T: 'A', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', V: 'B', D: 'H', H: 'D', N: 'N' };

function revcomp(s) { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] || 'N'; return o; }
function siteRegex(site) {
  let p = '';
  for (const ch of site) p += ch === 'N' ? '.' : IUPAC[ch].length === 1 ? ch : '[' + IUPAC[ch] + ']';
  return new RegExp('(?=' + p + ')', 'g');
}

function buildEnzymes() {
  const out = [];
  for (const group of Object.keys(ENZYME_TABLE)) {
    for (const tok of ENZYME_TABLE[group].split(/\s+/)) {
      if (!tok) continue;
      const [name, rawSite, temp, flag] = tok.split('|');
      let site, c1, c2;
      const m = rawSite.match(/^([A-Z]+)\((-?\d+)\/(-?\d+)\)$/);
      if (m) { site = m[1]; c1 = +m[2]; c2 = +m[3]; }
      else { c1 = rawSite.indexOf('^'); site = rawSite.replace('^', ''); c2 = site.length - c1; }
      const rc = revcomp(site);
      out.push({
        name, site, c1, c2, group, temp: temp ? +temp : 37,
        dam: flag === 'dam' ? 'blocked' : flag === 'damo' ? 'overlap' : null,
        palindromic: rc === site, display: rawSite,
        fwd: siteRegex(site), rev: rc === site ? null : siteRegex(rc)
      });
    }
  }
  return out.sort((a, b) => a.name.localeCompare(b.name));
}
const ENZYMES = buildEnzymes();
const ENZYME_BY_NAME = Object.fromEntries(ENZYMES.map(e => [e.name, e]));

// ---------------------------------------------------------------- buffers
// Activity data comes from src/buffers.js (parsed from NEB's published chart). Absent = no entry, never guessed.
const BUF_DATA = root.GELART_BUFFERS || { bufs: [], e: {} };
const BUFFERS = BUF_DATA.bufs.length ? BUF_DATA.bufs : ['r1.1', 'r2.1', 'r3.1', 'rCutSmart'];
ENZYMES.forEach(e => {
  const d = BUF_DATA.e[e.name];
  e.sup = d ? d[0] : null;            // buffer the enzyme ships in ("U" = its own unique buffer)
  e.act = d ? d[1] : null;            // percent activity in each of the four standard buffers
  e.starIn = d ? d[2] : null;         // star-activity risk flagged in that buffer
  e.product = d && d[3] ? d[3] : e.name;  // the catalogue item actually charted (HF / v2 variant)
  e.supplement = (BUF_DATA.supp && BUF_DATA.supp[e.name]) || null;  // SAM etc; omitting it fails the digest
});
// Bit per buffer in which this enzyme reaches minAct. No data means no bits: an enzyme whose behaviour
// is unknown is never silently paired with another under a buffer constraint.
function bufMask(e, minAct) {
  if (!e.act) return 0;
  let m = 0;
  for (let i = 0; i < 4; i++) if (e.act[i] >= minAct) m |= 1 << i;
  return m;
}
// Best single buffer for a digest: maximise the WORST enzyme's activity, break ties toward rCutSmart.
function sharedBuffer(names, minAct) {
  const es = names.map(n => ENZYME_BY_NAME[n]).filter(Boolean);
  if (!es.length || es.some(e => !e.act)) return null;
  let best = null;
  for (let i = 0; i < 4; i++) {
    let worst = 100, star = false;
    for (const e of es) { if (e.act[i] < worst) worst = e.act[i]; if (e.starIn && e.starIn[i]) star = true; }
    if (worst < (minAct || 0)) continue;
    const score = worst * 10 + (i === 3 ? 1 : 0);
    if (!best || score > best.score) best = { buffer: BUFFERS[i], idx: i, worst, star, score };
  }
  return best;
}

// ---------------------------------------------------------------- panels
// A digest the freezer cannot perform is not a recipe. Restrict the search to enzymes actually on hand.
// The VICTORA_LAB SnapGene enzyme set (2021) UNION the enzymes that appear in 2023-2026 bench notes but
// are missing from it, since that file has gone stale. 88 of the 95 lab names exist in the table above;
// absent and therefore unsearchable: BfuAI, BmgBI, BsaAI, DpnI, PluTI, PspXI, SgrAI.
const LAB_PANEL = ['Acc65I','AccI','AfeI','AgeI','AhdI','AleI','ApaI','AscI','AseI','AvrII','BamHI','BbsI','BglI','BglII','BmtI','BpmI','BsaI','BsaWI','BseRI','BseYI','BsgI','BsiWI','BslI','BsmBI','BsmFI','BsmI','BspEI','BspHI','BspMI','BsrBI','BsrDI','BsrGI','BssHII','BstEII','BstXI','BtgI','ClaI','EagI','EcoRI','EcoRV','Fnu4HI','FseI','FspI','HincII','HindIII','HpaI','KasI','KpnI','MfeI','MluI','MscI','NaeI','NcoI','NdeI','NgoMIV','NheI','NotI','NruI','NsiI','PacI','PmeI','PmlI','PpuMI','PshAI','PsiI','PspOMI','PstI','PvuI','SacI','SacII','SalI','SapI','SbfI','ScaI','SexAI','SfoI','SnaBI','SpeI','SphI','SrfI','SspI','SwaI','XbaI','XcmI','XhoI','XmaI','XmnI','ZraI'];
const LAB_SET = {}; LAB_PANEL.forEach(n => (LAB_SET[n] = 1));

// Supplier synonyms with IDENTICAL recognition site AND cut positions. Enzymes that merely share a site
// but cut it differently (PluTI vs KasI, BfuAI vs BspMI) are deliberately NOT aliased - the fragments differ.
const ALIASES = { kpn2i: 'BspEI', bspdi: 'ClaI', tspmi: 'XmaI' };

const PANELS = {
  lab:      { label: 'Victora lab stock', test: e => LAB_SET[e.name] },
  all:      { label: 'Every enzyme in the table', test: null },
  common:   { label: '6+ bp and degenerate', test: e => e.group === 'six' || e.group === 'deg' },
  six:      { label: '6+ bp cutters only', test: e => e.group === 'six' },
  frequent: { label: '4-5 bp frequent cutters', test: e => e.group === 'four' },
  iis:      { label: 'Type IIS only', test: e => e.group === 'iis' }
};
// spec: a panel key, or an explicit array of enzyme names.
function panelEnzymes(spec) {
  if (Array.isArray(spec)) { const s = new Set(spec); return ENZYMES.filter(e => s.has(e.name)); }
  const p = PANELS[spec];
  return !p || !p.test ? ENZYMES : ENZYMES.filter(p.test);
}
// Parse a pasted freezer list: any separator, ignores unknown names, reports them so nothing is silently dropped.
function parsePanel(text) {
  const want = String(text || '').split(/[^A-Za-z0-9.]+/).filter(Boolean);
  const byLower = {};
  ENZYMES.forEach(e => (byLower[e.name.toLowerCase()] = e.name));
  const found = [], missing = [];
  want.forEach(w => {
    const k = w.toLowerCase();
    const hit = byLower[k] || byLower[(ALIASES[k] || '').toLowerCase()];
    if (hit) { if (found.indexOf(hit) < 0) found.push(hit); } else missing.push(w);
  });
  return { names: found, missing: missing };
}

// ---------------------------------------------------------------- sequences
function cleanSeq(s) { return s.toUpperCase().replace(/U/g, 'T').replace(/[^ACGTRYSWKMBDHVN]/g, ''); }

function parseSequences(text, fallbackName) {
  const out = [];
  if (/^LOCUS\s/m.test(text)) {
    for (const rec of text.split(/^\/\/\s*$/m)) {
      const locus = rec.match(/^LOCUS\s+(\S+)(.*)$/m), origin = rec.search(/^ORIGIN/m);
      if (!locus || origin < 0) continue;
      const seq = cleanSeq(rec.slice(origin + 6).replace(/[\d\s]/g, ''));
      if (seq) out.push({ name: locus[1], seq, circular: /circular/i.test(locus[2]) ? true : /linear/i.test(locus[2]) ? false : null });
    }
  } else if (/^\s*>/.test(text)) {
    for (const rec of text.split(/^\s*>/m)) {
      if (!rec.trim()) continue;
      const nl = rec.indexOf('\n'), head = (nl < 0 ? rec : rec.slice(0, nl)).trim();
      const seq = cleanSeq(nl < 0 ? '' : rec.slice(nl));
      if (seq) out.push({ name: head.split(/\s+/)[0] || fallbackName, seq, circular: /circular/i.test(head) ? true : /linear/i.test(head) ? false : null });
    }
  } else {
    const seq = cleanSeq(text);
    if (seq) out.push({ name: fallbackName || 'sequence', seq, circular: null });
  }
  return out;
}

// SnapGene .dna: packets of [type:1][length:4 BE][data]; packet 0 = topology flags byte + ASCII sequence.
function parseSnapGene(buffer, fallbackName) {
  const dv = new DataView(buffer), bytes = new Uint8Array(buffer);
  if (bytes[0] !== 0x09) return [];
  let off = 0;
  while (off + 5 <= bytes.length) {
    const type = bytes[off], len = dv.getUint32(off + 1);
    if (type === 0) {
      let s = '';
      for (let i = off + 6; i < off + 5 + len; i++) s += String.fromCharCode(bytes[i]);
      return [{ name: fallbackName || 'snapgene', seq: cleanSeq(s), circular: (bytes[off + 5] & 1) === 1 }];
    }
    off += 5 + len;
  }
  return [];
}

// ---------------------------------------------------------------- digestion
// A cut at k sits between bases k-1 and k (0-based, top strand).
function findCuts(seq, enz, circular) {
  const n = seq.length, ext = circular ? seq + seq.slice(0, 40) : seq, cuts = new Set();
  const scan = (re, offset) => {
    re.lastIndex = 0;
    let m;
    while ((m = re.exec(ext)) !== null) {
      if (m.index >= n) break;
      let k = m.index + offset;
      if (circular) cuts.add(((k % n) + n) % n);
      else if (k > 0 && k < n) cuts.add(k);
      re.lastIndex = m.index + 1;
    }
  };
  scan(enz.fwd, enz.c1);
  if (enz.rev) scan(enz.rev, enz.site.length - enz.c2);
  return Array.from(cuts).sort((a, b) => a - b);
}

function fragmentsFromCuts(cuts, n, circular) {
  const f = [];
  if (!cuts.length) return circular ? [] : [n];
  if (circular) {
    for (let i = 1; i < cuts.length; i++) f.push(cuts[i] - cuts[i - 1]);
    f.push(n - cuts[cuts.length - 1] + cuts[0]);
  } else {
    let prev = 0;
    for (const c of cuts) { f.push(c - prev); prev = c; }
    f.push(n - prev);
  }
  return f;
}

function mergeCuts(lists) {
  const all = [].concat(...lists).sort((a, b) => a - b), out = [];
  for (let i = 0; i < all.length; i++) if (i === 0 || all[i] !== all[i - 1]) out.push(all[i]);
  return out;
}

function digest(plasmid, enzymeNames) {
  const lists = enzymeNames.map(nm => findCuts(plasmid.seq, ENZYME_BY_NAME[nm], plasmid.circular));
  return fragmentsFromCuts(mergeCuts(lists), plasmid.seq.length, plasmid.circular).sort((a, b) => b - a);
}

function cutCounts(plasmid, enzymes) {
  const o = {};
  for (const e of enzymes || ENZYMES) o[e.name] = findCuts(plasmid.seq, e, plasmid.circular).length;
  return o;
}

// ---------------------------------------------------------------- migration
// Van Winkle, Beheshti & Rill (Electrophoresis 2002):  1/mu(L) = 1/muL - (1/muL - 1/muS) * exp(-L/gamma)
// mu in 1e-4 cm^2/V/s. Values are approximate (TAE, ~5 V/cm) -- calibrate with a ladder for a real run.
const AGAROSE_MODELS = [
  { pct: 0.7, muS: 3.1, muL: 0.85, gamma: 8000 },
  { pct: 1.0, muS: 2.9, muL: 0.55, gamma: 5000 },
  { pct: 1.5, muS: 2.5, muL: 0.30, gamma: 3000 },
  { pct: 2.0, muS: 2.1, muL: 0.15, gamma: 2000 }
];
function modelFor(pct) {
  const M = AGAROSE_MODELS;
  if (pct <= M[0].pct) return M[0];
  if (pct >= M[M.length - 1].pct) return M[M.length - 1];
  for (let i = 1; i < M.length; i++) if (pct <= M[i].pct) {
    const a = M[i - 1], b = M[i], t = (pct - a.pct) / (b.pct - a.pct);
    return { pct, muS: a.muS + t * (b.muS - a.muS), muL: a.muL + t * (b.muL - a.muL), gamma: a.gamma + t * (b.gamma - a.gamma) };
  }
}
function mobility(len, m) { return 1 / (1 / m.muL - (1 / m.muL - 1 / m.muS) * Math.exp(-len / m.gamma)); }
// The field is set by the ELECTRODE GAP, not the gel length -- a mini tank runs ~15 cm of buffer over a 7 cm gel,
// so computing V/cm from gel length overstates the field roughly twofold.
function fieldVcm(volts, electrodeGapCm) { return volts / electrodeGapCm; }
// E.t_eff, with t_eff = tref.(t/tref)^p absorbing buffer exhaustion / heating drift. p = 1 is strictly linear in time.
function travel(m, vcm, minutes) {
  const tref = m.tref || 60, p = m.p == null ? 1 : m.p;
  return vcm * tref * Math.pow(Math.max(0, minutes) / tref, p) * 60 * 1e-3;
}
function distanceMm(len, m, vcm, minutes) { return (m.d0 || 0) + mobility(len, m) * travel(m, vcm, minutes); }
function visibility(massNg, gel) { return Math.max(0, Math.min(1, (massNg - gel.visMin) / (gel.visFull - gel.visMin))); }

function resolveGel(g) {
  const gel = Object.assign({ pct: 1, vcm: 5, time: 60, top: 6, len: 60, rows: 16, visMin: 2, visFull: 20, loadNg: 500 }, g);
  gel.model = g.model || modelFor(gel.pct);
  gel.rowH = (gel.len - gel.top) / gel.rows;
  return gel;
}

function laneBands(fragments, g) {
  const gel = resolveGel(g), total = fragments.reduce((a, b) => a + b, 0);
  return fragments.map(len => {
    const mass = len / total * gel.loadNg, dist = distanceMm(len, gel.model, gel.vcm, gel.time);
    return { len, mass, dist, vis: visibility(mass, gel), offGel: dist > gel.len, aboveImage: dist < gel.top };
  });
}

function splat(row, vis, dist, gel) {
  const pos = (dist - gel.top) / gel.rowH - 0.5, i = Math.floor(pos), fr = pos - i;
  if (i >= 0 && i < gel.rows) row[i] += vis * (1 - fr);
  if (i + 1 >= 0 && i + 1 < gel.rows) row[i + 1] += vis * fr;
}
function rowIntensity(bands, g) {
  const gel = resolveGel(g), row = new Float32Array(gel.rows);
  for (const b of bands) if (b.vis > 0 && !b.offGel && !b.aboveImage) splat(row, b.vis, b.dist, gel);
  return Array.from(row, v => Math.min(1, v));
}

// ---------------------------------------------------------------- ladder calibration
function nelderMead(f, x0, iters) {
  const n = x0.length;
  let S = [x0.slice()];
  for (let i = 0; i < n; i++) { const p = x0.slice(); p[i] += 0.25; S.push(p); }
  let F = S.map(f);
  for (let it = 0; it < iters; it++) {
    const ord = F.map((v, i) => i).sort((a, b) => F[a] - F[b]);
    S = ord.map(i => S[i]); F = ord.map(i => F[i]);
    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) c[j] += S[i][j] / n;
    const at = t => c.map((cj, j) => cj + t * (S[n][j] - cj));
    const xr = at(-1), fr = f(xr);
    if (fr < F[0]) { const xe = at(-2), fe = f(xe); if (fe < fr) { S[n] = xe; F[n] = fe; } else { S[n] = xr; F[n] = fr; } }
    else if (fr < F[n - 1]) { S[n] = xr; F[n] = fr; }
    else {
      const xc = at(0.5), fc = f(xc);
      if (fc < F[n]) { S[n] = xc; F[n] = fc; }
      else for (let i = 1; i <= n; i++) { S[i] = S[i].map((v, j) => S[0][j] + 0.5 * (v - S[0][j])); F[i] = f(S[i]); }
    }
  }
  const b = F.indexOf(Math.min(...F));
  return { x: S[b], f: F[b] };
}
// Timecourse fit. points: [{len, t, mm}] -- t = POWER-ON minutes, mm = well to band centre, one row per band per photo.
// Fits in two stages so the extra parameters have to earn their place: stage 1 is the strict linear-in-time model,
// stage 2 adds drift (p) and a well/measurement offset (d0). Huber loss so one mis-measured band cannot drag the fit.
function fitTimecourse(points, vcm, opts) {
  opts = opts || {};
  const tref = opts.tref || 60, hub = opts.huber == null ? 1.5 : opts.huber, init = modelFor(opts.pct || 1);
  const resid = m => points.map(pt => distanceMm(pt.len, m, pt.vcm == null ? vcm : pt.vcm, pt.t) - pt.mm);
  const loss = r => { let s = 0; for (const e of r) { const a = Math.abs(e); s += a <= hub ? 0.5 * e * e : hub * (a - 0.5 * hub); } return s; };
  const rmse = m => { const r = resid(m); let s = 0; for (const e of r) s += e * e; return Math.sqrt(s / r.length); };

  const mk3 = x => ({ muS: Math.exp(x[0]), muL: Math.exp(x[1]), gamma: Math.exp(x[2]), p: 1, d0: 0, tref, custom: true });
  const c3 = x => { const m = mk3(x); return m.muL >= m.muS ? 1e9 : loss(resid(m)); };
  const r3 = nelderMead(c3, [Math.log(init.muS), Math.log(init.muL), Math.log(init.gamma)], 900);
  const simple = mk3(r3.x);

  const mk5 = x => ({ muS: Math.exp(x[0]), muL: Math.exp(x[1]), gamma: Math.exp(x[2]), p: Math.exp(x[3]), d0: x[4], tref, custom: true });
  const c5 = x => { const m = mk5(x); return (m.muL >= m.muS || m.p < 0.3 || m.p > 2 || Math.abs(m.d0) > 15) ? 1e9 : loss(resid(m)); };
  let bestR = null;
  for (const seed of [r3.x.concat([0, 0]), r3.x.concat([Math.log(0.85), 0])]) {
    const r = nelderMead(c5, seed, 1500);
    if (!bestR || r.f < bestR.f) bestR = r;
  }
  const full = mk5(bestR.x);
  simple.rmse = rmse(simple); full.rmse = rmse(full);
  simple.n = full.n = points.length;
  // per-timepoint mean residual: a systematic trend here is the signature of drift the linear model cannot hold
  const byTime = {};
  resid(simple).forEach((e, i) => { const t = points[i].t; (byTime[t] = byTime[t] || []).push(e); });
  const drift = Object.keys(byTime).map(Number).sort((a, b) => a - b)
    .map(t => ({ t, meanResid: byTime[t].reduce((a, b) => a + b, 0) / byTime[t].length, n: byTime[t].length }));
  return { simple, full, drift, recommended: full.rmse < simple.rmse * 0.9 ? full : simple };
}

// Detection floor from a digest of known total load: points [{mass (ng), vis (0 = invisible ... 1 = strong)}].
// The threshold is only identifiable if the bands actually straddle it -- all-bright bands fit any floor below
// the faintest one, at zero residual, so the flags below matter more than the rmse.
function fitVisibility(points) {
  const lo = points.filter(p => p.vis <= 0.02).length;
  const mid = points.filter(p => p.vis > 0.02 && p.vis < 0.98).length;
  const hi = points.filter(p => p.vis >= 0.98).length;
  const cost = x => {
    const a = Math.exp(x[0]), b = a + Math.exp(x[1]);
    let s = 0;
    for (const p of points) s += (Math.max(0, Math.min(1, (p.mass - a) / (b - a))) - p.vis) ** 2;
    return s;
  };
  const r = nelderMead(cost, [Math.log(2), Math.log(18)], 800);
  return {
    visMin: Math.exp(r.x[0]), visFull: Math.exp(r.x[0]) + Math.exp(r.x[1]),
    rmse: Math.sqrt(r.f / points.length), span: { lo, mid, hi },
    visMinDetermined: lo >= 1 && mid >= 1, visFullDetermined: hi >= 1 && mid >= 1
  };
}

// How many rows the grid can honestly carry: adjacent rows must be separated by more than the combined
// positional error and band thickness, or two "pixels" are the same smear.
function maxRows(gel, rmseMm, bandMm) {
  const g = resolveGel(gel), pitch = 2 * Math.sqrt(rmseMm * rmseMm + ((bandMm || 1.5) / 2) ** 2);
  return { rows: Math.max(1, Math.floor((g.len - g.top) / pitch)), pitch, usable: g.len - g.top, currentPitch: g.rowH };
}

// points: [{len, mm}] measured from the wells on a gel run at vcm for `minutes`.
function fitModel(points, vcm, minutes, pct) {
  const init = modelFor(pct || 1);
  const toM = x => ({ muS: Math.exp(x[0]), muL: Math.exp(x[1]), gamma: Math.exp(x[2]) });
  const cost = x => { const m = toM(x); if (m.muL >= m.muS) return 1e9; let s = 0; for (const p of points) s += (distanceMm(p.len, m, vcm, minutes) - p.mm) ** 2; return s; };
  const r = nelderMead(cost, [Math.log(init.muS), Math.log(init.muL), Math.log(init.gamma)], 600);
  return Object.assign(toM(r.x), { rmse: Math.sqrt(r.f / points.length), custom: true });
}

// ---------------------------------------------------------------- solver
// For every lane (pixel column) find the enzyme set whose digest best paints the target rows.
async function solve(opts) {
  const { plasmids, target, maxEnz = 3, budget = 250000, topB = 8, onProgress, yieldFn } = opts;
  const gel = resolveGel(opts.gel);
  const times = opts.times && opts.times.length ? opts.times : [gel.time];
  const rows = gel.rows, maxCuts = opts.maxCutsPerEnzyme || 30;
  const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());

  // unique non-empty target columns
  const keyToU = new Map(), colsU = [], weightU = [], laneToU = [];
  target.forEach(col => {
    if (!col.some(v => v)) { laneToU.push(-1); return; }
    const key = col.map(v => v ? 1 : 0).join('');
    if (!keyToU.has(key)) { keyToU.set(key, colsU.length); colsU.push(Uint8Array.from(col, v => v ? 1 : 0)); weightU.push(0); }
    const u = keyToU.get(key); weightU[u]++; laneToU.push(u);
  });
  const nU = colsU.length;
  const best = times.map(() => colsU.map(() => []));          // best[ti][u] = sorted [{err,p,ids,key}]
  const insert = (list, err, p, ids) => {
    if (list.length === topB && err >= list[topB - 1].err) return;
    const key = p + ':' + ids.slice().sort((a, b) => a - b).join(',');
    for (const e of list) if (e.key === key) return;
    let i = list.length; while (i > 0 && list[i - 1].err > err) i--;
    list.splice(i, 0, { err, p, ids: ids.slice(), key });
    if (list.length > topB) list.pop();
  };

  // per-plasmid cutters grouped by incubation temperature
  const minAct = opts.minAct || 0;
  const pools = plasmids.map(pl => {
    const n = pl.seq.length, cutters = [];
    for (const e of opts.enzymes) {
      const cuts = findCuts(pl.seq, e, pl.circular);
      if (cuts.length && cuts.length <= maxCuts) cutters.push({ e, cuts, mask: minAct ? bufMask(e, minAct) : 15 });
    }
    const groups = {};
    cutters.forEach((c, i) => (groups[c.e.temp] = groups[c.e.temp] || []).push(i));
    return { pl, n, cutters, groups: Object.values(groups), muCache: new Float32Array(n + 1) };
  });

  const row = new Float32Array(rows), mu = new Float32Array(4096), vis = new Float32Array(4096);
  function evalCombo(p, ids, tis) {
    const pool = pools[p], n = pool.n;
    // A double digest has to run in ONE buffer. Reject any combination with no buffer where every
    // enzyme clears the activity threshold, before doing any of the expensive work.
    if (minAct && ids.length > 1) {
      let m = 15;
      for (const i of ids) m &= pool.cutters[i].mask;
      if (!m) return;
    }
    const cuts = ids.length === 1 ? pool.cutters[ids[0]].cuts : mergeCuts(ids.map(i => pool.cutters[i].cuts));
    const frags = fragmentsFromCuts(cuts, n, pool.pl.circular), nf = Math.min(frags.length, 4096);
    for (let f = 0; f < nf; f++) {
      const L = frags[f];
      mu[f] = pool.muCache[L] || (pool.muCache[L] = mobility(L, gel.model));
      vis[f] = visibility(L / n * gel.loadNg, gel);
    }
    const pen = 0.02 * ids.length;
    for (const ti of tis) {
      const Et = travel(gel.model, gel.vcm, times[ti]), D0 = gel.model.d0 || 0;
      row.fill(0); let stray = 0;
      for (let f = 0; f < nf; f++) {
        if (vis[f] <= 0) continue;
        const d = D0 + mu[f] * Et;
        if (d > gel.len) continue;
        if (d < gel.top) { stray += vis[f]; continue; }
        splat(row, vis[f], d, gel);
      }
      for (let u = 0; u < nU; u++) {
        const T = colsU[u]; let err = 0.5 * stray + pen;
        for (let r = 0; r < rows; r++) { const v = row[r] > 1 ? 1 : row[r]; err += T[r] ? 1 - v : v; }
        insert(best[ti][u], err, p, ids);
      }
    }
  }

  // Same arithmetic as evalCombo for a single run time, but writing the per-column errors out instead of
  // inserting them. Kept separate so the hot search loop does not pay for the extra indirection.
  function errorsInto(p, ids, ti, out) {
    const pool = pools[p], n = pool.n;
    const cuts = ids.length === 1 ? pool.cutters[ids[0]].cuts : mergeCuts(ids.map(i => pool.cutters[i].cuts));
    const frags = fragmentsFromCuts(cuts, n, pool.pl.circular), nf = Math.min(frags.length, 4096);
    for (let f = 0; f < nf; f++) {
      const L = frags[f];
      mu[f] = pool.muCache[L] || (pool.muCache[L] = mobility(L, gel.model));
      vis[f] = visibility(L / n * gel.loadNg, gel);
    }
    const Et = travel(gel.model, gel.vcm, times[ti]), D0 = gel.model.d0 || 0, pen = 0.02 * ids.length;
    row.fill(0);
    let stray = 0;
    for (let f = 0; f < nf; f++) {
      if (vis[f] <= 0) continue;
      const d = D0 + mu[f] * Et;
      if (d > gel.len) continue;
      if (d < gel.top) { stray += vis[f]; continue; }
      splat(row, vis[f], d, gel);
    }
    for (let u = 0; u < nU; u++) {
      const T = colsU[u];
      let err = 0.5 * stray + pen;
      for (let r = 0; r < rows; r++) { const v = row[r] > 1 ? 1 : row[r]; err += T[r] ? 1 - v : v; }
      out[u] = err;
    }
  }

  const choose = (n, k) => { let c = 1; for (let i = 0; i < k; i++) c = c * (n - i) / (i + 1); return n < k ? 0 : c; };
  let last = now(), done = 0, totalWork = 1;
  const tick = async label => {
    if (yieldFn && now() - last > 30) { if (onProgress) onProgress(Math.min(0.99, done / totalWork), label); await yieldFn(); last = now(); }
  };
  async function exhaustive(k, tis, label) {
    for (let p = 0; p < pools.length; p++) for (const G of pools[p].groups) {
      if (G.length < k) continue;
      const idx = Array.from({ length: k }, (_, i) => i);
      for (;;) {
        evalCombo(p, idx.map(i => G[i]), tis);
        if ((++done & 255) === 0) await tick(label);
        let i = k - 1; while (i >= 0 && idx[i] === G.length - k + i) i--;
        if (i < 0) break;
        idx[i]++; for (let j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
      }
    }
  }
  async function beam(k, tis, label) {
    for (const ti of tis) for (let u = 0; u < nU; u++) {
      for (const parent of best[ti][u].filter(e => e.ids.length === k - 1)) {
        const pool = pools[parent.p], G = pool.groups.find(g => g.includes(parent.ids[0]));
        for (const i of G) if (!parent.ids.includes(i)) { evalCombo(parent.p, parent.ids.concat(i), [ti]); done++; }
        await tick(label);
      }
    }
  }
  const layerCount = k => pools.reduce((s, pool) => s + pool.groups.reduce((a, G) => a + choose(G.length, k), 0), 0);
  const totalErr = ti => laneToU.reduce((s, u) => s + (u < 0 ? 0 : best[ti][u].length ? best[ti][u][0].err : rows), 0);

  if (nU) {
    const allT = times.map((_, i) => i);
    totalWork = 0;
    for (let k = 1; k <= maxEnz; k++) totalWork += Math.min(layerCount(k), k <= 2 ? Infinity : budget);
    for (let k = 1; k <= Math.min(2, maxEnz); k++) await exhaustive(k, allT, k === 1 ? 'single digests' : 'double digests');
    // deeper layers only on the most promising run times
    const active = opts.deepAllTimes ? allT : allT.slice().sort((a, b) => totalErr(a) - totalErr(b)).slice(0, 3);
    for (let k = 3; k <= maxEnz; k++) {
      if (layerCount(k) <= budget) await exhaustive(k, active, k + '-enzyme digests');
      else await beam(k, active, k + '-enzyme digests (beam)');
    }
  }

  let bestTi = 0;
  times.forEach((_, ti) => { if (totalErr(ti) < totalErr(bestTi) - 1e-9) bestTi = ti; });
  const optionFor = (p, ids, err) => {
    const pool = pools[p], enz = ids.map(i => pool.cutters[i].e);
    // sort name and product together, or a two-enzyme digest can end up naming the wrong tubes
    const pairs = enz.map(e => [e.name, e.product]).sort((a, b) => a[0].localeCompare(b[0]));
    const names = pairs.map(x => x[0]);
    const sb = sharedBuffer(names, minAct);
    return {
      err, plasmid: p, enzymes: names, temp: enz[0].temp, products: pairs.map(x => x[1]),
      buffer: sb ? sb.buffer : null, bufWorst: sb ? sb.worst : null, bufStar: sb ? sb.star : false
    };
  };
  const lanesAt = ti => laneToU.map(u => u < 0 ? { empty: true, options: [] } : {
    empty: false, options: best[ti][u].map(o => optionFor(o.p, o.ids, o.err))
  });

  // ---- digest reuse
  // Identical columns already collapse to one unique column, so they always share a reaction. This pass
  // goes further: columns that are merely SIMILAR can share one too, trading a bounded fidelity loss for
  // fewer tubes. Weighted set cover over the candidates the search already surfaced.
  // tol  = how much worse than ITS OWN best a column may be made, in error units (relative).
  // minQ = the quality no column may be pushed below, as a fraction of its on-pixels (absolute).
  // Both are needed. tol alone is unbounded: a column whose best error is 0.3 would accept 2.3 at tol=2,
  // and since an EMPTY lane already scores exactly on[u], any error at or above that is worse than
  // loading buffer only - which is how the reported quality used to go negative.
  function reusePass(ti, tol, minQ) {
    minQ = minQ == null ? 0.6 : Math.max(0, Math.min(1, minQ));
    const cands = [], seen = {};
    for (let u = 0; u < nU; u++) for (const o of best[ti][u]) {
      const k = o.p + ':' + o.ids.join(',');
      if (!seen[k]) { seen[k] = 1; cands.push({ p: o.p, ids: o.ids.slice() }); }
    }
    if (!cands.length) return null;
    const E = cands.map(c => { const a = new Float64Array(nU); errorsInto(c.p, c.ids, ti, a); return a; });

    const floor = new Float64Array(nU), onU = new Float64Array(nU), ceil = new Float64Array(nU);
    for (let u = 0; u < nU; u++) {
      let m = Infinity;
      for (let i = 0; i < E.length; i++) if (E[i][u] < m) m = E[i][u];
      floor[u] = m;
      let s = 0; const T = colsU[u];
      for (let r = 0; r < rows; r++) s += T[r];
      onU[u] = s || 1;
      // never below minQ of this column's own pixels, and never worse than its unconstrained best
      ceil[u] = Math.max(floor[u], Math.min(floor[u] + tol, onU[u] * (1 - minQ)));
    }
    const fits = (i, u) => E[i][u] <= ceil[u];

    const assigned = new Int32Array(nU).fill(-1);
    for (;;) {
      let pick = -1, bestCover = 0, bestSum = Infinity;
      for (let i = 0; i < cands.length; i++) {
        let cover = 0, sum = 0;
        for (let u = 0; u < nU; u++) if (assigned[u] < 0 && fits(i, u)) { cover += weightU[u]; sum += E[i][u]; }
        if (cover > bestCover || (cover === bestCover && cover > 0 && sum < bestSum)) { bestCover = cover; bestSum = sum; pick = i; }
      }
      if (pick < 0 || bestCover === 0) break;
      for (let u = 0; u < nU; u++) if (assigned[u] < 0 && fits(pick, u)) assigned[u] = pick;
      let left = 0; for (let u = 0; u < nU; u++) if (assigned[u] < 0) left++;
      if (!left) break;
    }
    // a column the cover could not serve within its ceiling keeps its own best digest
    for (let u = 0; u < nU; u++) if (assigned[u] < 0) {
      let m = Infinity, mi = 0;
      for (let i = 0; i < E.length; i++) if (E[i][u] < m) { m = E[i][u]; mi = i; }
      assigned[u] = mi;
    }
    let err = 0, on = 0, worstQ = 1; const used = {};
    for (let u = 0; u < nU; u++) {
      const e = E[assigned[u]][u];
      err += e * weightU[u]; on += onU[u] * weightU[u];
      worstQ = Math.min(worstQ, 1 - e / onU[u]);
      used[assigned[u]] = 1;
    }
    return {
      digests: Object.keys(used).length, err, quality: on ? 1 - err / on : 0, worstColumnQuality: worstQ,
      tol, minQ,
      lanes: laneToU.map(u => u < 0 ? { empty: true, options: [] } : {
        empty: false, options: [optionFor(cands[assigned[u]].p, cands[assigned[u]].ids, E[assigned[u]][u])]
      })
    };
  }
  const laneErrAt = (ti, lane) => { const u = laneToU[lane]; return u < 0 ? 0 : best[ti][u].length ? best[ti][u][0].err : rows; };
  if (onProgress) onProgress(1, 'done');
  return {
    time: times[bestTi], lanes: lanesAt(bestTi), times, lanesAt, laneErrAt, combos: done,
    reuse: opts.reuseTol > 0 && nU > 1 ? reusePass(bestTi, opts.reuseTol, opts.reuseMinQuality) : null,
    reuseAt: (ti, tol, minQ) => reusePass(ti, tol, minQ),
    timeScan: times.map((t, ti) => ({ time: t, err: totalErr(ti) }))
  };
}

// ---------------------------------------------------------------- shape library + catalog
const BITMAPS = {
  'Smiley': ['..XXXXXX..', '.X......X.', 'X..X..X..X', 'X..X..X..X', 'X........X', 'X.X....X.X', 'X..XXXX..X', 'X........X', '.X......X.', '..XXXXXX..'],
  'Heart': ['..XX...XX..', '.X..X.X..X.', 'X....X....X', 'X.........X', 'X.........X', '.X.......X.', '..X.....X..', '...X...X...', '....X.X....', '.....X.....'],
  'Heart (filled)': ['.XX...XX.', 'XXXX.XXXX', 'XXXXXXXXX', 'XXXXXXXXX', '.XXXXXXX.', '..XXXXX..', '...XXX...', '....X....'],
  'Invader': ['..X.....X..', '...X...X...', '..XXXXXXX..', '.XX.XXX.XX.', 'XXXXXXXXXXX', 'X.XXXXXXX.X', 'X.X.....X.X', '...XX.XX...'],
  'Arrow up': ['....X....', '...X.X...', '..X...X..', '.X.....X.', 'XXXX.XXXX', '...X.X...', '...X.X...', '...XXX...'],
  'House': ['....X....', '...X.X...', '..X...X..', '.X.....X.', 'XXXXXXXXX', 'X.......X', 'X..XXX..X', 'X..X.X..X', 'XXXXXXXXX']
};
// parametric outlines: (W lanes, H rows) -> H strings of W chars
const PARAMETRIC = {
  'Square': (W, H) => grid(W, H, (x, y) => x === 0 || y === 0 || x === W - 1 || y === H - 1),
  'Triangle': (W, H) => grid(W, H, (x, y) => y === H - 1 || Math.round(Math.abs(x - (W - 1) / 2) * (H - 1) / ((W - 1) / 2)) === y),
  'Diamond': (W, H) => grid(W, H, (x, y) => { const d = Math.abs(x - (W - 1) / 2) / ((W - 1) / 2), r = Math.round((1 - d) * (H - 1) / 2); return y === Math.round((H - 1) / 2) - r || y === Math.round((H - 1) / 2) + r; }),
  'Circle': (W, H) => grid(W, H, (x, y) => { const u = (x - (W - 1) / 2) / (W / 2), s = Math.sqrt(Math.max(0, 1 - u * u)) * (H - 1) / 2, c = (H - 1) / 2; return y === Math.round(c - s) || y === Math.round(c + s); }),
  'Cross \u2715': (W, H) => grid(W, H, (x, y) => { const d = Math.round(x * (H - 1) / (W - 1)); return y === d || y === H - 1 - d; }),
  'Chevron': (W, H) => grid(W, H, (x, y) => Math.round(Math.abs(x - (W - 1) / 2) * (H - 1) / ((W - 1) / 2)) === y),
  'Diagonal': (W, H) => grid(W, H, (x, y) => Math.round(x * (H - 1) / (W - 1)) === y),
  'Double helix': (W, H) => grid(W, H, (x, y) => { const s = Math.sin(x * 2 * Math.PI / Math.min(W, 8)) * (H - 1) / 2, c = (H - 1) / 2; return y === Math.round(c + s) || y === Math.round(c - s); }),
  'Wave': (W, H) => grid(W, H, (x, y) => y === Math.round((H - 1) / 2 + Math.sin(x * 2 * Math.PI / Math.min(W, 8)) * (H - 1) / 2)),
  'Equals': (W, H) => grid(W, H, (x, y) => y === 0 || y === H - 1)
};
function grid(W, H, on) { const out = []; for (let y = 0; y < H; y++) { let s = ''; for (let x = 0; x < W; x++) s += on(x, y) ? 'X' : '.'; out.push(s); } return out; }

// Place a bitmap on a lanesxrows target (columns of 0/1), centred horizontally, top edge at row y0. Clips.
function placeBitmap(bmp, lanes, rows, y0) {
  const W = bmp[0].length, x0 = Math.floor((lanes - W) / 2);
  const t = Array.from({ length: lanes }, () => new Array(rows).fill(0));
  bmp.forEach((line, y) => { for (let x = 0; x < W; x++) if (line[x] === 'X') { const X = x + x0, Y = y + y0; if (X >= 0 && X < lanes && Y >= 0 && Y < rows) t[X][Y] = 1; } });
  return t;
}
function shapeVariants(lanes, rows, opts) {
  opts = opts || {};
  const oStride = opts.offsetStride || 1, hStride = opts.heightStride || 2;
  const widths = opts.widths || [lanes, Math.min(lanes, 9), Math.min(lanes, 7)];
  const out = [];
  const add = (name, bmp, size) => { const H = bmp.length; if (bmp[0].length > lanes || H > rows) return; for (let y0 = 0; y0 + H <= rows; y0 += oStride) out.push({ name, size, y0, bmp }); };
  for (const [name, bmp] of Object.entries(BITMAPS)) add(name, bmp, bmp[0].length + '\u00d7' + bmp.length);
  for (const [name, fn] of Object.entries(PARAMETRIC)) {
    const seen = new Set();
    for (const W of widths) for (let H = 5; H <= rows; H += hStride) {
      if (seen.has(W + 'x' + H)) continue; seen.add(W + 'x' + H);
      add(name, fn(W, H), W + '\u00d7' + H);
    }
  }
  return out;
}

// Score every shape variant against the DNA in one solver pass; keep each shape's best variant + run time.
async function catalog(opts) {
  const gel = resolveGel(opts.gel), lanes = opts.lanes, variants = shapeVariants(lanes, gel.rows, opts.variantOpts);
  const wide = [];
  variants.forEach(v => { v.first = wide.length; v.target = placeBitmap(v.bmp, lanes, gel.rows, v.y0); v.target.forEach(c => wide.push(c)); });
  const res = await solve(Object.assign({}, opts, { target: wide, budget: opts.budget || 30000, deepAllTimes: true }));
  const byShape = new Map();
  for (const v of variants) {
    const on = v.target.reduce((s, c) => s + c.reduce((a, b) => a + b, 0), 0);
    res.times.forEach((t, ti) => {
      let err = 0; for (let l = 0; l < lanes; l++) err += res.laneErrAt(ti, v.first + l);
      const quality = Math.max(0, 1 - err / on);
      const cur = byShape.get(v.name + ' ' + v.size);
      if (!cur || quality > cur.quality) byShape.set(v.name + ' ' + v.size, {
        name: v.name, size: v.size, y0: v.y0, time: t, quality, err, on, target: v.target,
        lanes: res.lanesAt(ti).slice(v.first, v.first + lanes)
      });
    });
  }
  // best size per shape first, then the other sizes
  const all = Array.from(byShape.values()).sort((a, b) => b.quality - a.quality), seen = new Set();
  all.forEach(e => { e.bestOfShape = !seen.has(e.name); seen.add(e.name); });
  return { entries: all, variants: variants.length, combos: res.combos };
}

// ---------------------------------------------------------------- ladders (bp; bold = reference bands)
const LADDERS = {
  'NEB 1 kb Plus': { sizes: [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], bold: [3000, 1000, 500] },
  'NEB 1 kb': { sizes: [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500], bold: [3000] },
  'NEB 100 bp': { sizes: [1517, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], bold: [1000, 500] },
  'GeneRuler 1 kb Plus': { sizes: [20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300, 200, 75], bold: [5000, 1500, 500] }
};

root.GelArt = {
  ENZYMES, ENZYME_BY_NAME, GROUP_LABELS, LADDERS, AGAROSE_MODELS, PANELS, LAB_PANEL, panelEnzymes, parsePanel,
  BUFFERS, bufMask, sharedBuffer,
  revcomp, cleanSeq, parseSequences, parseSnapGene, findCuts, fragmentsFromCuts, mergeCuts, digest, cutCounts,
  modelFor, mobility, travel, distanceMm, fieldVcm, resolveGel, laneBands, rowIntensity,
  fitModel, fitTimecourse, fitVisibility, maxRows, solve,
  BITMAPS, PARAMETRIC, placeBitmap, shapeVariants, catalog
};
})(typeof globalThis !== 'undefined' ? globalThis : this);
