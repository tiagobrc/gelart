/* Gel Art rendering engine -- turns fragment lists into a gel-doc image on a millimetre grid.
 * Everything is laid out in MILLIMETRES and scaled to pixels at the end, so the drawn ruler
 * measures the predicted gel exactly the way a ruler in the photograph measures the real one.
 */
(function (root) {
'use strict';
const G = root.GelArt;

// EtBr under UV is orange-red; a gel doc clips it toward white where the band is thick.
const PALETTES = {
  uv:   [[0, [10, 9, 13]], [0.10, [36, 17, 21]], [0.32, [126, 48, 21]], [0.58, [219, 119, 47]], [0.84, [255, 199, 129]], [1, [255, 250, 238]]],
  gray: [[0, [9, 9, 10]], [0.15, [40, 40, 42]], [0.50, [140, 140, 142]], [0.85, [225, 225, 228]], [1, [255, 255, 255]]],
  blue: [[0, [6, 10, 14]], [0.20, [14, 50, 70]], [0.50, [30, 140, 170]], [0.80, [140, 225, 235]], [1, [240, 255, 255]]]
};

function ramp(stops, t) {
  t = t < 0 ? 0 : t > 1 ? 1 : t;
  for (let i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      const a = stops[i - 1], b = stops[i], u = (t - a[0]) / ((b[0] - a[0]) || 1);
      return [a[1][0] + u * (b[1][0] - a[1][0]), a[1][1] + u * (b[1][1] - a[1][1]), a[1][2] + u * (b[1][2] - a[1][2])];
    }
  }
  return stops[stops.length - 1][1];
}

// Band thickness: a freshly loaded band is sharp and broadens by diffusion as the run goes on.
// Calibrated to the house rule that a band reads as 2-3 mm wide (FWHM = 2.355 sigma).
function sigmaMm(o) { return (o.sigma0 || 0.75) * (1 + (o.time || 60) / 300); }

function layout(o) {
  const gutter = o.showRuler === false ? 3 : 13;
  const laneW = o.laneWidthMm || 5.5;
  const wellTop = 5;                       // gel above the wells
  const foot = o.labelBelow ? 9 : 4;
  return {
    gutter, laneW, wellTop, foot, nl: o.nLanes,
    widthMm: gutter + o.nLanes * laneW + 3,
    heightMm: wellTop + o.gelLen + foot
  };
}

function renderGel(canvas, opts) {
  const gel = G.resolveGel(opts.gel);
  const L = layout({ nLanes: opts.lanes.length, gelLen: gel.len, laneWidthMm: opts.laneWidthMm, showRuler: opts.showRuler, labelBelow: opts.labelBelow });
  const s = opts.pxPerMm || 6;
  const W = Math.round(L.widthMm * s), H = Math.round(L.heightMm * s);
  const dpr = opts.dpr || 1;
  // The photographic layer is built pixel by pixel, so it must be built at DEVICE resolution --
  // putImageData ignores the canvas transform and would otherwise fill only the top-left 1/dpr of the canvas.
  const DW = Math.round(W * dpr), DH = Math.round(H * dpr), sd = s * dpr;
  canvas.width = DW; canvas.height = DH;
  canvas.style.width = W + 'px'; canvas.style.height = H + 'px';
  const ctx = canvas.getContext('2d');
  ctx.setTransform(1, 0, 0, 1, 0, 0);

  const pal = PALETTES[opts.palette] || PALETTES.uv;
  const sig = sigmaMm({ sigma0: opts.sigma0, time: gel.time });

  // vertical intensity profile per lane, one sample per device pixel row
  const profiles = opts.lanes.map(lane => {
    const prof = new Float32Array(DH);
    if (!lane || !lane.bands) return prof;
    for (const b of lane.bands) {
      if (!(b.vis > 0) || b.dist < 0 || b.dist > gel.len) continue;
      const yc = (L.wellTop + b.dist) * sd, r = Math.ceil(4 * sig * sd);
      const y0 = Math.max(0, Math.floor(yc - r)), y1 = Math.min(DH, Math.ceil(yc + r));
      for (let y = y0; y < y1; y++) { const dz = (y - yc) / sd / sig; prof[y] += b.vis * Math.exp(-0.5 * dz * dz); }
    }
    return prof;
  });

  // deterministic speckle so the image is stable across redraws
  let seed = 20260919;
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };

  const img = ctx.createImageData(DW, DH), d = img.data;
  const gelX0 = L.gutter - 1.5, gelX1 = L.gutter + L.nl * L.laneW + 1.5;
  const gelY0 = 1.5, gelY1 = L.wellTop + gel.len + 2;
  const edge = 0.55;
  for (let y = 0; y < DH; y++) {
    const ymm = y / sd;
    const inGelY = ymm > gelY0 && ymm < gelY1;
    for (let x = 0; x < DW; x++) {
      const xmm = x / sd;
      let I = 0;
      const li = Math.floor((xmm - L.gutter) / L.laneW);
      if (li >= 0 && li < L.nl && ymm > L.wellTop - 1.2) {
        const within = (xmm - L.gutter) - li * L.laneW;
        let hx = 1;
        if (within < edge) hx = 0.5 - 0.5 * Math.cos(Math.PI * within / edge);
        else if (within > L.laneW - edge) hx = 0.5 - 0.5 * Math.cos(Math.PI * (L.laneW - within) / edge);
        I = profiles[li][y] * hx;
      }
      const bg = (inGelY && xmm > gelX0 && xmm < gelX1) ? 0.05 : 0;
      const v = Math.min(1.35, I) + bg + (rnd() - 0.5) * 0.024;
      const c = ramp(pal, v);
      const o = (y * DW + x) * 4;
      d[o] = c[0]; d[o + 1] = c[1]; d[o + 2] = c[2]; d[o + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);   // overlays from here on are drawn in mm-scaled CSS pixels

  const mm = v => v * s;
  const ink = opts.ink || 'rgba(255,238,215,0.78)';
  ctx.font = `500 ${Math.max(7, 3.0 * s)}px "JetBrains Mono", ui-monospace, monospace`;
  ctx.textBaseline = 'middle';

  // wells
  ctx.fillStyle = 'rgba(0,0,0,0.55)';
  ctx.strokeStyle = 'rgba(255,225,190,0.22)';
  ctx.lineWidth = Math.max(0.5, 0.18 * s);
  for (let i = 0; i < L.nl; i++) {
    const x = mm(L.gutter + i * L.laneW + L.laneW * 0.16), w = mm(L.laneW * 0.68);
    ctx.fillRect(x, mm(L.wellTop - 1.6), w, mm(1.6));
    ctx.strokeRect(x, mm(L.wellTop - 1.6), w, mm(1.6));
  }

  // art grid: the rows the solver is actually aiming at
  if (opts.showGrid) {
    ctx.strokeStyle = 'rgba(120,200,255,0.20)';
    ctx.lineWidth = Math.max(0.5, 0.12 * s);
    ctx.setLineDash([mm(0.8), mm(0.8)]);
    for (let r = 0; r <= gel.rows; r++) {
      const y = mm(L.wellTop + gel.top + r * gel.rowH);
      ctx.beginPath(); ctx.moveTo(mm(L.gutter), y); ctx.lineTo(mm(L.gutter + L.nl * L.laneW), y); ctx.stroke();
    }
    ctx.setLineDash([]);
  }

  // ruler -- origin at the wells, which is where migration is measured from
  if (opts.showRuler !== false) {
    ctx.strokeStyle = ink; ctx.fillStyle = ink;
    ctx.lineWidth = Math.max(0.6, 0.16 * s);
    const rx = mm(L.gutter - 2.2);
    ctx.beginPath(); ctx.moveTo(rx, mm(L.wellTop)); ctx.lineTo(rx, mm(L.wellTop + gel.len)); ctx.stroke();
    ctx.textAlign = 'right';
    for (let v = 0; v <= gel.len; v++) {
      const y = mm(L.wellTop + v);
      const major = v % 10 === 0, med = v % 5 === 0;
      ctx.globalAlpha = major ? 0.95 : med ? 0.6 : 0.35;
      ctx.beginPath(); ctx.moveTo(rx, y); ctx.lineTo(rx - mm(major ? 2.6 : med ? 1.8 : 1.0), y); ctx.stroke();
      if (major) { ctx.globalAlpha = 0.95; ctx.fillText(String(v / 10), rx - mm(3.2), y); }
    }
    ctx.globalAlpha = 1;
    ctx.save();
    ctx.translate(mm(1.6), mm(L.wellTop + gel.len / 2));
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.font = `500 ${Math.max(6, 2.5 * s)}px "Barlow Semi Condensed", system-ui, sans-serif`;
    ctx.fillText('cm from wells', 0, 0);
    ctx.restore();
  }

  // ladder band sizes
  ctx.textAlign = 'left';
  ctx.font = `500 ${Math.max(6, 2.4 * s)}px "JetBrains Mono", ui-monospace, monospace`;
  opts.lanes.forEach((lane, i) => {
    if (!lane || !lane.ladder || !lane.bands) return;
    const x = mm(L.gutter + i * L.laneW + L.laneW + 0.6);
    ctx.fillStyle = 'rgba(255,240,220,0.72)';
    // Where the ladder compresses at the top, neighbouring labels would overprint each other -- drop those.
    let lastY = -99;
    for (const b of lane.bands) {
      if (b.dist > gel.len || !b.label || b.dist - lastY < 2.7) continue;
      lastY = b.dist;
      ctx.fillText(b.label, x, mm(L.wellTop + b.dist));
    }
  });

  // lane numbers under the gel
  if (opts.labelBelow) {
    ctx.textAlign = 'center';
    ctx.font = `600 ${Math.max(6, 2.6 * s)}px "Barlow Semi Condensed", system-ui, sans-serif`;
    // Number the sample lanes 1..n so the gel matches the recipe table; the ladder is not a sample.
    let n = 0;
    for (let i = 0; i < L.nl; i++) {
      const lane = opts.lanes[i], isLadder = !!(lane && lane.ladder);
      if (!isLadder) n++;
      ctx.fillStyle = isLadder ? 'rgba(255,220,170,0.85)' : 'rgba(255,240,220,0.55)';
      ctx.fillText(isLadder ? 'L' : String(n), mm(L.gutter + i * L.laneW + L.laneW / 2), mm(L.wellTop + gel.len + 4.5));
    }
  }
  return { widthMm: L.widthMm, heightMm: L.heightMm, layout: L, sigma: sig };
}

function laneFromDigest(plasmid, enzymes, gel) {
  const frags = G.digest(plasmid, enzymes);
  return { bands: G.laneBands(frags, gel), fragments: frags, enzymes, plasmid: plasmid.name, label: enzymes.join(' + ') };
}

// Commercial ladders are formulated so the reference bands are brighter; band MASSES are not published
// per lot here, so brightness is rendered as a two-level convention rather than invented nanograms.
function ladderLane(name, gel) {
  const spec = G.LADDERS[name], g = G.resolveGel(gel);
  return {
    ladder: true, label: name,
    bands: spec.sizes.map(len => ({
      len, dist: G.distanceMm(len, g.model, g.vcm, g.time),
      vis: spec.bold.indexOf(len) >= 0 ? 1 : 0.5,
      label: len >= 1000 ? (len / 1000) + 'k' : String(len)
    })).filter(b => b.dist <= g.len)
  };
}

root.GelRender = { renderGel, laneFromDigest, ladderLane, PALETTES, sigmaMm, layout };
})(typeof globalThis !== 'undefined' ? globalThis : this);
