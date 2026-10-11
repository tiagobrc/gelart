---
Title:  Gel migration calibration protocol
Author: Tiago Castro
Lab: Victora
Project_driver: Tc
Experiment_Type: Calibration
Created: Monday 21st September 2026
Aliases:
  - Standing protocol for fitting DNA migration and detection floor to a rig, buffer and agarose percentage
Tags:
  - protocol
  - calibration
  - electrophoresis
  - gelart
Status: Active
Success: 1
---

> **A gel predicts nothing until it has been measured.** Fragment sizes are exact — they come from the
> sequence. Everything between a fragment size and a position on the gel runs through a model whose
> constants are literature values for somebody else's tank. This protocol replaces them with measured ones
> for one rig, one buffer and one agarose percentage, and in the same run measures what the imager can see.

Companion to [[Regular DNA Electrophoresis Agarose Gel Protocol]], which covers pouring and running.
This document covers only what is different about a run whose purpose is measurement.

---

## The shape

A calibration run produces four artefacts, named on the `Exp<MMDDYY>` date code:

| Artefact | What it is | Lives in |
| -------- | ---------- | -------- |
| `Exp<MMDDYY>.md` | notebook record — design, reasoning, results, corrections | vault, `Experiments/` |
| `Exp<MMDDYY>_protocol.html` | full run sheet, with the reasoning behind every number | vault, embedded in the note |
| `Exp<MMDDYY>_quicksheet.html` | two pages: page 1 pours and loads, page 2 photographs and measures | printed, taped to the bench |
| `Exp<MMDDYY>_model.json` | the fitted constants, and the only output anything downstream consumes | run folder + the tool that uses it |

A calibration without the fourth file has not finished. A `model.json` without the note that produced it
cannot be trusted a year later.

---

## What is being fitted

```
1/mu(L) = 1/mu_L - (1/mu_L - 1/mu_S) * exp(-L/gamma)
d       = d0 + mu(L) * E * t_ref * (t/t_ref)^p
```

| Parameter | Controls | What constrains it |
| --------- | -------- | ------------------ |
| `mu_S` | speed of small fragments | the small ladder bands, any one photo |
| `mu_L` | speed of large fragments | the large ladder bands, any one photo |
| `gamma` | where the curve bends between them | the full ladder, any one photo |
| `p` | whether distance stays linear in time | **the timepoints, and nothing else** |
| `d0` | constant offset in the measuring convention | all photos pooled |
| `visMin`, `visFull` | the detection floor and saturation point | the dilution series |

`p` is why this is not a one-photograph experiment. Every other parameter is pinned by a single endpoint
image, which is what an ordinary sizing gel gives you. If the buffer exhausts or the gel warms, migration
slows, and a model fitted at one endpoint is wrong everywhere except that endpoint, with nothing on the gel
to reveal it.

---

## The rules

These are the ones that stop the recurring errors. Everything else is ordinary electrophoresis.

### Field comes from the electrode gap, never the gel length

`volts = V/cm x gap`. Measure the gap platinum-to-platinum and **record it in the note**. A mini tank runs
~15 cm of buffer over a 7 cm gel; computing V/cm from gel length overstates the field about twofold. A run
sheet that says "147 V" without the gap it was computed for is not reusable on another box.

Stay at or below **5 V/cm**, per [[Genotyping protocol standard]]. Heat closes small gaps and a calibration
run is long.

### One agarose percentage per run

Mobility constants are percentage-specific. A model fitted at 1.2 % does not transfer to 1.5 %. Calibrate
the percentage you actually intend to use, all the way to an accepted `model.json`, before starting a
second. Do not pour two percentages and read them at one stop point.

### Five timepoints on one gel, photographed serially

Not a single endpoint, and **not staggered loading**. Staggered loading gives five timepoints on one image
for far less handling and is the obvious shortcut, but it answers a different question: a lane loaded at
t = 20 min spends its whole life in an already-warm gel, so it does not represent a run of that duration
*starting from zero*, which is what any prediction has to be. Re-photographing one gel follows the same
molecules cumulatively.

Timepoints at 0.2, 0.4, 0.6, 0.8 and 1.0 of `t_max`.

### Power-on minutes, never wall clock

Pause the stopwatch whenever the supply is off. Five gel removals cost 5–10 minutes of field time; counted
as run time they are indistinguishable from the drift being measured, and the fit will quietly absorb them
into `p`.

### A ruler in every frame, lying in the plane of the gel

Each photograph is scaled independently and the results pooled, so all frames must resolve to the same
millimetre scale. Gels swell and the camera will not return to the same height. Gel dimensions are not a
substitute. Without a ruler the frames cannot be combined and the run yields nothing.

### Photo 1 is a measurement, not a photograph

Measure it, refit the time scale, and rewrite the remaining timepoints before the gel runs on. An
uncalibrated prior can be wrong by a factor of two — it was, here, by comparison with [[Gen091226]] — and a
schedule built on it will either run fragments off the end or stop with everything bunched at the top.

If the observed/predicted ratio **varies across fragment sizes**, the model's shape is wrong and not merely
its speed. A single rescale will not fix that; take all five photos regardless.

### Hold out a real digest

At least one plasmid digest, loaded and measured, is excluded from the fit. Fit residuals always flatter
themselves; the held-out RMSE is the only honest error estimate and the only number acceptance is judged on.

Choose it by enumeration, not habit: **4–8 fragments, every band above the detection floor, every band
still on the gel at the last timepoint, and no two bands within 5 mm** at the final timepoint. A second
digest on a different plasmid with a different enzyme is the independent check.

### Sizes come from the maps, never from older notes

Inherited from [[Genotyping protocol standard]] and it applies unchanged. Every expected fragment is
computed in silico from the curated map. Record any correction in the note **with both values** — a silent
fix is indistinguishable from a typo next year.

### The dilution series has to go dark

A complete digest is equimolar, so band mass is exactly proportional to fragment length and no per-lot
table is needed — the masses follow from the sequence. Run a two-fold series until bands **disappear**.

A series in which every band stays visible fits *any* floor below the faintest one at zero residual, and
reports a perfect fit that means nothing. The threshold is identifiable only if some bands fail to appear.
Add loading dye **after** the dilutions are made, or it dilutes along with them.

### Score intensity only on early frames

If stain is cast into the gel, it migrates. EtBr carries a positive charge and moves toward the **cathode**,
so by the last timepoint the lower gel is stain-depleted and faint bands there are an artefact of the stain
rather than of the DNA. Band *positions* in late frames remain valid; *intensities* do not. Score the
dilution series from photo 1 or 2.

### Judge the endpoint by the ladder, not the dye

Also inherited. Give the run a **first check point** rather than a fixed duration, and a stop criterion in
terms of where a named ladder band sits. The smallest fragment of interest running off the end is the
failure mode that costs the whole run; [[Gen090926]] records a 260 bp product lost that way.

### Three ladders, not one

Left, centre and right of the comb. One ladder fits the model; three also measure how much the gel smiles,
which decides whether the full comb is usable or only the central lanes.

---

## Choosing the parameters

Nothing below is a fixed value. These are the rules for deriving them for a given run.

**Agarose percentage** — from the fragment range that matters, using the existing house assignments:
0.8 % above 3 kb, 1.0–1.2 % for plasmid digests and wide ranges, 1.5 % for 300–900 bp, 2–3 % for small
products and tight pairs.

**Running room** — measure the gel below the comb and record it. Every derived timepoint scales with it.
1.2–1.5 % wants at least 70 mm; 100 mm is comfortable; a 3 % gel resolving a tight pair needs 80 mm minimum.

**`t_max`** — the time at which the **smallest fragment you care about** reaches the end of the imageable
gel, less a 5 mm margin:

```
t_max = usable_mm / ( mu(L_small) * E )
```

computed from the current prior, then corrected after photo 1. Everything smaller than `L_small` is
expected to run off and is not measured.

**Load** — enough that every band of the validation digest clears the expected floor with margin. Work it
out from `mass_band = len_band / len_plasmid * load`, not by habit.

**Grid, if the model feeds a band-placement tool** — the row pitch you intend to use sets the acceptance
threshold, and the achievable pitch follows from the result:

```
max_rows = usable_mm / ( 2 * sqrt( rmse^2 + (band_width/2)^2 ) )
```

with `band_width` 2–3 mm per the house rule that "3.2 mm is about one band width — too tight to call".

---

## Measuring

- Scale **each photo separately** from its own ruler. Never carry a scale between frames.
- Prefer a lane profile (Fiji: rectangle down the lane, `Analyze > Plot Profile`) to eyeballing. It gives
  reproducible peak centres, and intensity for the dilution series at no extra cost.
- Measure **leading (anode) edge of the well to centre of band**, consistently. A constant offset is fitted
  out as `d0`; an inconsistent one is not.
- Leave a row blank for a band that ran off or cannot be read. Never estimate a position.

---

## Analysis and acceptance

Fit in two stages so the extra parameters must earn their place: stage 1 strictly linear in time, stage 2
adding drift and offset, with stage 2 adopted only if it is more than 10 % better. Use a robust loss, so a
single mis-measured band surfaces in an outlier list instead of dragging the model.

| Output | How to read it |
| ------ | -------------- |
| drift diagnostic | flat = migration linear in time; a monotonic walk early-negative to late-positive = it is not |
| **held-out RMSE** | **the number acceptance is judged on** |
| smiling | > 1.5 mm at the 90th percentile: use central lanes only, or drop the voltage |
| outliers | > 3 mm off; re-measure those bands before believing anything else |
| visibility | must report a value, not UNDETERMINED |

**Accept when the held-out RMSE is below half the row pitch you intend to use.** Below that a band lands
where it was predicted to. For a 20-row grid over 100 mm the pitch is 5 mm and the bar is 2.5 mm; on
simulated data with 0.25 mm measurement noise the fit recovers a known drifting gel to 0.29 mm, so that bar
is comfortable rather than strained.

**Reject and repeat** if visibility comes back UNDETERMINED (extend the dilution series downward), if
outliers survive re-measurement (suspect the digest — check for a partial on its own gel), or if smiling
exceeds the lane span you need.

---

## When to recalibrate

| Change | Recalibrate? |
| ------ | ------------ |
| Different agarose percentage | **Yes** — constants are percentage-specific |
| Different tank or electrode gap | **Yes** — and record the new gap |
| Different buffer system, or TAE made to a different recipe | **Yes** |
| New stain, or switching cast-in to post-stain | Positions no, **visibility yes** |
| New imager or camera setup | Positions no, **visibility yes** |
| Different V/cm on the same rig | No — `E` is an explicit term in the model |
| Different run time | No — that is the point of fitting `p` |
| New agarose lot, same supplier and percentage | No, but note it; if residuals grow, suspect it first |

---

## Common failures

| Symptom | Cause |
| ------- | ----- |
| Held-out RMSE good, predictions still wrong in use | calibrated at one percentage, applied at another |
| Drift term large and implausible | wall-clock time recorded instead of power-on minutes |
| Frames disagree with each other | scale carried between photos instead of read from each ruler |
| Visibility UNDETERMINED | dilution series never went dark |
| Faint bands only at the bottom of late frames | stain depletion, not DNA — rescore from photo 1 |
| Every band shifted by a constant | measuring convention changed mid-run; `d0` absorbs one offset, not two |
| Smallest fragments missing entirely | ran past `t_max`; the stop criterion was the dye, not the ladder |

---

Worked example: [[Exp092026]] — 1.2 % TAE, 4 V/cm on a 36.8 cm box, 20-well comb, `A_5 / SphI` and
`A_16 / BsrBI` held out, two-fold dilution series to 0.37 ng in the smallest band.

Built to the pattern of [[Genotyping protocol standard]] — three artefacts, V/cm not volts, a per-gel check
point, and sizes from maps rather than notes.
