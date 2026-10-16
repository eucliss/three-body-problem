# Run outputs

What each committed file in `runs/` is, the command that made it, and what it supports. Raw scan
outputs (`*_N800_T100.jsonl`, `*_clustered.jsonl`, `ladder_box_*.jsonl`, `test_m1.jsonl`,
`pilot_m025.jsonl`) are ignored; the commands below regenerate them.

## Grid searches (one per mass)

| Files | Command | Notes |
|---|---|---|
| `m01_N800_T100.log`, `m01_N800_T100_orbits.jsonl` | `node search.js 0.1 800 100 runs/m01_N800_T100.jsonl 0.25` then `node process.js runs/m01_N800_T100.jsonl runs/m01_N800_T100_orbits.jsonl` | threshold 0.25, no clustering, process.js defaults (dmax 0.25, 2 hits per point) |
| `m025_N800_T100.log`, `m025_N800_T100_orbits.jsonl` | same with 0.25 in place of 0.1 | as above |
| `m015_*`, `m02_*`, `m03_*`, `m04_*`, `m05_*` (`.log`, `_orbits.jsonl`) | `./run_masses.sh 0.15 0.2 0.3 0.4 0.5` | threshold 0.15, `prefilter.js` clustering (best 3 per cell, cap 3,000), `PATIENCE=16` |

The two settings were not mixed by accident: m₃ = 0.1 and 0.25 were run first, before the clustering
step existed, and the looser threshold is what found rungs 5–7 at m₃ = 0.25 and rung 3 at m₃ = 0.1
directly in the grid. Rerunning those two masses with `run_masses.sh` finds fewer rungs.

## Calibration and pilot

- `test_m1_orbits.jsonl`: m₃ = 1 calibration, `node search.js 1 120 40 runs/test_m1.jsonl 0.25` then `node process.js runs/test_m1.jsonl runs/test_m1_orbits.jsonl`. 18 orbits: 10 with catalogued words, and 8 whose words are powers of shorter words (repeat traversals and satellites); no primitive word is flagged as new.
- `pilot_m025_orbits.jsonl`: first look at m₃ = 0.25 (`node search.js 0.25 300 60 runs/pilot_m025.jsonl 0.25`, then `process.js`); exploratory, superseded by the 800 × 800 run.

## Continuation branches (read by `web/build_data.js`)

`node contT.js <name> <m3> <v1> <v2> <T> <dT> <Tstop> <out>` follows a branch in m₃ with the period
as parameter; `dT` is the initial step (the step adapts).

| File | Seed | dT | Tstop |
|---|---|---|---|
| `contT_k1.json` | L1@0.25 | 0.05 | 60 |
| `contT_k2.json` | L2@0.25 | 0.05 | 80 |
| `contT_k3.json` | L3@0.25 | 0.5 | not recorded (stalled at T = 70.2) |
| `contT_k4.json`, `contT_k5.json` | L4@0.25, L5@0.25 | 0.5 | not recorded (stalled) |
| `contT_k6.json` | L6@0.25 | 0.5 (first step halved) | not recorded (stalled) |
| `contT_k3_hi.json`, `contT_k3_hi_rev.json` | L3@0.75, forward and backward | +0.05, −0.05 (first step halved) | not recorded |
| `contT_k{1,2,3,4}_down.json` | L1–L4@0.25 | −0.05 | 8.6, 16.5, 25.2, 33.8 |
| `contM_k2.json`, `contM_k3.json` | last points of contT_k2/k3 | `node cont_m3.js runs/contT_kN.json 0.004 1.05 runs/contM_kN.json` (natural continuation in m₃) | |

Rungs 7–11 were not continued away from m₃ = 0.25. `atmass_extra.jsonl` holds the exact rungs at
m₃ = 0.15, 0.2, 0.3 (`node atmass.js runs/contT_k1_down.json 0.15 0.2`, and likewise for k2_down,
k4_down, k1 and k2 at 0.3).

## Ladder extension and stability

- `ladder_m025.json`, `ladder_seeds.txt` (plain-text list of rungs 6–11 as found): `node ladder.js 0.25 12 runs/ladder_m025.json runs/m025_N800_T100_orbits.jsonl`. Rungs 8–11 came from its local box scans.
- `claims_monodromy.json`, `claims_extra_monodromy.json`: `node stability_all.js results/claims.json runs/claims_monodromy.json` (and `claims_extra`), then `python3 stability.py runs/claims_monodromy.json > results/claims_stability.json`.
- `ladder_m025_monodromy.json`, `k3_branch_monodromy.json`: monodromy matrices along the m₃ = 0.25 ladder and along the k = 3 branch, used to check the natural-continuation stops for a multiplier at +1 (there was none, so those stops were not folds).
- `map_m025_N400.json`, `map_m05_N400.json`: `node map.js 0.25 400 60 runs/map_m025_N400.json` and `node map.js 0.5 400 60 runs/map_m05_N400.json` (about a minute each); the page's search maps.

## Exploratory (not used by any result)

`arc_k3.json` (abandoned arclength continuation), `branch_k3_up.json`, `branch_k3b_up.json`,
`cont_k*_up.json`, `cont_k*_down.json` (natural continuation; its apparent branch ends turned out to
be artefacts, see the T-continuation instead), `k8_seed.json`, `ladder_k8_box.jsonl`,
`ladder_m05.jsonl` (rungs 4–5 at m₃ = 0.5 before they were refined into `results/claims.json`).
