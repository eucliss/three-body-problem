#!/bin/sh
# full grid searches (800 x 800, T <= 100), clustering of near-returns, refinement and classification
for m3 in "$@"; do
  tag=$(echo "$m3" | tr -d .)
  base="runs/m${tag}_N800_T100"
  [ -s "$base.jsonl" ] || node search.js "$m3" 800 100 "$base.jsonl" 0.15 > "$base.log" 2>&1
  node prefilter.js "$base.jsonl" "${base}_clustered.jsonl" 0.15 0.004 0.5 3 3000 >> "$base.log" 2>&1
  PATIENCE=16 node process.js "${base}_clustered.jsonl" "${base}_orbits.jsonl" 0.15 1 >> "$base.log" 2>&1
done
