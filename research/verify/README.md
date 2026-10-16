# Independent checks, as run

These scripts were written and run by separate reviewers, who did not use the project's solver or
topology code. They are kept as a record of what was run, not as polished tools. Their results are
in `../results/verification_raw.json` (every reviewer's full output, including the skeptic's
judgement) and are summarised in `../results/verification.json`.

| Folder | Reviewer | Main scripts | Output |
|---|---|---|---|
| `numerical/` | periodicity and stability | `nb.py` (unregularized planar equations), `check.py` (SciPy DOP853 at 10⁻¹³, closure at three tolerances, minimal period, closest approach), `altmeth.py` (Radau IIA), `stab.py` (own variational equations), `shaperet.py` | `claims.log`, `radau.log`, `claims_stab.log`, `classic.log` (figure-eight, butterfly I, moth I, goggles, yin-yang Ia), `extra.log` (the 7 rungs added later), `L3_m05_recheck.log` |
| `topology/` | words | `topo.py` (own shape-sphere word code, two crossing-reading methods), `calibrate.py`, `catval.py`, `claims_run.py`, `catsearch.py`, `b4.py` (butterfly IV) | in `verification_raw.json` |
| `novelty/` | catalogue comparison | `parse_cat.py`, `compare_words.py`, `numeric_cmp.py`, `cmp3d.py` | in `verification_raw.json` |
| `skeptic/` | final judge | `words.py` (word matching under the physical symmetry group), `half.py` (the T/2 return to the slice) | in `verification_raw.json` |

`numerical/check.py` runs from anywhere: `python3 -I check.py ../../results/claims_extra.json out.json`.

The `topology/` and `skeptic/` scripts hard-code the reviewers' scratch directory and the
checkout path `/home/user/three-body-problem`, and several read the Li & Liao catalogue from a local
clone of https://github.com/sjtu-liao/three-body (commit b96b019). To rerun them elsewhere, edit those
paths at the top of each script.
