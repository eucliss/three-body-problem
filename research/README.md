# The Three-Body Ladder: new periodic orbits

A systematic search for periodic orbits of three gravitating bodies, two of mass 1 and one lighter body of mass m₃, found a family of orbits built like a ladder. On **rung k** the two heavy bodies wind around each other **k times**, then unwind the same k turns while the light body swings out to the other side. Rungs 1 and 2 were already known. **Rungs 3 to 11 are new.**

**Findings page:** [The Three-Body Ladder](https://claude.ai/artifact/8jFqA4o72NgABXAKKGE3kR). It animates every rung, spells each orbit's topological word as it happens, and shows the search maps, the continuation branches and the checks. You can also open [`web/dist/ladder.html`](web/dist/ladder.html) directly; it is one self-contained file.

## Result

| | |
|---|---|
| New periodic orbits | **23**, at m₃ = 0.1, 0.15, 0.2, 0.25, 0.3, 0.5 and 0.75 |
| New topologies (free-group words) | **9**: (Ba)ⁿ(bA)ⁿ for n = 4 … 12, i.e. rungs k = 3 … 11 |
| Independent check | every orbit closes in a separately written integrator (DOP853 at 10⁻¹³; also Radau IIA for the first 16) to between 10⁻¹¹ and 8×10⁻⁹ |
| Stability | all linearly unstable; a perturbation grows by ×10 to ×1,200 per period |
| Search | 3.5 million initial conditions on 800 × 800 grids at seven masses, T ≤ 100 |
| Compared against | the 2,044 orbits of this slice published by Li & Liao (2017) and Li, Jing & Liao (2018), plus the free-fall and Šuvakov–Dmitrašinović orbits, by word under every physical symmetry and by starting values |

All orbits start in the standard isosceles slice with zero momentum and zero angular momentum (G = 1):

$$
\mathbf r_1=(-1,0),\ \mathbf r_2=(1,0),\ \mathbf r_3=(0,0),\qquad \dot{\mathbf r}_1=\dot{\mathbf r}_2=(v_1,v_2),\ \ \dot{\mathbf r}_3=-\tfrac{2}{m_3}(v_1,v_2),\qquad T^*=T\,|E|^{3/2}.
$$

The word of rung k, read on Montgomery's shape sphere in the conventions of Šuvakov & Dmitrašinović and Li & Liao, is

$$
w_k=(\mathtt{Ba})^{k+1}(\mathtt{bA})^{k+1}\ \sim\ \mathtt{A}(\mathtt{Ba})^{k}\mathtt{B}\,\mathtt{a}(\mathtt{bA})^{k}\mathtt{b},\qquad |w_k|=4k+4 .
$$

The (Ba) letters are written while the heavy pair winds and the (bA) letters while it unwinds. Half a period in, the light body passes back between the heavy pair and the three bodies line up as at the start, mirrored.

## All orbits

`data/ladder_orbits.json` has the full-precision values, words and stability. "Known word, new mass" means the word is catalogued at another mass but no orbit at this mass is published.

| rung k | m₃ | status | v₁ | v₂ | T | T* | closure (independent) | growth / period |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.1 | known word, new mass | 0.0569054662 | 0.0451246978 | 16.80791928 | 7.602 | – | – |
| 3 | 0.1 | **new** | 0.0650481454 | 0.0448812695 | 25.41340058 | 10.903 | 2.7e-11 | ×39.2 |
| 4 | 0.1 | **new** | 0.0686904544 | 0.0450225438 | 33.97243289 | 14.174 | 4.3e-10 | ×52.7 |
| 1 | 0.15 | known word, new mass | 0.0842024422 | 0.0541169070 | 9.10518365 | 4.842 | – | – |
| 2 | 0.15 | known word, new mass | 0.0912964179 | 0.0643434018 | 16.90843445 | 8.278 | – | – |
| 3 | 0.15 | **new** | 0.1010486152 | 0.0649052017 | 25.49973769 | 11.652 | 5.2e-11 | ×31.5 |
| 4 | 0.15 | **new** | 0.1058502641 | 0.0653281241 | 34.07565637 | 14.983 | 6.5e-11 | ×44.7 |
| 1 | 0.2 | known word, new mass | 0.1162683612 | 0.0716332483 | 9.37427550 | 5.430 | – | – |
| 2 | 0.2 | known word, new mass | 0.1282304973 | 0.0824316659 | 17.32303712 | 8.961 | – | – |
| 3 | 0.2 | **new** | 0.1388849414 | 0.0838744500 | 26.00917234 | 12.405 | 2.7e-11 | ×24.4 |
| 4 | 0.2 | **new** | 0.1445104318 | 0.0846291870 | 34.71762302 | 15.794 | 4.8e-10 | ×36.4 |
| 1 | 0.25 | known word, new mass | 0.1501626452 | 0.0891386899 | 9.74811317 | 6.024 | 1.8e-11 | ×12.1 |
| 2 | 0.25 | known word, new mass | 0.1668063963 | 0.1000193371 | 18.01071441 | 9.647 | 1.5e-11 | ×22.9 |
| 3 | 0.25 | **new** | 0.1781047333 | 0.1021982939 | 26.92164705 | 13.159 | 2.5e-11 | ×41.8 |
| 3' | 0.25 | **new** | 0.1580012072 | 0.0307487614 | 17.45166832 | 11.718 | 1.3e-09 | ×1218.8 |
| 4 | 0.25 | **new** | 0.1843204672 | 0.1032724455 | 35.88336424 | 16.604 | 2.1e-10 | ×69.4 |
| 5 | 0.25 | **new** | 0.1883005557 | 0.1039298198 | 44.86446896 | 20.006 | 1.5e-10 | ×110.9 |
| 6 | 0.25 | **new** | 0.1910975128 | 0.1043789281 | 53.85681740 | 23.377 | 6.8e-10 | ×169.9 |
| 7 | 0.25 | **new** | 0.1931869270 | 0.1047077084 | 62.85680119 | 26.726 | 1.0e-09 | ×249.7 |
| 8 | 0.25 | **new** | 0.1948165338 | 0.1049601857 | 71.86241721 | 30.056 | 1.8e-09 | ×353.8 |
| 9 | 0.25 | **new** | 0.1961288651 | 0.1051609964 | 80.87240857 | 33.372 | 3.2e-09 | ×486.0 |
| 10 | 0.25 | **new** | 0.1972121324 | 0.1053250678 | 89.88592298 | 36.676 | 4.6e-09 | ×650.1 |
| 11 | 0.25 | **new** | 0.1981240641 | 0.1054620032 | 98.90235043 | 39.970 | 7.6e-09 | ×850.1 |
| 1 | 0.3 | known word, new mass | 0.1857873813 | 0.1067434041 | 10.24065247 | 6.625 | – | – |
| 2 | 0.3 | known word, new mass | 0.2064780342 | 0.1174974427 | 18.96154518 | 10.336 | – | – |
| 3 | 0.3 | **new** | 0.2182936703 | 0.1202263282 | 28.23683145 | 13.914 | 3.0e-11 | ×50.4 |
| 4 | 0.3 | **new** | 0.2249356648 | 0.1215729077 | 37.58076077 | 17.413 | 7.0e-11 | ×90.4 |
| 5 | 0.3 | **new** | 0.2292343433 | 0.1223933941 | 46.95168116 | 20.862 | 1.5e-10 | ×148.5 |
| 1 | 0.5 | catalogued | 0.3420307308 | 0.1809369236 | 13.91533395 | 9.063 | – | ×7.2 |
| 2 | 0.5 | catalogued | 0.3697718457 | 0.1910065395 | 25.94419529 | 13.095 | – | ×10.1 |
| 3 | 0.5 | **new** | 0.3822826778 | 0.1946168976 | 38.30089073 | 16.915 | 4.8e-10 | ×9.9 |
| 4 | 0.5 | **new** | 0.3894814993 | 0.1965007869 | 50.75665222 | 20.616 | 4.2e-10 | ×12.0 |
| 5 | 0.5 | **new** | 0.3942248152 | 0.1976716358 | 63.26075912 | 24.241 | 5.7e-10 | ×17.1 |
| 1 | 0.75 | catalogued | 0.5337490177 | 0.3041674607 | 25.04603824 | 12.114 | – | – |
| 2 | 0.75 | known word, new mass | 0.5593298498 | 0.3101186237 | 46.61954804 | 16.495 | – | – |
| 3 | 0.75 | **new** | 0.5708736577 | 0.3121779086 | 68.71629418 | 20.585 | 3.7e-10 | ×119.1 |
| 1 | 1 | catalogued | 0.6150407229 | 0.5226158544 | 37.32052360 | 15.048 | – | – |

3′ at m₃ = 0.25 is a second orbit with the rung-3 word: a different orbit (T* 11.72 against 13.16) and much more unstable. Half a period in, it crosses the slice again at a different scale (0.8905), so the same orbit also starts from (v₁, v₂) = (0.1586412867, 0.1003753072) with T = 20.7687610 (closure 8×10⁻¹¹). Count it once.

## What is new, and what is not

- **New:** the orbits at m₃ < 0.5, where we found no published search of this slice; rungs 3 to 5 at m₃ = 0.5 and rung 3 at m₃ = 0.75, which lie inside the window that Li, Jing & Liao (2018) searched but are not in their catalogue; and the words (Ba)ⁿ(bA)ⁿ for n ≥ 4, which appear in none of the catalogues we could check.
- **Known:** rung 1's word (Li & Liao's I.A i.c. 3 at equal masses, and I.A i.c. 2 at m₃ = 0.5 and 0.75); rung 2's word (I.A i.c. 3 at m₃ = 0.5, and at equal masses in the high-precision set of Hristov et al. 2022). Rung 0's word would be that of the figure-eight.
- **Not new:** periods that grow linearly with word length are the topological dependence of Kepler's third law (Dmitrašinović & Šuvakov 2015). Ladders of words exist elsewhere, for example (BA)ⁿ(ba)ⁿ.
- **Not checked:** the Sofia group's equal-mass databases (Hristov, Hristova and coauthors, 2022–2025) could not be opened from this environment. If the higher rungs continue to m₃ = 1, some of these words may be listed there. Our continuation of rungs 2 and 3 stops at m₃ ≈ 0.85 and 0.76; rung 1 continues to m₃ = 1.
- **Numerical, not proven.** Existence rests on double-precision shooting confirmed by an independent integrator, not on interval arithmetic or high-precision (CNS) computation.

## Checks

1. **Pipeline calibration.** At m₃ = 1 the search recovers known orbits and flags none as new. At m₃ = 0.5, a mass Li, Jing & Liao (2018) catalogued, the full 800 × 800 scan refines 8 distinct orbits and all 8 are in their catalogue, so the comparison raises no false alarms there. That scan also misses the ladder at m₃ = 0.5, including the catalogued rungs 1 and 2. The m₃ = 0.5 rungs come from continuation. Grid scans miss orbits, which is also how a survey can miss rungs 3 to 5. The topology code reproduces the published words of 12 Šuvakov–Dmitrašinović orbits, including the 44-letter bumblebee.
2. **End-to-end anchors.** Following rung 1 by continuation from m₃ = 0.25 lands on Li, Jing & Liao's I.A i.c. 2 at m₃ = 0.75 and Li & Liao's equal-mass I.A i.c. 3, both to all ten published digits (`results/anchors.json`). Continuation of rungs 1 and 2 downward lands on orbits the grid found independently at m₃ = 0.15, 0.2 and 0.3.
3. **Independent review** (`verify/`, `results/verification.json`). Separate reviewers rewrote the integrator (SciPy DOP853 and Radau IIA, unregularized) and the topology code from scratch, and compared words against all catalogue words under rotation, inversion, mirror image and relabelling of the equal bodies. Every claim held: the stated period is minimal, there is no collision (closest approach 0.025), and energy and angular momentum are conserved to 10⁻¹². A skeptic pass corrected two over-claims by the reviewers: a "generous" word matching that allowed non-physical case swaps, and a symmetry-class statement about 3′.
4. **Stability** from the full 12 × 12 monodromy matrix, after removing the eight trivial multipliers exactly (restricting to perturbations that conserve the integrals, then quotienting out time shift and rotation). The method marks the figure-eight and moth I stable and butterfly I, yin-yang and goggles unstable, as published.

## How the search works

| Step | Script |
|---|---|
| Grid scan: integrate each (v₁, v₂) on an 800 × 800 grid to T = 100, record local minima of the phase-space return distance below 0.15 | `search.js` |
| Cluster the near-returns per small cell (keep the best 3), cap at 3,000 | `prefilter.js` |
| Levenberg–Marquardt shooting on (v₁, v₂, T) to a closure of 10⁻¹¹ (exact-time landing through the regularized step) | `refine.js`, `process.js` |
| Word on the shape sphere; canonical form up to rotation and inversion; comparison with the catalogues | `topology.js`, `catalog.js`, `calibrate_topology.js` |
| Continuation in m₃, parametrized by the period (passes folds in m₃) | `contT.js`, `atmass.js`, `cont_m3.js` |
| Next-rung prediction from the even period spacing, with local box scans | `ladder.js` |
| Monodromy matrix and Floquet multipliers | `monodromy.js`, `stability.py` |
| Return-distance maps for the page | `map.js` |
| Page data and page | `web/build_data.js`, `web/build_page.js` → `web/dist/ladder.html` |

The integrator is the repository's regularized Bulirsch–Stoer LogH solver (`../web/solver.js`), bit-identical across JavaScript engines.

Reproduce the main search at one mass (Node 18+, 4 cores, about 15 minutes):

```sh
cd research
./run_masses.sh 0.25            # search, cluster, refine, classify → runs/m025_N800_T100_orbits.jsonl
node contT.js k3 0.25 0.178104733251427 0.10219829385020422 26.92164705317384 0.05 70 runs/contT_k3.json
node web/build_data.js && node web/build_page.js
```

Large raw scan outputs (`runs/*_N800_T100.jsonl`, the maps) are not committed; the scripts regenerate them.

## References

- M. Šuvakov & V. Dmitrašinović, *Phys. Rev. Lett.* 110, 114301 (2013).
- X. Li & S. Liao, *Sci. China Phys. Mech. Astron.* 60, 129511 (2017).
- X. Li, Y. Jing & S. Liao, *Publ. Astron. Soc. Japan* 70, 64 (2018).
- X. Li & S. Liao, *New Astron.* 70, 22 (2019).
- V. Dmitrašinović & M. Šuvakov, *Phys. Lett. A* 379, 1939 (2015).
- V. Dmitrašinović, A. Hudomal, M. Shibayama & A. Sugita, *J. Phys. A* 51, 315101 (2018).
- I. Hristov et al., New families of periodic orbits for the planar three-body problem computed with high precision (2022), and later equal-mass searches (2023–2026).
- R. Montgomery, *Nonlinearity* 11, 363 (1998).
