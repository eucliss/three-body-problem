import json, sys, time
import numpy as np
sys.path.insert(0, "/tmp/claude-0/-home-user-elementFi-takehome/abe89a27-8b28-5869-a9f7-5ac0ed0d524d/scratchpad/verify/topology")
from topo import *

D = json.load(open("/home/user/three-body-problem/data/periodic_orbits.json"))
PUB = {
    "figure_eight_SD": "abAB",
    "butterfly_I": "(ab)2(AB)2",
    "moth_I": "baBABabABA",
    "goggles": "(ab)2ABBA(ba)2BAAB",
    "yin_yang_Ia": "(ab)2ABAba(BAB)",
}
# extra published (Table a1) used only as a hold-out check after calibration
EXTRA = {
    "butterfly_II": "(ab)2(AB)2",
    "butterfly_III": "(ab)2(ABA)(ba)2(BAB)",
    "moth_II": "(abAB)2A(baBA)2B",
    "moth_III": "ABabaBABabaBAbabABAbab",
    "yin_yang_Ib": "(ab)2ABAba(BAB)",
    "dragonfly": "b2(ABabAB)a2(BAbaBA)",
    "yarn": "baBABabABA",  # satellite of moth I, same family/word per Li&Liao
    "bumblebee": "b2(ABab)2A2(baBA)2baB2(abAB)2a2(BAba)2BA",
    "butterfly_IV": "((ba)2(AB)2)6A((ba)2(BA)2)6B",
}
which = sys.argv[1] if len(sys.argv) > 1 else "pub"
tab = PUB if which == "pub" else EXTRA
ok = {c: 0 for c in CONVENTIONS}
for name, pw in tab.items():
    o = D[name]
    t0 = time.time()
    sol = integrate(1.0, o["p1"], o["p2"], o["T"])
    ts, cr = crossings(sol, o["T"])
    y0 = sol.sol(ts); y1 = sol.sol(ts + o["T"])
    clos = np.linalg.norm(y1 - y0)
    target = canon(expand(pw))
    mids = "".join(str(c[1]) for c in cr)
    print(f"{name}: T={o['T']:.4f} closure={clos:.2e} ncross={len(cr)} mids={mids} pub={expand(pw)} ({time.time()-t0:.1f}s)")
    mp = min_puncture_dist(sol, ts, o["T"])
    print("   min dist to punctures:", {k: round(v, 4) for k, v in mp.items()})
    for c in CONVENTIONS:
        wa = free_reduce_cyclic(word_A(cr, c))
        wb = free_reduce_cyclic(word_B(sol, ts, o["T"], c))
        match = canon(wa) == target
        ok[c] += match
        flag = "MATCH" if match else ""
        agree = "" if canon(wa) == canon(wb) else "  A/B DISAGREE: B=" + wb
        print(f"   conv pole={c[0]} a-mid={c[1]} low={c[2]:+d}: {wa} {flag}{agree}")
print("convention scores (out of %d):" % len(tab))
for c in CONVENTIONS:
    print("  ", c, ok[c])
