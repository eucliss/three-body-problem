import json, sys, time
import numpy as np
sys.path.insert(0, "/tmp/claude-0/-home-user-elementFi-takehome/abe89a27-8b28-5869-a9f7-5ac0ed0d524d/scratchpad/verify/topology")
from topo import *
C = json.load(open("/home/user/three-body-problem/research/results/claims.json"))
WIN = [((1, 2), 1, +1), ((1, 2), 2, -1)]
sel = sys.argv[1:] if len(sys.argv) > 1 else None
res = []
for c in C:
    if sel and c["id"] not in sel:
        continue
    t0 = time.time()
    m3, v1, v2, T = c["m3"], c["v1"], c["v2"], c["T"]
    E0 = energy(initial_state(m3, v1, v2), (1, 1, m3))
    out = {"id": c["id"], "claimed": c["word"], "status": c["status"], "E_chk": E0 - c["E"],
           "Lf_ok": len(c["word"]) == c["Lf"] and free_reduce_cyclic(c["word"]) == c["word"]}
    words = {}
    for tol in (1e-11, 1e-13):
        sol = integrate(m3, v1, v2, T, rtol=tol, atol=tol)
        ts, cr = crossings(sol, T)
        cs = crossing_states(sol, ts, T)
        clos = float(np.linalg.norm(sol.sol(ts + T) - sol.sol(ts)))
        for conv in WIN:
            wa = free_reduce_cyclic(word_A(cr, conv))
            wb = free_reduce_cyclic(word_B_events(cs, conv))
            words[(tol, conv)] = (wa, canon(wa) == canon(wb))
        out[f"closure_{tol:g}"] = clos
        out[f"ncross_{tol:g}"] = len(cr)
    mp = min_puncture_dist(sol, ts, T)
    out["min_punct"] = {str(k): round(v, 5) for k, v in mp.items()}
    out["min_slope"] = float(min(x[3] for x in cr))
    # sub-period check
    out["subperiod_closure"] = {n: float(np.linalg.norm(sol.sol(ts + T / n) - sol.sol(ts))) for n in (2, 3, 4)}
    tgt = canon(c["word"])
    out["words"] = {f"{k[0]:g}|{k[1][1]}{k[1][2]:+d}": [v[0], canon(v[0]) == tgt, v[1]] for k, v in words.items()}
    out["all_match"] = all(canon(v[0]) == tgt for v in words.values())
    out["AB_agree"] = all(v[1] for v in words.values())
    out["sec"] = round(time.time() - t0, 1)
    res.append(out)
    print(json.dumps(out), flush=True)
json.dump(res, open("/tmp/claude-0/-home-user-elementFi-takehome/abe89a27-8b28-5869-a9f7-5ac0ed0d524d/scratchpad/verify/topology/claims_result%s.json" % ("" if not sel else "_sel"), "w"), indent=1)
