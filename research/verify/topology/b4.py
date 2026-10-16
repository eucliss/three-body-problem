import sys, json
sys.path.insert(0, "/tmp/claude-0/-home-user-elementFi-takehome/abe89a27-8b28-5869-a9f7-5ac0ed0d524d/scratchpad/verify/topology")
from topo import *
D = json.load(open("/home/user/three-body-problem/data/periodic_orbits.json"))
o = D["butterfly_IV"]
pub = expand("((ba)2(AB)2)6A((ba)2(BA)2)6B")
def maps(w):
    import itertools
    out = {}
    for name, m in {"id":{'a':'a','b':'b'}, "a<->b":{'a':'b','b':'a'}, "a<->A":{'a':'A','b':'B'}, "a<->B":{'a':'B','b':'A'}}.items():
        mm = dict(m); mm.update({k.upper(): v.swapcase() for k, v in m.items()})
        out[name] = "".join(mm[c] for c in w)
    return out
for tol in (1e-12, 1e-13, 1e-14):
    sol = integrate(1.0, o["p1"], o["p2"], o["T"], rtol=tol, atol=tol)
    ts, cr = crossings(sol, o["T"])
    wa = free_reduce_cyclic(word_A(cr, ((1,2),1,1)))
    print(tol, len(cr), len(wa), "closure", np.linalg.norm(sol.sol(ts+o["T"])-sol.sol(ts)))
    for k, v in maps(wa).items():
        print("   ", k, canon(v) == canon(pub))
print("pub len", len(pub))
print("ours", wa)
print("pub ", pub)
