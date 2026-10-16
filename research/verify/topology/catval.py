import re, sys, json, time
import numpy as np
sys.path.insert(0, "/tmp/claude-0/-home-user-elementFi-takehome/abe89a27-8b28-5869-a9f7-5ac0ed0d524d/scratchpad/verify/topology")
from topo import *
CAT = "/tmp/claude-0/-home-user-elementFi-takehome/abe89a27-8b28-5869-a9f7-5ac0ed0d524d/scratchpad/catalogs/liao/"
key_re = re.compile(r'<th>([IV]+\.[A-D])<sup>.*?<sub[^>]*>(\d+)</sub>\s*&nbsp;\(([\d.]+)\)</th>')
def parse_words(fn):
    txt = open(CAT + fn).read()
    out = {}
    for m in re.finditer(r'<tr>\s*' + key_re.pattern + r'\s*<td>([aAbB]+)</td>', txt):
        out[(m.group(1), int(m.group(2)), m.group(3))] = m.group(4)
    return out
def parse_pics(fn):
    txt = open(CAT + fn).read()
    out = {}
    pat = r'<tr>\s*' + key_re.pattern + r'\s*<th>\s*([-\d.eE]+)</th>\s*<th>\s*([-\d.eE]+)</th>\s*<th>\s*([-\d.eE]+)</th>\s*<th>\s*([-\d.eE]+)</th>\s*<th>(\d+)</th>'
    for m in re.finditer(pat, txt):
        out[(m.group(1), int(m.group(2)), m.group(3))] = dict(v1=float(m.group(4)), v2=float(m.group(5)), T=float(m.group(6)), Ts=float(m.group(7)), Lf=int(m.group(8)))
    return out
W = parse_words("three-body-unequal-mass-free-group-word.md")
P = parse_pics("three-body-unequal-mass-pictures.md")
print("parsed words", len(W), "pics", len(P), "joined", len(set(W) & set(P)))
if __name__ == "__main__":
    nper = int(sys.argv[1]); masses = sys.argv[2].split(",")
    WIN = [((1, 2), 1, +1), ((1, 2), 2, -1)]
    tot = good = 0
    bad = []
    for ms in masses:
        keys = sorted([k for k in set(W) & set(P) if k[2] == ms], key=lambda k: P[k]["T"])
        # spread over classes: take the nper shortest per class
        chosen = []
        for cl in sorted(set(k[0] for k in keys)):
            chosen += [k for k in keys if k[0] == cl][:nper]
        for k in chosen:
            p = P[k]; m3 = float(ms)
            try:
                sol = integrate(m3, p["v1"], p["v2"], p["T"], rtol=1e-13, atol=1e-13)
            except Exception as e:
                try:
                    sol = integrate(m3, p["v1"], p["v2"], p["T"], rtol=1e-13, atol=1e-13, tmax=1.02 * p["T"])
                except Exception as e:
                    print(k, "INTEGRATION FAILED", e); continue
            ts, cr = crossings(sol, p["T"])
            clos = np.linalg.norm(sol.sol(ts + p["T"]) - sol.sol(ts))
            tgt = canon(W[k])
            ws = [free_reduce_cyclic(word_A(cr, c)) for c in WIN]
            ok = [canon(w) == tgt for w in ws]
            tot += 1; good += all(ok)
            if not all(ok):
                bad.append((k, W[k], ws, ok, clos))
            print(k, "T=%.3f" % p["T"], "Lf", p["Lf"], len(W[k]), "match", ok, "closure %.1e" % clos, flush=True)
    print("TOTAL", tot, "all-conv match", good)
    for b in bad:
        print("MISMATCH", b)
