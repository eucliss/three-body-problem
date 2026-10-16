import re, sys, json
sys.path.insert(0, "/tmp/claude-0/-home-user-elementFi-takehome/abe89a27-8b28-5869-a9f7-5ac0ed0d524d/scratchpad/verify/topology")
from topo import canon
CAT = "/tmp/claude-0/-home-user-elementFi-takehome/abe89a27-8b28-5869-a9f7-5ac0ed0d524d/scratchpad/catalogs/liao/"
def swap(w, m):
    mm = dict(m); mm.update({k.upper(): v.swapcase() for k, v in m.items()})
    return "".join(mm[c] for c in w)
SYM = {"id": {'a': 'a', 'b': 'b'}, "a<->b": {'a': 'b', 'b': 'a'}, "a<->A": {'a': 'A', 'b': 'B'}, "a<->B": {'a': 'B', 'b': 'A'}}
entries = []
for fn in ["three-body-free-group-word.md", "three-body-unequal-mass-free-group-word.md", "free-fall-3b-free-group-word.md"]:
    txt = open(CAT + fn).read()
    for m in re.finditer(r'<th>(.*?)</th>\s*<td>([aAbB]+)</td>', txt, re.S):
        lab = re.sub(r'<[^>]+>', ' ', m.group(1)).replace('&nbsp;', ' ')
        entries.append((fn.split('-free')[0], " ".join(lab.split()), m.group(2)))
for fn in ["topological-sequences-of-3D-periodic-orbits.txt"]:
    pass
print("catalogue words parsed:", len(entries))
idx = {}
for e in entries:
    idx.setdefault(canon(e[2]), []).append(e)
C = json.load(open("/home/user/three-body-problem/research/results/claims.json"))
for c in C:
    hits = {}
    for nm, m in SYM.items():
        h = idx.get(canon(swap(c["word"], m)), [])
        if h:
            hits[nm] = [(x[0], x[1]) for x in h]
    print(c["id"], c["status"], "->", hits if hits else "no catalogue word equal (rot/inv, nor under a<->b, a<->A, a<->B)")
