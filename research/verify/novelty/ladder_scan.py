"""Scan catalogue for words of the ladder shapes X^n Y^n with X,Y two-letter syllables,
strictly up to rotation only (no relabelling), and report which (X,Y,n) occur at which mass."""
import sys, json, re
cat = json.load(open(sys.argv[1]))
nums = {(x['file'].replace('pictures','W'), x['label']): x['nums'] for x in cat if 'nums' in x}
def rots(w): return {w[i:]+w[:i] for i in range(len(w))}
syl = [x+y for x in 'aAbB' for y in 'aAbB' if x.lower()!=y.lower()]
out = []
for r in cat:
    if 'word' not in r: continue
    w = r['word']; L = len(w)
    if L % 4: continue
    n = L//4
    R = rots(w)
    for X in syl:
        for Y in syl:
            if X==Y: continue
            if X*n+Y*n in R:
                key=(r['file'].replace('free-group-word','W'), r['label'])
                out.append((X,Y,n,r['label'],r['file'][:22],w, nums.get(key)))
for o in sorted(set((a,b,c,d,e,f,tuple(g) if g else None) for a,b,c,d,e,f,g in out)):
    print(o)
