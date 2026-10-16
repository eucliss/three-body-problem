import re, sys, json
D='/tmp/claude-0/-home-user-elementFi-takehome/abe89a27-8b28-5869-a9f7-5ac0ed0d524d/scratchpad/catalogs/liao/'
files=['three-body-free-group-word.md','three-body-unequal-mass-free-group-word.md','free-fall-3b-free-group-word.md']
def red(w):
    s=[]
    for c in w:
        if s and s[-1]!=c and s[-1].lower()==c.lower(): s.pop()
        else: s.append(c)
    while len(s)>1 and s[0]!=s[-1] and s[0].lower()==s[-1].lower(): s=s[1:-1]
    return ''.join(s)
def mp(w,m): return ''.join(m[c] for c in w)
ID={c:c for c in 'aAbB'}
CASE={'a':'A','A':'a','b':'B','B':'b'}
SWAP={'a':'b','b':'a','A':'B','B':'A'}
SWC={'a':'B','B':'a','b':'A','A':'b'}
PHYS=[ID,CASE,SWAP,SWC]
GEN=PHYS+[{'a':'A','A':'a','b':'b','B':'B'},{'a':'a','A':'A','b':'B','B':'b'},
          {'a':'b','b':'A','A':'B','B':'a'},{'a':'B','B':'A','A':'b','b':'a'}]
def key(w,maps):
    w=red(w); best=None
    for m in maps:
        for v in (mp(w,m), mp(w,m)[::-1]):
            for i in range(len(v)):
                r=v[i:]+v[:i]
                if best is None or r<best: best=r
    return best
cat=[]
for f in files:
    t=open(D+f).read()
    for lab,w in re.findall(r'<th>(.*?)</th>\s*<td>([aAbB]+)</td>',t,re.S):
        lab=re.sub(r'<[^>]+>','',lab).replace('&nbsp;',' ')
        cat.append((f.split('-')[0]+':'+' '.join(lab.split()),w))
print('catalogue words',len(cat))
P=lambda n:'Ba'*n+'bA'*n
Q=lambda n:'BA'*n+'ba'*n
for name,maps in (('PHYS',PHYS),('GEN8',GEN)):
    kc={}
    for lab,w in cat: kc.setdefault(key(w,maps),[]).append(lab)
    for n in range(1,14):
        print(name,'P%d'%n, kc.get(key(P(n),maps),[])[:6])
claims=json.load(open(sys.argv[1]))
kc={}
for lab,w in cat: kc.setdefault(key(w,PHYS),[]).append(lab)
for e in claims:
    w=e['word']; n=len(red(w))//4
    print(e['id'], 'Lf',len(red(w)), 'isP%d'%n, key(w,PHYS)==key(P(n),PHYS), 'catPHYS:',kc.get(key(w,PHYS),[])[:4])
