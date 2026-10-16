"""For pairs of words, list exactly which (reverse?, letter map) transformations relate them up to rotation."""
INV={'a':'A','A':'a','b':'B','B':'b'}
def maps():
    out=[]
    for sab in (0,1):
        for ca in (0,1):
            for cb in (0,1):
                m={}
                for x in 'aAbB':
                    y=x
                    if ca and y.lower()=='a': y=INV[y]
                    if cb and y.lower()=='b': y=INV[y]
                    if sab: y={'a':'b','A':'B','b':'a','B':'A'}[y]
                    m[x]=y
                name=('swap_ab+' if sab else '')+('case_a+' if ca else '')+('case_b' if cb else '')
                out.append((name or 'id', m))
    return out
def rots(w): return {w[i:]+w[:i] for i in range(len(w))}
def relate(w1,w2):
    res=[]
    for rev in (0,1):
        u=w1[::-1] if rev else w1
        for name,m in maps():
            v=''.join(m[c] for c in u)
            if v in rots(w2): res.append(('rev+' if rev else '')+name)
    return res
P=lambda n:'Ba'*n+'bA'*n
Q=lambda n:'BA'*n+'ba'*n
pairs=[('L5 claim (k=5)','ABaBaBaBaBaBabAbAbAbAbAb','I.A.2(4)','BABABAbabababababaBABABA'),
       ('L1 claim','ABaBabAb','I.A.2(0.5)','BaBabAbA'),
       ('L2 claim','ABaBaBabAbAb','I.A.3(0.5)','BaBaBabAbAbA'),
       ('L2 claim','ABaBaBabAbAb','I.A.2(2)/I.A.1(4)','BABabababABA'),
       ('L1 claim','ABaBabAb','I.A.2 eq / I.A.1(2)','BAbabaBA'),
       ('eq I.A.3 (P2)','BaBabAbA','eq I.A.2 (Q2)','BAbabaBA'),
       ('P1','BabA','Q1','BAba'),
       ('X3 claim','aBaBaBabAbAbAbAB','L3 claim','ABaBaBaBabAbAbAb'),
       ('L3@0.1 claim','aBaBaBabAbAbAbAB','P4','BaBaBaBabAbAbAbA'),
       ('L4@0.1 claim','aBaBaBaBabAbAbAbAbAB','P5','BaBaBaBaBabAbAbAbAbA'),
       ('SD butterfly I','ababABAB','Q2 (LL I.A.2)','BAbabaBA')]
for a,wa,b,wb in pairs:
    print(f'{a:18s} vs {b:22s}:', relate(wa,wb) or 'NOT related')
for n in range(1,13):
    print('P%d vs Q%d:'%(n,n), relate(P(n),Q(n)))
