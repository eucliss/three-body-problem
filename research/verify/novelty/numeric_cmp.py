"""Numeric comparison: for each claim, list catalogue entries at the same m3
with |T*-T*c|/T*c < 1% or max(|dv1|,|dv2|) < 1e-3 (also report nearest by v and T*).
Also recompute E and T* of each claim independently from (m3, v1, v2, T)."""
import sys, json, math
cat = json.load(open(sys.argv[1])); claims = json.load(open(sys.argv[2]))
words = {(x['file'].replace('free-group-word','X'), x['label']): x['word'] for x in cat if 'word' in x}
tab = []
for x in cat:
    if 'nums' in x and x['file'] == 'three-body-unequal-mass-pictures.md':
        v1, v2, T, Ts, Lf = x['nums'][:5]
        tab.append(dict(label=x['label'], m3=float(x['mass']), v1=v1, v2=v2, T=T, Ts=Ts, Lf=int(Lf),
                        word=words[(x['file'].replace('pictures','X'), x['label'])]))
eq = []
for x in cat:
    if 'nums' in x and x['file'] == 'three-body-pictures.md':
        v1, v2, T, Ts, Lf = x['nums'][:5]
        eq.append(dict(label=x['label'], m3=1.0, v1=v1, v2=v2, T=T, Ts=Ts, Lf=int(Lf),
                       word=words[(x['file'].replace('pictures','X'), x['label'])]))
def energy(m3, v1, v2):
    # r1=(-1,0), r2=(1,0), r3=(0,0); v1=v2=(v1,v2), v3=-2(v1,v2)/m3; G=1
    s = v1*v1 + v2*v2
    K = 0.5*(1*s + 1*s + m3*(4*s/(m3*m3)))
    U = -(1*1/2.0 + 1*m3/1.0 + 1*m3/1.0)
    return K + U
# sanity: catalogue T* reproduces from T and E
bad = 0
for e in tab + eq:
    E = energy(e['m3'], e['v1'], e['v2'])
    if E < 0:
        ts = e['T']*abs(E)**1.5
        if abs(ts - e['Ts']) > 2e-3*max(1, e['Ts']): bad += 1
    else: bad += 1
print('catalogue entries', len(tab), 'equal', len(eq), 'T* mismatches vs recomputed', bad)
for c in claims:
    E = energy(c['m3'], c['v1'], c['v2']); ts = c['T']*abs(E)**1.5
    print('\n%s m3=%g v=(%.6f,%.6f) T=%.4f T*=%.4f (recomputed E=%.10f vs %.10f, T*=%.6f)' %
          (c['id'], c['m3'], c['v1'], c['v2'], c['T'], c['Tstar'], E, c['E'], ts))
    same = [e for e in (eq if abs(c['m3']-1) < 1e-12 else tab) if abs(e['m3'] - c['m3']) < 1e-9]
    if not same:
        print('   no catalogue entries at this mass (masses in catalogue: %s)' % sorted({e['m3'] for e in tab}))
        continue
    hitsT = [e for e in same if abs(e['Ts'] - c['Tstar'])/e['Ts'] < 0.01]
    hitsV = [e for e in same if max(abs(e['v1']-c['v1']), abs(e['v2']-c['v2'])) < 1e-3]
    for e in hitsT: print('   T* within 1%%: %-22s v=(%.6f,%.6f) T=%.4f T*=%.3f Lf=%d word=%s' % (e['label'], e['v1'], e['v2'], e['T'], e['Ts'], e['Lf'], e['word']))
    for e in hitsV: print('   v within 1e-3: %-22s v=(%.6f,%.6f) T=%.4f T*=%.3f Lf=%d word=%s' % (e['label'], e['v1'], e['v2'], e['T'], e['Ts'], e['Lf'], e['word']))
    near = sorted(same, key=lambda e: math.hypot(e['v1']-c['v1'], e['v2']-c['v2']))[:3]
    for e in near: print('   nearest v: %-22s dv=%.2e v=(%.6f,%.6f) T=%.4f T*=%.3f Lf=%d word=%s' % (e['label'], math.hypot(e['v1']-c['v1'], e['v2']-c['v2']), e['v1'], e['v2'], e['T'], e['Ts'], e['Lf'], e['word'][:40]))
    # also: same T (period) within 1%? and same Lf?
    hitsP = [e for e in same if abs(e['T'] - c['T'])/e['T'] < 0.01]
    print('   entries with T within 1%%: %s' % [(e['label'], e['Lf'], round(e['Ts'],3)) for e in hitsP])
    print('   entries with Lf == %d: %s' % (c['Lf'], [(e['label'], round(e['Ts'],3), e['word']) for e in same if e['Lf'] == c['Lf']][:12]))
