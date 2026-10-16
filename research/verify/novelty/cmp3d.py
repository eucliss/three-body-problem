import re, sys, json, math
rows=[]
for fn in sys.argv[1:3]:
    for line in open(fn):
        m=re.match(r'\s*O_\{(\d+)\}\(([\d.]+)\)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S)', line)
        if m:
            rows.append(dict(n=int(m.group(1)), m3=float(m.group(2)), z0=float(m.group(3)), vx=float(m.group(4)), vy=float(m.group(5)), vz=float(m.group(6)), T=float(m.group(7)), st=m.group(8), f=fn.split('/')[-1][:20]))
claims=json.load(open(sys.argv[3]))
near_planar=[r for r in rows if abs(r['z0'])<0.02 and abs(r['vz'])<0.02]
print('3D rows', len(rows), 'near-planar (|z0|,|vz|<0.02):', len(near_planar))
for r in sorted(near_planar, key=lambda r:(r['m3'],r['n']))[:40]: print('  ', r)
for c in claims:
    same=[r for r in rows if abs(r['m3']-c['m3'])<1e-9]
    if not same:
        print(c['id'], 'no 3D entries at m3', c['m3']); continue
    best=sorted(same, key=lambda r: math.hypot(r['vx']-c['v1'], r['vy']-c['v2']))[:3]
    print(c['id'], 'm3', c['m3'], 'n3D', len(same), 'nearest (vx,vy):', [(r['n'], round(math.hypot(r['vx']-c['v1'], r['vy']-c['v2']),4), round(r['z0'],4), round(r['vz'],5), round(r['T'],3), round(r['T']/c['T'],3)) for r in best])
