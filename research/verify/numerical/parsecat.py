"""Parse Li&Liao picture tables into JSON. usage: python3 -I parsecat.py <md> <out.json>"""
import re, sys, json
s = open(sys.argv[1], encoding='utf-8', errors='replace').read()
rows = re.split(r'<tr>', s)
out = []
for r in rows:
    cells = re.findall(r'<th[^>]*>(.*?)</th>', r, flags=re.S)
    if len(cells) < 6:
        continue
    name = re.sub(r'<[^>]+>', '', cells[0]).replace('&nbsp;', ' ').strip()
    try:
        v1, v2, T, Ts = (float(cells[i].strip()) for i in range(1, 5))
        Lf = int(cells[5].strip())
    except Exception:
        continue
    mm = re.search(r'\(([\d.]+)\)', name)
    out.append({'name': ' '.join(name.split()), 'm3': float(mm.group(1)) if mm else 1.0,
                'v1': v1, 'v2': v2, 'T': T, 'Tstar': Ts, 'Lf': Lf})
json.dump(out, open(sys.argv[2], 'w'), indent=0)
print(len(out), 'rows; masses:', sorted(set(o['m3'] for o in out))[:40])
