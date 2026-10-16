"""Independent parser for the Li & Liao catalogue markdown/HTML tables.

Usage: python3 -I parse_cat.py <catalog_dir> <out_json>
Produces a JSON list of {file, label, cls, num, mass, word, (v1, v2, T, Tstar, Lf)}.
"""
import sys, re, json, os

catdir, out = sys.argv[1], sys.argv[2]

TAG = re.compile(r'<[^>]+>')

def clean(s):
    s = TAG.sub(' ', s)
    s = s.replace('&nbsp;', ' ')
    return ' '.join(s.split())

def rows(path):
    txt = open(path, encoding='utf-8').read()
    # split on <tr>
    parts = re.split(r'<tr>', txt)
    for p in parts[1:]:
        cells = re.findall(r'<(th|td)[^>]*>(.*?)</\1>', p, flags=re.S)
        yield [c[1] for c in cells]

LABEL = re.compile(r'^\s*(.*?)\s*\(([\d.,\s]+)\)\s*$')

def parse_label(raw):
    # raw html of first cell, e.g. 'I.A<sup><i><i>i.c.</i></i></sup><sub ...>1</sub> &nbsp;(0.5)'
    c = clean(raw)
    m = LABEL.match(c)
    mass = None
    head = c
    if m:
        head, mass = m.group(1), m.group(2).replace(' ', '')
    return c, head, mass

records = []
for fn, kind in [('three-body-free-group-word.md', 'word'),
                 ('three-body-unequal-mass-free-group-word.md', 'word'),
                 ('free-fall-3b-free-group-word.md', 'word'),
                 ('three-body-pictures.md', 'num'),
                 ('three-body-unequal-mass-pictures.md', 'num'),
                 ('free-fall-3b-pictures.md', 'num')]:
    path = os.path.join(catdir, fn)
    n = 0
    for cells in rows(path):
        if not cells:
            continue
        label, head, mass = parse_label(cells[0])
        if 'Class' in label or label.startswith('F i (m') or 'Free group' in label:
            continue
        rec = {'file': fn, 'label': label, 'head': head, 'mass': mass}
        if kind == 'word':
            if len(cells) < 2:
                continue
            w = clean(cells[1]).replace(' ', '')
            if not w:
                continue
            rec['word'] = w
        else:
            vals = [clean(c) for c in cells[1:]]
            nums = []
            for v in vals:
                try:
                    nums.append(float(v))
                except ValueError:
                    pass
            rec['nums'] = nums
        records.append(rec)
        n += 1
    print(fn, n, file=sys.stderr)

json.dump(records, open(out, 'w'), indent=0)
