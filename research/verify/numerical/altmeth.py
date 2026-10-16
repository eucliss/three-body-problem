"""Cross-check closure with a different integrator family (implicit Radau IIA / LSODA)."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
from nb import ic, make_rhs
claims = json.load(open(sys.argv[1])); meth = sys.argv[2]; tol = float(sys.argv[3]); ids = sys.argv[4:]
for c in claims:
    if ids and c['id'] not in ids: continue
    m = (1.0, 1.0, c['m3']); Y0 = ic(c['m3'], c['v1'], c['v2'])
    t0 = time.time()
    sol = solve_ivp(make_rhs(m), (0, c['T']), Y0, method=meth, rtol=tol, atol=tol)
    print(f"{c['id']:<28} {meth} tol={tol:g} res={np.max(np.abs(sol.y[:,-1]-Y0)):.2e} nsteps={len(sol.t)} t={time.time()-t0:.0f}s", flush=True)
