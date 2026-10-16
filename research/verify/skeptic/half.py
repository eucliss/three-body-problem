import json, sys, numpy as np
from scipy.integrate import solve_ivp
C=json.load(open(sys.argv[1]))
def rhs(t,y,m):
    r=y[:6].reshape(3,2); v=y[6:]; a=np.zeros((3,2))
    for i in range(3):
        for j in range(3):
            if i!=j:
                d=r[j]-r[i]; a[i]+=m[j]*d/np.linalg.norm(d)**3
    return np.concatenate([v,a.ravel()])
for e in C:
    if e['id'] not in ('X3@0.25','L3@0.25','L3@0.5','L5@0.5','L4@0.1'): continue
    m3=e['m3']; m=[1,1,m3]; v=np.array([e['v1'],e['v2']]); T=e['T']
    y0=np.concatenate([[-1,0,1,0,0,0],v,v,-2*v/m3])
    s=solve_ivp(rhs,(0,T),y0,args=(m,),method='DOP853',rtol=1e-13,atol=1e-13,dense_output=True)
    yT=s.y[:,-1]; yh=s.sol(T/2)
    r=yh[:6].reshape(3,2); vv=yh[6:].reshape(3,2)
    # Fix(R): r3=0, r1=-r2, v1=v2
    fixres=max(np.abs(r[2]).max(),np.abs(r[0]+r[1]).max(),np.abs(vv[0]-vv[1]).max())
    scale=np.linalg.norm(r[1]-r[0])/2
    print(e['id'],'closure %.1e'%np.abs(yT-y0).max(),'T/2 in Fix(R) residual %.1e'%fixres,'half-size %.4f'%scale)
