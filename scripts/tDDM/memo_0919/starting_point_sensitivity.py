# Sensitivity of the probability of acceptance to the starting point, per trial,
# from each participant's best-run 4j fit (common random numbers, w = 0.5 +/- 0.05).
# Converts the group difference in acceptance rate (and its 95% CI) into a
# starting-point difference. Run starting_point_export.R first. 2026-09-27.
import numpy as np, pandas as pd, sys, time
S='results/tDDM/memo_analyses_0919/starting_point_sensitivity_20260927'; NP=int(sys.argv[1]) if len(sys.argv)>1 else 400; DT=float(sys.argv[2]) if len(sys.argv)>2 else 0.002
f=pd.read_csv('results/tDDM/canonical_4j_accurate_bestof10/fits_4j_accurate_best10.csv')
t=pd.read_csv(S+'/trials.csv').merge(f[['subject','d_mean','d_contrast','thres','nDT','tau']],on='subject')
SG=1.4; MX=10.0
v1=(t.d_mean*t.mean_sig).values; v2=(t.d_mean*t.mean_sig+t.d_contrast*t.contrast_sig).values
a=(2*t.thres).values; tau=t.tau.values; n=len(t)
def pup(w,seed):
    rng=np.random.default_rng(seed); res=np.zeros(n)
    for s in range(0,n,500):
        idx=np.arange(s,min(n,s+500)); m=len(idx)
        A=np.repeat(a[idx],NP); V1=np.repeat(v1[idx],NP); V2=np.repeat(v2[idx],NP); TAU=np.repeat(tau[idx],NP)
        x=A*w; up=np.zeros(m*NP); done=np.zeros(m*NP,bool); tt=0.0
        while tt<MX and not done.all():
            v=np.where((TAU>0)&(tt<TAU),V1,V2)
            x=x+np.where(done,0,v*DT+SG*np.sqrt(DT)*rng.standard_normal(m*NP)); tt+=DT
            hitu=(~done)&(x>=A); hitl=(~done)&(x<=0)
            up[hitu]=1; done|=hitu|hitl
        up[~done]=(x[~done]>A[~done]/2)
        res[idx]=up.reshape(m,NP).mean(1)
    return res
t0=time.time(); h=0.05
pp=pup(0.5+h,1); pm=pup(0.5-h,1); p0=pup(0.5,1)
t['sens']=(pp-pm)/(2*h); t['p0']=p0
ps=t.groupby(['subject','group']).agg(sens=('sens','mean'),p0=('p0','mean')).reset_index()
fbar=ps.sens.mean()
print(f'NP={NP} dt={DT} time={time.time()-t0:.0f}s')
print('mean predicted P(accept) at w=.5:',round(ps.p0.mean(),3),' observed:',round((t.choice_ddm==1).mean(),3))
print('sensitivity dP/dw: mean',round(fbar,4),' HC',round(ps[ps.group=="HC"].sens.mean(),4),' DEP',round(ps[ps.group=="DEP"].sens.mean(),4), ' sd across subj',round(ps.sens.std(),3))
for D in (0.0159,0.0663):
    print(f'acceptance diff {D:.4f} -> starting-point diff {D/fbar:.4f}')
