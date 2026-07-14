"""
RESPONDER-THRESHOLD SENSITIVITY.
The 10% true-responder definition is a modeling choice. Here we hold the protocol
(2x3 adaptive, CV=0.15, canonical seed 777 -- identical to Table 2) fixed and only
change the threshold theta that defines a "true responder" (tau >= theta), then
recompute operating characteristics. This shows how much the reported sens/spec
depend on the (arbitrary) responder definition.
"""
import numpy as np

N_PAT, CV_STD, N_S1, N_S2 = 1000, 0.15, 6, 3

# canonical cohort (seed 777, identical to reproduce_manuscript_numbers.py)
rng = np.random.default_rng(777)
egfr = np.clip(rng.normal(22, 6, N_PAT), 10, 35)
gut = np.clip(rng.normal(1.0, 0.35, N_PAT), 0.3, 1.8)
slope = np.clip(rng.normal(-2.0, 0.8, N_PAT), -5.0, -0.5)
bis = (25.0/egfr)**1.2 * np.clip(rng.normal(1.0, 0.30, N_PAT), 0.3, 2.5) * 5.4
tau = np.clip(rng.normal(0.30 * gut, 0.14), 0, 0.70)
nr = rng.random(N_PAT) < 0.18
tau[nr] = np.clip(rng.normal(0.03, 0.02, nr.sum()), 0, 0.08)

MDE_S1 = 1.645 * CV_STD * np.sqrt(2.0/N_S1)
MDE_S2 = 1.645 * CV_STD * np.sqrt(2.0/(N_S1+N_S2))

# run the adaptive protocol ONCE (identical to Table 2)
rng_run = np.random.default_rng(777)
cls = np.full(N_PAT, 'N', dtype='U1')
for p in range(N_PAT):
    b, t = bis[p], tau[p]
    A1 = b*(1+rng_run.normal(0,CV_STD,N_S1//2)); B1 = b*(1-t)*(1+rng_run.normal(0,CV_STD,N_S1//2))
    A2 = b*(1+rng_run.normal(0,CV_STD,N_S1//2)); B2 = b*(1-t)*(1+rng_run.normal(0,CV_STD,N_S1//2))
    As, Bs = np.concatenate([A1,A2]), np.concatenate([B1,B2])
    obs1 = (As.mean()-Bs.mean())/As.mean() if As.mean()>0 else 0
    if obs1 > MDE_S1: cls[p]='R'
    elif obs1 < 0:    cls[p]='N'
    else:
        A3 = b*(1+rng_run.normal(0,CV_STD,N_S2)); B3 = b*(1-t)*(1+rng_run.normal(0,CV_STD,N_S2))
        Ac, Bc = np.concatenate([As,A3]), np.concatenate([Bs,B3])
        obs2 = (Ac.mean()-Bc.mean())/Ac.mean() if Ac.mean()>0 else 0
        cls[p] = 'R' if obs2 > MDE_S2 else 'N'

print("="*72)
print("RESPONDER-THRESHOLD SENSITIVITY (protocol fixed; truth relabeled)")
print("="*72)
print(f"{'theta':>7}{'true resp':>11}{'true NR':>9}{'sens %':>9}{'spec %':>9}")
for theta in [0.05, 0.10, 0.15, 0.20]:
    tr = tau >= theta
    tp = ((cls=='R') & tr).sum(); fn = ((cls=='N') & tr).sum()
    fp = ((cls=='R') & ~tr).sum(); tn = ((cls=='N') & ~tr).sum()
    sens = tp/(tp+fn)*100 if (tp+fn)>0 else 0
    spec = tn/(tn+fp)*100 if (tn+fp)>0 else 0
    mark = "  <- reference" if theta==0.10 else ""
    print(f"{theta*100:>6.0f}%{tr.sum():>11}{(~tr).sum():>9}{sens:>9.1f}{spec:>9.1f}{mark}")

print("\nRead: sensitivity rises and specificity falls as the responder bar is raised,")
print("because a higher theta moves near-threshold patients from the responder to the")
print("non-responder pool. The protocol's detection threshold (MDE=14%) is unchanged.")
