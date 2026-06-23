"""
ITERATION 4 (convergence test). Iter 3 showed: for a 74yo, the elaborate gut/novel
stack adds ~1 month of dialysis-free survival because COMPETING RISK dominates.
So the real question: does reducing the COMPETING HAZARD (AKI prevention, CV/overall
health, fall prevention, sick-day rules) move the needle more than toxin drugs?
This decides where the 1-year effort should actually go.
"""
import numpy as np
rng = np.random.default_rng(41)
N = 40000

CAT = {
 "base_medical":    dict(slope=(2.0,0.25), adh=0.10, feas=0.97),
 "gut_stack":       dict(slope=(1.0,0.5),  adh=0.20, feas=0.60),   # lumped creative gut stack
}
def eff_slope(name):
    s=CAT[name]; g=max(rng.normal(s["slope"][0], s["slope"][0]*s["slope"][1]),0)
    g*=np.clip(rng.normal(1-s["adh"],0.08),0.3,1.0)
    return g
def regimen_slope(names):
    gains=[]
    for n in names:
        if rng.random()>CAT[n]["feas"]: continue
        gains.append(eff_slope(n))
    gains.sort(reverse=True)
    return -4.0 + sum(g*(0.8**i) for i,g in enumerate(gains))

# competing hazard, scaled by a prevention factor (1.0=usual care, <1 = good AKI/CV prevention)
def simulate(names, hazard_scale=1.0, horizon_yr=25):
    reach=0; dfree=[]
    for _ in range(N):
        eGFR=25.0; sl=min(regimen_slope(names),-0.2); t=0.0; outcome=None
        for yr in range(horizon_yr):
            h=(0.09+0.012*yr)*hazard_scale
            if rng.random()<h: t=yr+rng.random(); outcome="comp"; break
            eGFR+=sl; t=yr+1
            if eGFR<=8: outcome="dial"; break
        if outcome=="dial": reach+=1
        dfree.append(t*12)
    return reach/N, np.median(dfree), np.percentile(dfree,[5,95])

scenarios = {
 "Base only (usual care)":                 (["base_medical"], 1.0),
 "Base + gut stack (toxin route)":         (["base_medical","gut_stack"], 1.0),
 "Base + AKI/CV prevention (-30% hazard)": (["base_medical"], 0.70),
 "Base + gut + prevention (both)":         (["base_medical","gut_stack"], 0.70),
 "Base + AGGRESSIVE prevention (-45%)":    (["base_medical"], 0.55),
}
print("="*90)
print("ITER 4   P(reach dialysis) | dialysis-free survival median [90% CI] | gain vs usual care")
print("="*90)
ref=None
for name,(sel,hs) in scenarios.items():
    pr,df,ci=simulate(sel,hs)
    if ref is None: ref=df
    print(f"{name:<42} P(dial)={pr*100:>3.0f}%  dfree={df:>4.0f} mo [{ci[0]:.0f}-{ci[1]:.0f}]  (+{df-ref:>3.0f} mo)")

print("\nRead: compare the +months from the 'toxin route' vs the 'prevention route'.")
