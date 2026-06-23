"""
RED-TEAM REFINEMENT ENGINE. Repeats the critique->fix->rescore loop greedily:
each round an 'expert' attacks the WEAKEST component of the current design, the
best available fix is applied, robustness R is recomputed. Runs until convergence.

Design quality components (0..1): benefit (risk-adjusted realized), safety,
adherence, honesty (evidence integrity), and burden (lower better).
Robustness R rewards the WEAKEST link (multiplicative) -> forces breadth, not depth.
Each fix improves its target but may cost nominal potency or another component:
the honest engineering tradeoff.
"""
import numpy as np

state = dict(benefit=0.55, safety=0.45, adherence=0.50, honesty=0.40, burden=0.55, potency=0.62)

def R(s):
    return 100*(s["benefit"]**0.30 * s["safety"]**0.25 * s["adherence"]**0.20 *
                s["honesty"]**0.15 * (1-s["burden"])**0.10)

# fix pool: each targets a component (+improve toward 1), with optional side costs.
# (critique, fix, target, improve, side_costs[list of (comp, delta)])
FIXES = [
 ("Surrogate fallacy: AST-120 trial neg on hard endpoints","Reframe primary outcome to symptom/QoL+dialysis-prob; f_tox=0 base case","honesty",0.55,[]),
 ("AST-120 adsorbs her antihypertensives (#1 lever)","Drop AST-120 from base regimen","safety",0.45,[("potency",-0.12),("benefit",-0.02)]),
 ("Spores sacrifice documented urea function","Keep functional Lacto + bedtime dosing","benefit",0.20,[]),
 ("eGFR slope is non-linear; dialysis-date extrapolation unreliable","Report ranges only, no dialysis dates; not a clinical predictor","honesty",0.35,[]),
 ("Ketoacid Ca-load -> vascular calcification","Cap ketoacid dose; monitor serum Ca; diet-first","safety",0.30,[("potency",-0.03)]),
 ("Protein 0.69 g/kg risks frailty in 74yo","Relax to 0.75 g/kg; ketoacid optional","safety",0.25,[("potency",-0.04)]),
 ("Hyperkalemia: fiber+RAS+diet underestimated","K monitoring 2-4wk; cap fiber; low-K blend; K-binder if needed","safety",0.30,[]),
 ("Fiber bloating -> adherence collapse","Slow 3-wk titration; acacia-weighted; split AM/PM","adherence",0.35,[]),
 ("Parameters fabricated; no external validation","Tie every claim to n-of-1 measurement; Bayesian shrink to null; wide CIs","honesty",0.45,[("benefit",-0.03)]),
 ("Competing risk dominates at 74; low ceiling","Deprescribe low-yield arms; prioritize QoL + AKI prevention","burden",-0.30,[("benefit",0.04)]),
 ("Probiotic CFU/viability/storage variability","Specify validated product, refrigerated, CFU-verified","adherence",0.20,[]),
 ("Broad drug-supplement timing interactions","Full medication timing audit; 2h spacing rule","safety",0.20,[]),
 ("n-of-1 regression-to-mean / placebo","Washout + replicate measures + randomized intro order","honesty",0.25,[]),
 ("Cumulative pill/sachet burden unsustainable","Minimum-effective-set; quarterly deprescribe review","burden",-0.30,[("potency",-0.03)]),
 ("Self-experiment iatrogenic risk","Safety stop-rules + supervision triggers (K, Ca, weight loss)","safety",0.20,[]),
 ("AKI single-event can erase years of slope gain","Sick-day rules + nephrotoxin/contrast/dehydration guard","benefit",0.18,[]),
]

def apply(s, fix):
    s=dict(s); _,_,tgt,imp,costs=fix
    if tgt=="burden":
        s["burden"]=max(0.05, s["burden"]+imp)      # imp negative lowers burden
    else:
        s[tgt]=min(0.98, s[tgt]+imp*(1-s[tgt]))      # diminishing returns toward 1
    for comp,d in costs:
        s[comp]=float(np.clip(s[comp]+d,0.05,0.98))
    return s

print("="*92)
print("RED-TEAM REFINEMENT LOOP (greedy: attack weakest link each round)")
print("="*92)
print(f"{'rd':>3}{'R':>7}  {'component attacked':<20}{'critique -> fix'}")
print(f"{'--':>3}{R(state):>7.1f}  {'(start)':<20}")

used=[]; trajectory=[R(state)]
for rd in range(1,17):
    # weakest component (the expert's target)
    comps={k:state[k] for k in ["benefit","safety","adherence","honesty"]}
    comps["burden"]=1-state["burden"]
    weakest=min(comps,key=comps.get)
    # choose the available fix giving the biggest R gain
    best=None
    for i,fx in enumerate(FIXES):
        if i in used: continue
        cand=apply(state,fx); gain=R(cand)-R(state)
        if best is None or gain>best[0]: best=(gain,i,fx,cand)
    if best is None or best[0]<0.15:
        print(f"\nConverged at round {rd-1}: remaining fixes add <0.15 R.")
        break
    gain,i,fx,cand=best; used.append(i); state=cand; trajectory.append(R(state))
    print(f"{rd:>3}{R(state):>7.1f}  {fx[2]:<20}{fx[0][:36]} -> {fx[1][:34]}")

print("\n" + "="*92)
print(f"FINAL  R = {R(state):.1f}  (from {trajectory[0]:.1f})   rounds applied = {len(used)}")
print("Final components:", {k:round(state[k],2) for k in state})
print(f"\nNote the honest tradeoff: nominal potency {0.62:.2f} -> {state['potency']:.2f} (IS reduction)")
print("but robustness (safety/adherence/honesty/burden) rose sharply -> higher RISK-ADJUSTED value.")
print("Trajectory R:", [round(x,1) for x in trajectory])
