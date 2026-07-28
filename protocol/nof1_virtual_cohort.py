"""
n-of-1 VIRTUAL COHORT PROTOCOL SIMULATION (N=100).

Tests whether a repeated-crossover n-of-1 protocol can DETECT individual
treatment effects of the converged gut-clearance regimen, given realistic:
  - within-person IS measurement noise (CV=0.22)
  - heterogeneous true effects (15% non-responder sub-population)
  - eGFR drift during the protocol (crossover-controlled)

Patient: CKD stage 3b-4, sampled from plausible distributions.
Design: AB crossover (A=base medical only, B=base + gut regimen).
  Period 4wk on-treatment, washout 2wk, K measurements in final week.
Primary: serum indoxyl sulfate (IS) % reduction.
Decision: observed IS drop > DT threshold (one-sided 95%).
"""

import numpy as np
rng = np.random.default_rng(101)

N_PAT  = 100
CV     = 0.22               # within-person IS measurement CV (assay + biology)
WK_T   = 4                  # treatment period (weeks)
WK_W   = 2                  # washout between periods (weeks)
WK_CYC = 2*(WK_T + WK_W)   # 12 wk per AB cycle
N_REP  = 500                # Monte Carlo reps per patient per design

# =========================================================================
# (1) VIRTUAL COHORT: 100 patients from stage 3b-4 CKD
# =========================================================================
age     = np.clip(rng.normal(74, 5, N_PAT), 60, 85).astype(int)
eGFR    = np.clip(rng.normal(22, 6, N_PAT), 10, 35)
# baseline IS: rises as eGFR falls; inter-patient variability in gut production
base_IS = (25.0/eGFR)**1.2 * np.clip(rng.normal(1.0, 0.20, N_PAT), 0.5, 1.8) * 100
# gut microbiome response variability (patient-level random effect on treatment)
gut_var = np.clip(rng.normal(1.0, 0.30, N_PAT), 0.3, 1.8)
# eGFR slope after base medical treatment (mL/min/yr, negative)
slope   = np.clip(rng.normal(-2.0, 0.8, N_PAT), -5.0, -0.5)
# annual competing hazard (death/CV/AKI)
hazard  = 0.05 + 0.002*(age - 60) + 0.003*np.clip(35 - eGFR, 0, 25)
# f_tox: toxin-mortality coupling (uncertain; prior skewed near 0)
f_tox   = np.clip(rng.beta(1.5, 6.0, N_PAT), 0, 0.5)

print("="*82)
print("(1) VIRTUAL COHORT (N=100, CKD stage 3b-4)")
print("="*82)
for v, nm in [(age,"age"), (eGFR,"eGFR (mL/min)"), (base_IS,"baseline IS (rel%)"),
              (slope,"eGFR slope/yr"), (hazard,"annual hazard"), (f_tox,"f_tox")]:
    print(f"  {nm:<22} mean={np.mean(v):>6.1f}  sd={np.std(v):>5.1f}"
          f"  [{np.min(v):>5.1f} - {np.max(v):>5.1f}]")

# =========================================================================
# (2) TRUE INDIVIDUAL TREATMENT EFFECTS (hidden from protocol)
# =========================================================================
tau = np.clip(rng.normal(0.35 * gut_var, 0.12), 0, 0.70)
non_resp = rng.random(N_PAT) < 0.15
n_nr = non_resp.sum()
tau[non_resp] = np.clip(rng.normal(0.03, 0.02, n_nr), 0, 0.08)
true_resp = tau >= 0.10

print(f"\n{'='*82}")
print(f"(2) TRUE IS REDUCTION (hidden): mean={np.mean(tau)*100:.0f}%  "
      f"sd={np.std(tau)*100:.0f}%  non-resp={n_nr}/{N_PAT}")
print(f"{'='*82}")
for lo, hi, lb in [(0,.10,"<10% (non-resp)"), (.10,.25,"10-25% (weak)"),
                    (.25,.45,"25-45% (moderate)"), (.45,1.0,"45%+ (strong)")]:
    m = (tau >= lo) & (tau < hi); n = m.sum()
    if n: print(f"  {lb:<24} n={n:>3}  mean true τ={np.mean(tau[m])*100:.0f}%")

# =========================================================================
# (3) SIMULATION ENGINE
# =========================================================================
def run(nc, km, drift=True):
    """Simulate n-of-1 crossover.  nc=AB cycles, km=measures per period.
    Returns (obs[N_PAT,N_REP], total_weeks, dt_threshold)."""
    n_arm = nc * km
    obs = np.zeros((N_PAT, N_REP))
    for p in range(N_PAT):
        a_drift, b_drift = [], []
        wk = 0
        for _ in range(nc):
            wk_a = wk + WK_T
            if drift:
                eg = max(eGFR[p] + slope[p]*wk_a/52, 5)
                da = (25./eg)**1.2 / (25./eGFR[p])**1.2
            else:
                da = 1.0
            a_drift.extend([da]*km)
            wk = wk_a + WK_W

            wk_b = wk + WK_T
            if drift:
                eg = max(eGFR[p] + slope[p]*wk_b/52, 5)
                db = (25./eg)**1.2 / (25./eGFR[p])**1.2
            else:
                db = 1.0
            b_drift.extend([db]*km)
            wk = wk_b + WK_W

        ad = np.array(a_drift)[None, :]
        bd = np.array(b_drift)[None, :]
        A = base_IS[p] * ad * (1 + rng.normal(0, CV, (N_REP, n_arm)))
        B = base_IS[p] * bd * (1 - tau[p]) * (1 + rng.normal(0, CV, (N_REP, n_arm)))
        ma, mb = A.mean(axis=1), B.mean(axis=1)
        obs[p] = np.where(ma > 0, (ma - mb) / ma, 0)

    dt = 1.645 * CV * np.sqrt(2.0 / n_arm)
    return obs, nc * WK_CYC, dt

# =========================================================================
# (4) PROTOCOL COMPARISON: power × FP × duration
# =========================================================================
designs = [(1,2), (1,3), (2,2), (2,3), (3,2), (3,3)]
RES = {}
for nc, km in designs:
    o, wk, dt = run(nc, km)
    det = (o > dt).mean(axis=1)
    RES[(nc,km)] = dict(obs=o, wk=wk, dt=dt, det=det,
                        pw=det[true_resp].mean(), fp=det[~true_resp].mean())

print(f"\n{'='*82}")
print("(4) PROTOCOL COMPARISON  (power=TP among true responders, FP among non-resp)")
print(f"{'='*82}")
print(f"  {'design':<15}{'power':>8}{'FP':>8}{'DT':>8}{'weeks':>8}{'~months':>9}")
print(f"  {'-'*56}")
for nc, km in designs:
    r = RES[(nc,km)]
    tag = " <--" if (nc,km) == (2,3) else ""
    print(f"  {nc}cyc x {km}meas{r['pw']*100:>7.0f}%{r['fp']*100:>7.0f}%"
          f"{r['dt']*100:>7.0f}%{r['wk']:>7}{r['wk']/4.3:>8.0f}{tag}")

D = RES[(2,2)]   # default for downstream analysis

# =========================================================================
# (5) POWER BY TRUE EFFECT SIZE (2cyc x 2meas)
# =========================================================================
print(f"\n{'='*82}")
print(f"(5) DETECTION POWER BY TRUE EFFECT SIZE (2x2, DT={D['dt']*100:.0f}%)")
print(f"{'='*82}")
for lo, hi in [(0,.10), (.10,.20), (.20,.30), (.30,.40), (.40,.70)]:
    m = (tau >= lo) & (tau < hi); n = m.sum()
    if n == 0: continue
    pw = D['det'][m].mean()
    print(f"  tau {lo*100:>2.0f}-{hi*100:<3.0f}%  n={n:>3}  power={pw*100:>5.0f}%  {'#'*int(pw*30)}")

# =========================================================================
# (6) eGFR DRIFT IMPACT (analytical bias + power comparison)
# =========================================================================
# Analytical: B periods measured later -> eGFR lower -> IS higher -> effect underestimated
drift_bias = np.zeros(N_PAT)
for p in range(N_PAT):
    a_d, b_d = [], []
    wk = 0
    for _ in range(2):
        wk_a = wk + WK_T
        eg = max(eGFR[p] + slope[p]*wk_a/52, 5)
        a_d.append((25./eg)**1.2 / (25./eGFR[p])**1.2)
        wk = wk_a + WK_W
        wk_b = wk + WK_T
        eg = max(eGFR[p] + slope[p]*wk_b/52, 5)
        b_d.append((25./eg)**1.2 / (25./eGFR[p])**1.2)
        wk = wk_b + WK_W
    drift_bias[p] = (1 - tau[p]) * (1 - np.mean(b_d)/np.mean(a_d))

o_nd, _, dt_nd = run(2, 2, drift=False)
det_nd = (o_nd > dt_nd).mean(axis=1)
pw_drift = D['det'][true_resp].mean()
pw_nodr  = det_nd[true_resp].mean()

print(f"\n{'='*82}")
print(f"(6) eGFR DRIFT IMPACT (2x2 design, analytical)")
print(f"{'='*82}")
print(f"  drift bias (underestimates effect):")
print(f"    median = {np.median(drift_bias)*100:.1f} pp"
      f"   90% range: [{np.percentile(drift_bias,5)*100:.1f} to"
      f" {np.percentile(drift_bias,95)*100:.1f}] pp")
print(f"  power: with drift = {pw_drift*100:.0f}%  |  without drift = {pw_nodr*100:.0f}%"
      f"  (delta = {abs(pw_drift-pw_nodr)*100:.0f} pp)")
print(f"  -> A-B-A-B crossover keeps drift bias within ~1 pp for typical patients")
print(f"     improvement: A-B-B-A ordering cancels linear drift completely")

# =========================================================================
# (7) BAYESIAN SHRINKAGE: empirical Bayes, averaged over 100 protocol runs
# =========================================================================
se_i = CV * np.sqrt(2.0 / 4)    # 4 measures per arm in 2x2
N_EB = 100
eb_rmse_raw = np.zeros(N_EB)
eb_rmse_shr = np.zeros(N_EB)
eb_rmse_raw_grp = {g: np.zeros(N_EB) for g in ["nr","weak","mod+"]}
eb_rmse_shr_grp = {g: np.zeros(N_EB) for g in ["nr","weak","mod+"]}
grp_masks = {"nr": ~true_resp, "weak": (tau>=.10)&(tau<.25), "mod+": tau>=.25}

for rep in range(N_EB):
    obs_r = D['obs'][:, rep]
    mu_r = np.mean(obs_r)
    vp_r = max(np.var(obs_r) - se_i**2, 0.001)
    w_r = vp_r / (vp_r + se_i**2)
    sh_r = w_r * obs_r + (1 - w_r) * mu_r
    eb_rmse_raw[rep] = np.sqrt(np.mean((obs_r - tau)**2))
    eb_rmse_shr[rep] = np.sqrt(np.mean((sh_r - tau)**2))
    for g, mk in grp_masks.items():
        if mk.sum() > 0:
            eb_rmse_raw_grp[g][rep] = np.sqrt(np.mean((obs_r[mk] - tau[mk])**2))
            eb_rmse_shr_grp[g][rep] = np.sqrt(np.mean((sh_r[mk] - tau[mk])**2))

# keep first run for illustrative patients section
obs1 = D['obs'][:, 0]
mu_hat = np.mean(obs1)
var_pop = max(np.var(obs1) - se_i**2, 0.001)
w = var_pop / (var_pop + se_i**2)
shrunk = w * obs1 + (1 - w) * mu_hat

rr_avg = np.mean(eb_rmse_raw)
rs_avg = np.mean(eb_rmse_shr)

print(f"\n{'='*82}")
print(f"(7) EMPIRICAL BAYES SHRINKAGE (averaged over {N_EB} single-run executions)")
print(f"{'='*82}")
print(f"  shrinkage weight w ~ {w:.2f}  (0=trust pop mean only, 1=trust data only)")
print(f"  RMSE raw (avg):      {rr_avg*100:.1f}%")
print(f"  RMSE shrunk (avg):   {rs_avg*100:.1f}%  ({(1-rs_avg/rr_avg)*100:+.0f}%)")
grp_labels = {"nr": "non-resp", "weak": "weak 10-25%", "mod+": "moderate+"}
for g in ["nr","weak","mod+"]:
    mk = grp_masks[g]; n = mk.sum()
    if n == 0: continue
    rr = np.mean(eb_rmse_raw_grp[g])
    rs = np.mean(eb_rmse_shr_grp[g])
    chg = (1-rs/rr)*100 if rr > 0 else 0
    print(f"    {grp_labels[g]:<14} n={n:>3}  raw={rr*100:>5.1f}%  shrunk={rs*100:>5.1f}%  ({chg:+.0f}%)")

print(f"\n  interpretation: simple EB helps WEAK responders (near pop mean) but HURTS")
print(f"  non-responders (pulled up) and strong responders (pulled down).")
print(f"  -> for bimodal populations: use mixture-model EB, not simple EB.")
print(f"  -> for clinical decisions: threshold-based detection (section 4-5) is more robust.")

# =========================================================================
# (8) ILLUSTRATIVE PATIENTS (3 archetypes)
# =========================================================================
print(f"\n{'='*82}")
print(f"(8) THREE ILLUSTRATIVE PATIENTS (2x2 design)")
print(f"{'='*82}")
idxs = [("STRONG",   int(np.argmax(tau * true_resp.astype(float)))),
        ("MODERATE",  int(np.argmin(np.abs(tau - 0.28)))),
        ("NON-RESP",  int(np.argmin(np.abs(tau - 0.03))))]
for lb, i in idxs:
    med = np.median(D['obs'][i])
    q5, q95 = np.percentile(D['obs'][i], [5, 95])
    print(f"\n  [{lb}] pt#{i}: age={age[i]}, eGFR={eGFR[i]:.0f}, "
          f"baseline IS={base_IS[i]:.0f}, gut_var={gut_var[i]:.2f}")
    print(f"    true tau       = {tau[i]*100:.0f}%")
    print(f"    observed (500x): median={med*100:.0f}%  90%CI=[{q5*100:.0f}-{q95*100:.0f}%]")
    print(f"    single run:      raw={obs1[i]*100:.0f}%  shrunk={shrunk[i]*100:.0f}%")
    print(f"    detection power  = {D['det'][i]*100:.0f}%  (DT={D['dt']*100:.0f}%)")

# =========================================================================
# (9) COHORT-LEVEL CLASSIFICATION MATRIX (2x2 and 2x3 designs)
# =========================================================================
print(f"\n{'='*82}")
print(f"(9) CLASSIFICATION MATRIX: protocol correctly sorts responders vs non-resp")
print(f"{'='*82}")
for tag, key in [("2x2 (24 wk)", (2,2)), ("2x3 (24 wk)", (2,3))]:
    r = RES[key]
    n_tr = true_resp.sum(); n_nr = (~true_resp).sum()
    tp = (r['det'][true_resp]  > 0.50).sum()   # patient detected in >50% of reps -> "detected"
    fn = n_tr - tp
    fp = (r['det'][~true_resp] > 0.50).sum()
    tn = n_nr - fp
    ppv = tp/(tp+fp) if (tp+fp) > 0 else 0
    npv = tn/(tn+fn) if (tn+fn) > 0 else 0
    print(f"\n  {tag}  (measures/arm={key[0]*key[1]}, DT={r['dt']*100:.0f}%)")
    print(f"    true responders: {n_tr}   non-responders: {n_nr}")
    print(f"    TP={tp:>3}  FN={fn:>3}  |  FP={fp:>3}  TN={tn:>3}")
    print(f"    sensitivity={tp/n_tr*100:.0f}%  specificity={tn/n_nr*100:.0f}%"
          f"  PPV={ppv*100:.0f}%  NPV={npv*100:.0f}%")

# =========================================================================
# (10) PROTOCOL RECOMMENDATION
# =========================================================================
r23 = RES[(2,3)]
egfr_drop = abs(np.mean(slope)) * 24 / 52

print(f"\n{'='*82}")
print(f"(10) PROTOCOL RECOMMENDATION")
print(f"{'='*82}")
print(f"""
  RECOMMENDED DESIGN: 2 cycles x 3 measures/period (24 weeks, ~6 months)
    power = {r23['pw']*100:.0f}%  |  FP = {r23['fp']*100:.0f}%  |  DT = {r23['dt']*100:.0f}%

  PROTOCOL TIMELINE (A-B-A-B crossover):
    wk  0- 4 : A1 (control)        serum IS x3 at wk 3-4
    wk  4- 6 : washout
    wk  6-10 : B1 (gut regimen)    serum IS x3 at wk 9-10
    wk 10-12 : washout
    wk 12-16 : A2 (control)        serum IS x3 at wk 15-16
    wk 16-18 : washout
    wk 18-22 : B2 (gut regimen)    serum IS x3 at wk 21-22
    wk 22-24 : final washout + analysis

  DECISION RULE:
    IS drop > {r23['dt']*100:.0f}% -> RESPONDER  (continue regimen)
    IS drop < {r23['dt']*100:.0f}% -> NON-RESP   (deprescribe, eliminate burden)

  STEPWISE ADDITION (optional, extends to ~12 months):
    first  2 cycles: fiber sachet alone
    next   2 cycles: add probiotic (Lacto, bedtime dosing)
    -> each component independently evaluated via its own DT test

  eGFR COST: ~{egfr_drop:.1f} mL/min lost during 24-wk protocol (crossover-controlled)
  BLOOD DRAWS: 12 total (3 per period x 4 periods)

  VALUE:
    correct responder ID for {N_PAT} patients
    -> {true_resp.sum()} responders kept on effective regimen (power {r23['pw']*100:.0f}%)
    -> {(~true_resp).sum()} non-resp freed from useless burden (spec {(1-r23['fp'])*100:.0f}%)

  THE HINGE (from why_experts_failed_revalidation.py):
  this protocol is the SINGLE differentiator from the failed surrogate story.
  without n-of-1 measurement -> we are just another unvalidated guess.""")
