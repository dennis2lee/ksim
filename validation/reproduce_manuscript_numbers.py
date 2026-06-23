"""
SINGLE REPRODUCIBILITY SCRIPT: regenerates every number cited in manuscript.md.

Run this once to verify all manuscript tables and inline statistics.
Any mismatch between this output and the manuscript = a bug.
"""
import numpy as np

N_PAT = 1000
CV_STD = 0.15
CV_NAT = 0.25
N_S1 = 6
N_S2 = 3

# =========================================================================
# COHORT (same as robustness_experiments.py baseline, seed=777)
# =========================================================================
rng = np.random.default_rng(777)
egfr = np.clip(rng.normal(22, 6, N_PAT), 10, 35)
gut = np.clip(rng.normal(1.0, 0.35, N_PAT), 0.3, 1.8)
slope = np.clip(rng.normal(-2.0, 0.8, N_PAT), -5.0, -0.5)
bis = (25.0/egfr)**1.2 * np.clip(rng.normal(1.0, 0.30, N_PAT), 0.3, 2.5) * 5.4
tau = np.clip(rng.normal(0.30 * gut, 0.14), 0, 0.70)
nr = rng.random(N_PAT) < 0.18
tau[nr] = np.clip(rng.normal(0.03, 0.02, nr.sum()), 0, 0.08)

true_resp = tau >= 0.10
weak = true_resp & (tau < 0.20)
moderate = true_resp & (tau >= 0.20) & (tau < 0.35)
strong = tau >= 0.35
non_r = ~true_resp

print("="*70)
print("COHORT CHARACTERISTICS (manuscript section: Cohort characteristics)")
print("="*70)
print(f"  non-resp (<10%):  {non_r.sum()}")
print(f"  weak (10-20%):    {weak.sum()}")
print(f"  moderate (20-35%): {moderate.sum()}")
print(f"  strong (>=35%):   {strong.sum()}")
print(f"  total:            {N_PAT}")
print(f"  mean tau (resp):  {tau[true_resp].mean()*100:.0f}%")

# =========================================================================
# TABLE 1: MDE values
# =========================================================================
print(f"\n{'='*70}")
print("TABLE 1: MDE by design and CV")
print(f"{'='*70}")
for design, cv, n, wk in [(("2x3",0.25,6,24)),("2x3",0.15,6,24),
                            ("3x3",0.25,9,36),("4x3",0.25,12,48),("3x3",0.15,9,36)]:
    mde = 1.645 * cv * np.sqrt(2.0/n)
    print(f"  {design} CV={cv}  n={n}  MDE={mde*100:.0f}%  wk={wk}")

# =========================================================================
# TABLE 2: Single-run confusion matrix (full adaptive protocol)
# =========================================================================
print(f"\n{'='*70}")
print("TABLE 2: Single-run classification (seed=777)")
print(f"{'='*70}")

MDE_S1 = 1.645 * CV_STD * np.sqrt(2.0/N_S1)
MDE_S2 = 1.645 * CV_STD * np.sqrt(2.0/(N_S1+N_S2))

rng2 = np.random.default_rng(777)
# consume same state as robustness baseline
rng2_cohort = np.random.default_rng(777)
_ = np.clip(rng2_cohort.normal(22, 6, N_PAT), 10, 35)
_ = np.clip(rng2_cohort.normal(1.0, 0.35, N_PAT), 0.3, 1.8)
_ = np.clip(rng2_cohort.normal(-2.0, 0.8, N_PAT), -5.0, -0.5)
_ = (25.0/_)**1.2 * np.clip(rng2_cohort.normal(1.0, 0.30, N_PAT), 0.3, 2.5) * 5.4
_ = np.clip(rng2_cohort.normal(0.30 * np.clip(rng2_cohort.normal(1.0, 0.35, N_PAT), 0.3, 1.8), 0.14), 0, 0.70)

# Run single protocol execution
rng_run = np.random.default_rng(777)
cls = np.full(N_PAT, 'N', dtype='U1')
went_s2 = np.zeros(N_PAT, bool)

for p in range(N_PAT):
    b = bis[p]; t = tau[p]
    A1 = b * (1 + rng_run.normal(0, CV_STD, N_S1//2))
    B1 = b * (1-t) * (1 + rng_run.normal(0, CV_STD, N_S1//2))
    A2 = b * (1 + rng_run.normal(0, CV_STD, N_S1//2))
    B2 = b * (1-t) * (1 + rng_run.normal(0, CV_STD, N_S1//2))
    A_s1 = np.concatenate([A1, A2])
    B_s1 = np.concatenate([B1, B2])
    obs1 = (A_s1.mean()-B_s1.mean())/A_s1.mean() if A_s1.mean()>0 else 0
    if obs1 > MDE_S1:
        cls[p] = 'R'
    elif obs1 < 0:
        cls[p] = 'N'
    else:
        went_s2[p] = True
        A3 = b * (1 + rng_run.normal(0, CV_STD, N_S2))
        B3 = b * (1-t) * (1 + rng_run.normal(0, CV_STD, N_S2))
        Ac = np.concatenate([A_s1, A3])
        Bc = np.concatenate([B_s1, B3])
        obs2 = (Ac.mean()-Bc.mean())/Ac.mean() if Ac.mean()>0 else 0
        cls[p] = 'R' if obs2 > MDE_S2 else 'N'

tp = ((cls=='R') & true_resp).sum()
fn = ((cls=='N') & true_resp).sum()
fp = ((cls=='R') & non_r).sum()
tn = ((cls=='N') & non_r).sum()
sens = tp/(tp+fn)
spec = tn/(tn+fp)

print(f"  True resp: {true_resp.sum()},  True NR: {non_r.sum()}")
print(f"  TP={tp}  FN={fn}  FP={fp}  TN={tn}")
print(f"  Sensitivity: {sens*100:.1f}%")
print(f"  Specificity: {spec*100:.1f}%")
print(f"  Stage 2 rate: {went_s2.mean()*100:.0f}%")

# TABLE 3: subgroup detection
print(f"\n  Subgroup detection (single-run):")
for lb, mk in [("weak", weak), ("moderate", moderate), ("strong", strong)]:
    n = mk.sum()
    det = ((cls=='R') & mk).sum()
    print(f"    {lb}: {det}/{n} ({det/n*100:.0f}%)")

# =========================================================================
# KEY INLINE NUMBERS
# =========================================================================
print(f"\n{'='*70}")
print("INLINE NUMBERS CHECK")
print(f"{'='*70}")
avg_wk = 24 + went_s2.mean() * 12
avg_draws = 24 + went_s2.mean() * 12  # duplicate draws
print(f"  Avg duration: {avg_wk:.1f} wk")
print(f"  MDE at CV=0.15, n=6: {MDE_S1*100:.1f}%")
print(f"  MDE at CV=0.15, n=9: {MDE_S2*100:.1f}%")
print(f"  MDE at CV=0.25, n=6: {1.645*0.25*np.sqrt(2/6)*100:.1f}%")

print(f"\n  FP rate (single-run): {fp/(fp+tn)*100:.1f}%")
print(f"  Approx NR fraction: {non_r.sum()/N_PAT*100:.1f}%")
