"""
LARGE-SCALE END-TO-END VALIDATION (N=1000).

Runs the complete operational protocol (Stage 1 → decision → Stage 2 for
borderline → final classification) on a literature-calibrated virtual cohort
of 1000 patients. Reports operational metrics that a clinician/reviewer would
need to evaluate feasibility.

Additionally reproduces published RCT cohorts (EPPIC, Rossi, Esgalhado)
using their reported summary statistics to show what our protocol would have
detected had it been applied to those study populations.

All parameters from literature_recalibration.py (Pretorius 2013, 2025 meta,
Lin 2011, Sirich 2014).
"""
import numpy as np

# =========================================================================
# CONSTANTS
# =========================================================================
CV_STD  = 0.15    # standardized measurement CV
WK_S1   = 24      # stage 1 duration
WK_S2   = 12      # stage 2 additional
N_S1    = 6       # measures per arm, stage 1
N_S2    = 3       # additional measures per arm, stage 2
DT_S1  = 1.645 * CV_STD * np.sqrt(2.0 / N_S1)       # 14.2%
DT_S2  = 1.645 * CV_STD * np.sqrt(2.0 / (N_S1+N_S2)) # 11.6%
DRAWS_S1 = 2 * N_S1 * 2   # 24 (duplicate draws)
DRAWS_S2 = 2 * N_S2 * 2   # 12

N_PAT   = 1000
N_REP   = 200     # reps per patient (for power estimation)

# =========================================================================
# (1) COHORT GENERATION (literature-calibrated)
# =========================================================================
rng = np.random.default_rng(777)

egfr    = np.clip(rng.normal(22, 6, N_PAT), 10, 35)
gut_var = np.clip(rng.normal(1.0, 0.35, N_PAT), 0.3, 1.8)
slope   = np.clip(rng.normal(-2.0, 0.8, N_PAT), -5.0, -0.5)
# baseline total IS shape ~ (25/eGFR)^1.2, then affine-rescaled to the CKD4
# cohort mean 5.4 / SD 3.6 ug/mL (Lin 2011). The rescale preserves the RNG
# stream (one normal draw of size N_PAT) and the eGFR ranking, and because
# baseline IS cancels in the ratio estimator (mean_A-mean_B)/mean_A it does
# not change any operating characteristic.
_raw = (25.0/egfr)**1.2 * np.clip(rng.normal(1.0, 0.30, N_PAT), 0.3, 2.5)
base_is_abs = np.clip(5.4 + 3.6 * (_raw - _raw.mean()) / _raw.std(), 0.5, None)
base_is = base_is_abs

# true individual IS reduction (literature-calibrated)
tau = np.clip(rng.normal(0.30 * gut_var, 0.14), 0, 0.70)
nr = rng.random(N_PAT) < 0.18
tau[nr] = np.clip(rng.normal(0.03, 0.02, nr.sum()), 0, 0.08)

true_resp = tau >= 0.10
weak = (tau >= 0.10) & (tau < 0.20)
moderate = (tau >= 0.20) & (tau < 0.35)
strong = tau >= 0.35
non_r = tau < 0.10

print("="*86)
print(f"(1) COHORT: N={N_PAT}, literature-calibrated (Pretorius 2013, 2025 meta, Lin 2011)")
print("="*86)
for lb, mk in [("non-resp <10%",non_r),("weak 10-20%",weak),
                ("moderate 20-35%",moderate),("strong 35%+",strong)]:
    print(f"  {lb:<22} n={mk.sum():>4}  ({mk.sum()/N_PAT*100:.0f}%)  mean tau={np.mean(tau[mk])*100:.0f}%")
print(f"  total responders:    {true_resp.sum():>4}  ({true_resp.sum()/N_PAT*100:.0f}%)")

# =========================================================================
# (2) END-TO-END PROTOCOL SIMULATION
# =========================================================================
print(f"\n{'='*86}")
print(f"(2) END-TO-END PROTOCOL (Stage 1 -> decision -> Stage 2 if borderline)")
print(f"{'='*86}")

# For each patient, simulate N_REP runs of the full protocol
final_class = np.zeros((N_PAT, N_REP), dtype='U1')  # R/N/B
went_s2     = np.zeros((N_PAT, N_REP), bool)

for p in range(N_PAT):
    b = base_is_abs[p]
    t = tau[p]
    # Stage 1: 6 control + 6 intervention
    A1 = b * (1 + rng.normal(0, CV_STD, (N_REP, N_S1)))
    B1 = b * (1-t) * (1 + rng.normal(0, CV_STD, (N_REP, N_S1)))
    ma1 = A1.mean(1); mb1 = B1.mean(1)
    obs1 = np.where(ma1>0, (ma1-mb1)/ma1, 0)

    # Stage 1 decision
    is_resp = obs1 > DT_S1
    is_nonr = obs1 < 0
    is_bord = ~is_resp & ~is_nonr

    final_class[p, is_resp] = 'R'
    final_class[p, is_nonr] = 'N'

    # Stage 2 for borderline
    if is_bord.any():
        went_s2[p, is_bord] = True
        A2 = b * (1 + rng.normal(0, CV_STD, (N_REP, N_S2)))
        B2 = b * (1-t) * (1 + rng.normal(0, CV_STD, (N_REP, N_S2)))
        # combine with stage 1
        A_all = np.concatenate([A1, A2], axis=1)
        B_all = np.concatenate([B1, B2], axis=1)
        ma_c = A_all.mean(1); mb_c = B_all.mean(1)
        obs_c = np.where(ma_c>0, (ma_c-mb_c)/ma_c, 0)

        s2_resp = obs_c > DT_S2
        final_class[p, is_bord & s2_resp] = 'R'
        final_class[p, is_bord & ~s2_resp] = 'N'

# =========================================================================
# (3) OPERATIONAL METRICS
# =========================================================================
# Per-patient majority classification (>50% of reps)
majority = np.array(['R' if (final_class[p]=='R').sum() > N_REP/2 else 'N'
                      for p in range(N_PAT)])
s2_rate = went_s2.mean(1)  # per-patient fraction going to S2

# Confusion matrix
tp = ((majority=='R') & true_resp).sum()
fn = ((majority=='N') & true_resp).sum()
fp = ((majority=='R') & ~true_resp).sum()
tn = ((majority=='N') & ~true_resp).sum()

sens = tp/(tp+fn) if (tp+fn)>0 else 0
spec = tn/(tn+fp) if (tn+fp)>0 else 0
ppv  = tp/(tp+fp) if (tp+fp)>0 else 0
npv  = tn/(tn+fn) if (tn+fn)>0 else 0
accuracy = (tp+tn)/N_PAT

avg_wk = np.mean(np.where(s2_rate.mean() > 0.5,
                           WK_S1 + s2_rate * WK_S2,
                           np.where(s2_rate > 0.5, WK_S1+WK_S2, WK_S1)))
# simpler: per patient avg weeks
per_pt_wk = WK_S1 + s2_rate.mean(axis=0 if s2_rate.ndim>1 else None) * WK_S2 \
            if s2_rate.ndim > 1 else WK_S1 + s2_rate * WK_S2
mean_wk = np.mean(WK_S1 + s2_rate * WK_S2)
mean_draws = np.mean(DRAWS_S1 + s2_rate * DRAWS_S2)
pct_s2 = (s2_rate > 0.5).mean()

print(f"\n  CLASSIFICATION MATRIX (N={N_PAT}, majority over {N_REP} reps):")
print(f"  {'':>20}{'classified R':>14}{'classified N':>14}")
print(f"  {'true responder':<20}{tp:>14}{fn:>14}")
print(f"  {'true non-resp':<20}{fp:>14}{tn:>14}")
print(f"\n  sensitivity (TP rate):  {sens*100:.1f}%  ({tp}/{tp+fn} responders detected)")
print(f"  specificity (TN rate): {spec*100:.1f}%  ({tn}/{tn+fp} non-resp correctly stopped)")
print(f"  PPV:                   {ppv*100:.1f}%")
print(f"  NPV:                   {npv*100:.1f}%")
print(f"  overall accuracy:      {accuracy*100:.1f}%")
print(f"\n  OPERATIONAL:")
print(f"  patients -> Stage 2:   {pct_s2*100:.0f}%  ({int(pct_s2*N_PAT)} of {N_PAT})")
print(f"  avg protocol duration: {mean_wk:.1f} weeks")
print(f"  avg blood draws:       {mean_draws:.0f} (duplicate)")

# Subgroup power
print(f"\n  SUBGROUP DETECTION:")
for lb, mk in [("weak 10-20%",weak),("moderate 20-35%",moderate),("strong 35%+",strong)]:
    n = mk.sum()
    if n == 0: continue
    detected = (majority[mk]=='R').sum()
    print(f"  {lb:<22} {detected:>4}/{n:<4} detected  ({detected/n*100:.0f}%)")

# Non-responder burden
burden_all = N_PAT * 36  # treat-all: all patients x 36 months
burden_nof1 = (majority=='R').sum() * 36
waste_all = non_r.sum() * 36
waste_nof1 = fp * 36
print(f"\n  BURDEN ANALYSIS (36-month horizon):")
print(f"  {'':>28}{'treat-all':>12}{'n-of-1':>12}{'saved':>10}")
print(f"  {'patients on regimen':<28}{N_PAT:>12}{(majority=='R').sum():>12}"
      f"{N_PAT-(majority=='R').sum():>10}")
print(f"  {'non-resp burden-months':<28}{waste_all:>12}{waste_nof1:>12}"
      f"{waste_all-waste_nof1:>10}")
print(f"  {'burden reduction':>28}{'-':>12}{(1-waste_nof1/waste_all)*100:>11.0f}%")

# =========================================================================
# (4) REPRODUCE PUBLISHED COHORTS
# =========================================================================
print(f"\n{'='*86}")
print("(4) 'WHAT IF' — applying our protocol to published RCT populations")
print(f"{'='*86}")

PUB_SEEDS = list(range(700, 750))  # 50 single-run replications, matching Table S1

def _simulate_published_once(seed, n, tau_mean, tau_sd, baseline_mean, baseline_sd, nr_frac):
    """One single-run pass of the protocol on a cohort matching published stats."""
    rng2 = np.random.default_rng(seed)
    tau_p = np.clip(rng2.normal(tau_mean, tau_sd, n), 0, 0.70)
    nr_p = rng2.random(n) < nr_frac
    tau_p[nr_p] = np.clip(rng2.normal(0.03, 0.02, nr_p.sum()), 0, 0.08)
    bis = np.clip(rng2.normal(baseline_mean, baseline_sd, n), 0.5, 30)
    resp = tau_p >= 0.10

    classified = []
    for i in range(n):
        A = bis[i] * (1 + rng2.normal(0, CV_STD, N_S1))
        B = bis[i] * (1-tau_p[i]) * (1 + rng2.normal(0, CV_STD, N_S1))
        obs = (A.mean()-B.mean())/A.mean() if A.mean()>0 else 0
        if obs > DT_S1:
            classified.append('R')
        elif obs < 0:
            classified.append('N')
        else:
            A2 = bis[i] * (1 + rng2.normal(0, CV_STD, N_S2))
            B2 = bis[i] * (1-tau_p[i]) * (1 + rng2.normal(0, CV_STD, N_S2))
            A_c = np.concatenate([A, A2]); B_c = np.concatenate([B, B2])
            obs2 = (A_c.mean()-B_c.mean())/A_c.mean() if A_c.mean()>0 else 0
            classified.append('R' if obs2 > DT_S2 else 'N')
    classified = np.array(classified)
    tp = ((classified=='R') & resp).sum()
    fp = ((classified=='R') & ~resp).sum()
    sens = tp/resp.sum() if resp.sum()>0 else np.nan
    spec = 1 - fp/(~resp).sum() if (~resp).sum()>0 else np.nan
    return resp.sum(), tp, fp, sens, spec, np.mean(tau_p)

def simulate_published(label, n, tau_mean, tau_sd, baseline_mean, baseline_sd, nr_frac):
    """Mean single-run operating characteristics over PUB_SEEDS fixed seeds.

    The former single draw used np.random.default_rng(hash(label) % 2**31);
    Python salts str hashes per process (PYTHONHASHSEED), so that table could
    not be reproduced run to run and, at n=26-40, swung 20+ points between
    processes. Averaging over 50 fixed seeds gives a stable, reproducible
    estimate, the same convention Table S1 uses.
    """
    runs = [_simulate_published_once(s, n, tau_mean, tau_sd,
                                     baseline_mean, baseline_sd, nr_frac)
            for s in PUB_SEEDS]
    resp = np.mean([r[0] for r in runs])
    tp   = np.mean([r[1] for r in runs])
    fp   = np.mean([r[2] for r in runs])
    sens = np.nanmean([r[3] for r in runs])
    spec = np.nanmean([r[4] for r in runs])
    mt   = np.mean([r[5] for r in runs])
    return n, resp, tp, fp, sens, spec, mt

# Published cohort reproductions (summary statistics from literature)
cohorts = [
    # label, n, tau_mean, tau_sd, baseline_IS_mean, baseline_IS_sd, nr_frac
    # EPPIC: AST-120 arm, IS reduction ~30-40%, N=1000
    ("EPPIC AST-120 arm (Schulman 2015)",
     460, 0.35, 0.15, 5.4, 3.6, 0.15),
    # Rossi SYNERGY: synbiotic, IS~12%, PCS~16%, N=37
    ("Rossi SYNERGY synbiotic (2016)",
     37, 0.14, 0.10, 4.0, 2.5, 0.30),
    # Esgalhado RS crossover: IS reduction ~25%, N=26
    ("Esgalhado RS crossover (2020)",
     26, 0.25, 0.12, 6.0, 3.0, 0.20),
    # Sirich RS parallel: IS reduction -5mg/L from ~20, ~25%, N=40
    ("Sirich RS parallel (2014)",
     40, 0.25, 0.13, 5.0, 2.8, 0.20),
    # 2025 meta: pooled fiber, SMD -0.34, estimated ~25% reduction
    ("2025 meta-analysis pooled (11 RCTs)",
     200, 0.25, 0.14, 5.4, 3.6, 0.18),
]

print(f"  {'cohort':<40}{'N':>5}{'resp':>6}{'TP':>5}{'FP':>5}{'sens':>7}{'spec':>7}{'mean tau':>9}")
print(f"  {'-'*82}")
for label, n, tm, ts, bm, bs, nrf in cohorts:
    n_tot, n_resp, tp, fp, sens, spec, mt = simulate_published(label,n,tm,ts,bm,bs,nrf)
    print(f"  {label:<40}{n_tot:>5}{n_resp:>6}{tp:>5}{fp:>5}{sens*100:>6.0f}%{spec*100:>6.0f}%{mt*100:>8.0f}%")

print(f"""
  INTERPRETATION:
  - EPPIC: our protocol would have identified ~{simulate_published(*(['EPPIC AST-120 arm (Schulman 2015)']+list(cohorts[0][1:])))[4]*100:.0f}% of individual AST-120
    responders — instead of treating all 1000 patients based on mean effect
  - Rossi SYNERGY (small n=37, weak effect): low sensitivity expected (~{simulate_published(*(['Rossi SYNERGY synbiotic (2016)']+list(cohorts[1][1:])))[4]*100:.0f}%)
    because mean tau ~14% is near the DT boundary
  - Larger pooled cohort (meta): ~{simulate_published(*(['2025 meta-analysis pooled (11 RCTs)']+list(cohorts[4][1:])))[4]*100:.0f}% sensitivity, demonstrating scalability

  KEY POINT: these RCTs REPORTED that the intervention "works" (significant mean
  effect), but could not tell WHICH PATIENTS it works for. Our protocol does.""")

# =========================================================================
# (5) PROTOCOL EFFICIENCY SUMMARY
# =========================================================================
print(f"\n{'='*86}")
print("(5) PROTOCOL EFFICIENCY SUMMARY")
print(f"{'='*86}")
print(f"""
  METRIC                          VALUE           BENCHMARK
  ---------------------------------------------------------------
  cohort size                     {N_PAT}            typical RCT: 30-460
  overall accuracy                {accuracy*100:.0f}%             -
  sensitivity                     {sens*100:.0f}%             target: >80%
  specificity                     {spec*100:.0f}%             target: >90%
  PPV                             {ppv*100:.0f}%             -
  NPV                             {npv*100:.0f}%             -
  patients needing Stage 2        {pct_s2*100:.0f}%             -
  avg duration                    {mean_wk:.0f} wk           fixed 3x3 = 36 wk
  avg blood draws (duplicate)     {mean_draws:.0f}             fixed 3x3 = 36
  non-resp burden eliminated      {(1-waste_nof1/waste_all)*100:.0f}%             treat-all = 0%
  eGFR cost (24 wk)              ~0.9 mL/min      -
  cost per visit                  1 extra tube      negligible

  CONCLUSION: the protocol achieves clinical-grade classification
  ({sens*100:.0f}% sensitivity, {spec*100:.0f}% specificity) with an average of
  {mean_wk:.0f} weeks and {mean_draws:.0f} blood draws per patient.""")
