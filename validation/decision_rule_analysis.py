"""
Decision-rule analysis: aligning the responder definition with the statistical
null (review item R1) and controlling the two-stage overall type I error
(review item R2).

The original protocol labelled a patient a true responder at tau >= 10% but
tested H0: tau <= 0. Under that mismatch a patient whose true effect is, say,
6% is a "false positive" whenever the protocol correctly detects a real
non-zero reduction. Sharper measurement then makes the false-positive count go
up, which is what produced the apparent "optimal CV" near 0.12-0.15. That
optimum is an artifact of the mismatch, not a property of the measurement.

This script reports three rules side by side:

  zero-null      H0: tau <= 0     (the original manuscript rule)
  matched-null   H0: tau <= 0.10  (the null that matches the 10% definition)
  matched-null + alpha-spending, calibrated so the overall one-sided type I
                 error across both stages is 5%

Outputs feed manuscript Table 2, Table 4 and Tables S1, S8, S9.
"""

import numpy as np
from nof1_core import (make_cohort, run_protocol, critical_value, obf_alphas,
                       calibrate_boundaries, weighted_accuracy, baseline_moments)

N_PAT = 1000
COHORT_SEED = 777
MEAS_SEED = 777
THETA = 0.10
REP_SEEDS = range(700, 750)          # 50 replications, as elsewhere in the paper

cohort = make_cohort(N_PAT, COHORT_SEED, 0.18)


def replicated(cv, null_margin, alpha1=0.05, alpha2=0.05, seeds=REP_SEEDS,
               order=None, **kw):
    """Mean operating characteristics over replicate measurement draws.

    The cycle order is tied to the null: the zero-effect null is reported as it
    was originally specified, with the fixed control-first order, while the
    primary rule uses the within-patient randomization that the crossover-order
    analysis selects. Passing `order` overrides that pairing.
    """
    if order is None:
        order = 'fixed' if null_margin == 0.0 else 'randomized'
    runs = [run_protocol(cohort, cv=cv, measurement_seed=s, theta=THETA,
                         null_margin=null_margin, alpha1=alpha1, alpha2=alpha2,
                         order=order, **kw) for s in seeds]
    keys = ('sens', 'spec', 'fp_rate', 'stage2_rate', 'mean_weeks', 'ncc')
    out = {k: float(np.nanmean([r[k] for r in runs])) for k in keys}
    out['sens_sd'] = float(np.nanstd([r['sens'] for r in runs]))
    out['spec_sd'] = float(np.nanstd([r['spec'] for r in runs]))
    out['weak'] = float(np.nanmean([r['detect']['weak'] for r in runs]))
    out['tp'] = float(np.mean([r['tp'] for r in runs]))
    out['fp'] = float(np.mean([r['fp'] for r in runs]))
    out['fn'] = float(np.mean([r['fn'] for r in runs]))
    out['tn'] = float(np.mean([r['tn'] for r in runs]))
    return out


def banner(t):
    print('\n' + '=' * 78)
    print(t)
    print('=' * 78)


banner('COHORT (shared by every analysis below)')
m, s = baseline_moments(cohort)
print(f'  baseline IS realized mean {m:.3f} ug/mL, SD {s:.3f} ug/mL')
print('  ' + '  '.join(f'{k}={int(v.sum())}'
                       for k, v in cohort.strata(THETA).items()))

# ---------------------------------------------------------------------------
# R2: calibrate the two-stage boundaries under each null
# ---------------------------------------------------------------------------
banner('R2. ALPHA-SPENDING CALIBRATION (overall one-sided type I error = 5%)')
info = [6 / 9, 1.0]
print(f'  O\'Brien-Fleming nominal per-look alphas at information fractions '
      f'{info[0]:.3f}, {info[1]:.3f}: '
      + ', '.join(f'{a:.4f}' for a in obf_alphas(0.05, info)))

calib = {}
for margin, label in ((0.0, 'zero-null'), (THETA, 'matched-null')):
    cal = calibrate_boundaries(cohort, cv=0.15, null_margin=margin,
                               alpha_total=0.05, seeds=range(900, 920),
                               order='fixed' if margin == 0.0 else 'randomized')
    calib[margin] = cal
    print(f'\n  {label}: calibrated alpha1={cal["alpha1"]:.4f}, '
          f'alpha2={cal["alpha2"]:.4f} (scale {cal["scale"]:.3f})')
    print(f'    achieved overall type I error at tau = {margin:.2f}: '
          f'{cal["achieved"]*100:.1f}%')
    print(f'    Stage 1 critical value {critical_value(0.15, 6, cal["alpha1"], margin)*100:.1f}%, '
          f'Stage 2 {critical_value(0.15, 9, cal["alpha2"], margin)*100:.1f}%')

# uncorrected overall type I error, for the before/after comparison
for margin, label in ((0.0, 'zero-null'), (THETA, 'matched-null')):
    null_cohort = make_cohort(N_PAT, COHORT_SEED, 0.18)
    null_cohort.tau = np.full(N_PAT, margin)
    rates = []
    for sd in range(900, 950):
        r = run_protocol(null_cohort, cv=0.15, measurement_seed=sd,
                         theta=margin + 1e-9, null_margin=margin,
                         alpha1=0.05, alpha2=0.05,
                         order='fixed' if margin == 0.0 else 'randomized')
        rates.append((r['cls'] == 'R').mean())
    print(f'  {label}: UNCORRECTED overall type I error at tau = {margin:.2f}: '
          f'{np.mean(rates)*100:.1f}%')

# ---------------------------------------------------------------------------
# R1 + R2: the three rules compared at the reference CV
# ---------------------------------------------------------------------------
banner('R1/R2. THREE DECISION RULES AT CV = 0.15 '
       '(mean of 50 replications, seeds 700-749)')
rules = [
    ('zero-null, alpha 0.05 per look (original)', 0.0, 0.05, 0.05),
    ('zero-null, alpha-spending', 0.0,
     calib[0.0]['alpha1'], calib[0.0]['alpha2']),
    ('matched-null, alpha 0.05 per look', THETA, 0.05, 0.05),
    ('matched-null, alpha-spending (primary)', THETA,
     calib[THETA]['alpha1'], calib[THETA]['alpha2']),
]
print(f'  {"rule":<44}{"sens":>7}{"spec":>7}{"FP":>7}{"weak":>7}'
      f'{"S2":>6}{"wk":>7}')
print('  ' + '-' * 76)
rule_results = {}
for label, margin, a1, a2 in rules:
    r = replicated(0.15, margin, a1, a2)
    rule_results[label] = r
    print(f'  {label:<44}{r["sens"]*100:>6.1f}%{r["spec"]*100:>6.1f}%'
          f'{r["fp_rate"]*100:>6.1f}%{r["weak"]*100:>6.0f}%'
          f'{r["stage2_rate"]*100:>5.0f}%{r["mean_weeks"]:>7.1f}')

# ---------------------------------------------------------------------------
# R1: CV sweep under each null -- does the "optimal CV" survive?
# ---------------------------------------------------------------------------
banner('R1. CV SWEEP UNDER EACH NULL (mean of 50 replications)')
cvs = [0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.36]
print(f'  {"CV":>6} | {"zero-null":^30} | {"matched-null + spending":^30}')
print(f'  {"":>6} | {"sens":>7}{"spec":>7}{"NCC":>8}{"dNCC":>8} | '
      f'{"sens":>7}{"spec":>7}{"NCC":>8}{"dNCC":>8}')
print('  ' + '-' * 72)
sweep = {}
prev = {0.0: None, THETA: None}
for cv in cvs:
    row = {}
    for margin in (0.0, THETA):
        a1, a2 = ((0.05, 0.05) if margin == 0.0
                  else (calib[THETA]['alpha1'], calib[THETA]['alpha2']))
        row[margin] = replicated(cv, margin, a1, a2)
    sweep[cv] = row
z = sweep
for i, cv in enumerate(cvs):
    cells = []
    for margin in (0.0, THETA):
        r = z[cv][margin]
        d = ('' if i == 0 else
             f'{z[cvs[i-1]][margin]["ncc"] - r["ncc"]:+8.0f}')
        cells.append(f'{r["sens"]*100:>6.1f}%{r["spec"]*100:>6.1f}%'
                     f'{r["ncc"]:>8.0f}{d:>8}')
    print(f'  {cv:>6.2f} | {cells[0]} | {cells[1]}')
print('\n  dNCC is the change in net correct classifications relative to the')
print('  next-higher CV, i.e. the gain from tightening measurement by one step.')

banner('R1. FULL OPERATING CHARACTERISTICS UNDER THE PRIMARY RULE (Table S9)')
print(f'  {"CV":>6}{"sens":>8}{"spec":>8}{"FP rate":>10}{"weak":>7}'
      f'{"S2":>6}{"weeks":>8}')
print('  ' + '-' * 53)
for cv in cvs:
    r = z[cv][THETA]
    print(f'  {cv:>6.2f}{r["sens"]*100:>7.1f}%{r["spec"]*100:>7.1f}%'
          f'{r["fp_rate"]*100:>9.1f}%{r["weak"]*100:>6.0f}%'
          f'{r["stage2_rate"]*100:>5.0f}%{r["mean_weeks"]:>8.1f}')

# ---------------------------------------------------------------------------
# R1.3: NCC weights FP and FN equally -- check other cost ratios
# ---------------------------------------------------------------------------
banner('R1. SENSITIVITY TO THE FALSE-POSITIVE / FALSE-NEGATIVE COST RATIO')
print('  Expected loss per patient; the CV minimizing loss is starred.')
print(f'  {"FP:FN cost":>12} | ' + ''.join(f'{cv:>8.2f}' for cv in cvs)
      + '   best')
for margin, label in ((0.0, 'zero-null'), (THETA, 'matched-null')):
    print(f'  -- {label}')
    for ratio in (0.25, 0.5, 1.0, 2.0, 4.0):
        losses = [weighted_accuracy(
            {'tp': z[cv][margin]['tp'], 'fn': z[cv][margin]['fn'],
             'fp': z[cv][margin]['fp'], 'tn': z[cv][margin]['tn']},
            fp_cost=ratio, fn_cost=1.0) for cv in cvs]
        best = cvs[int(np.argmin(losses))]
        print(f'  {ratio:>10.2f}:1 | ' + ''.join(f'{v:>8.3f}' for v in losses)
              + f'   {best:.2f}')

banner('DONE')
print('  Reproduce with: python3 validation/decision_rule_analysis.py')
