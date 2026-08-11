"""
Independent validation of the sequential boundaries (review item R12).

Two problems with the previous version:

  1. The boundaries were calibrated on measurement seeds 900-919 and the
     achieved error rate was then measured on seeds 900-949. The calibration
     set was wholly contained in the validation set, so the reported 5.1% was
     partly a description of the data the boundaries were tuned on.

  2. No Monte Carlo uncertainty was reported, so there was no way to judge
     whether 5.1% differed from 5%.

Here the seed blocks are disjoint by construction: calibration uses 900-919 and
validation uses 3000-3199, and the script asserts they do not intersect. Every
rate is reported with a Monte Carlo standard error computed across seeds, each
seed being one independent replication of the whole cohort.

Type I error is evaluated on a null cohort, every tau set to the null margin,
which is the boundary of H0. That is the quantity the boundaries are supposed
to control. It is not the same as the classification specificity reported in
the robustness table, which has true non-responders in its denominator and
mixes patients whose true effect lies anywhere below the margin.
"""

import numpy as np
from nof1_core import (make_cohort, run_protocol, calibrate_boundaries,
                       critical_value, obf_alphas, sd_inflation)

N_PAT = 1000
COHORT_SEED = 777
THETA = 0.10
CV = 0.15

CALIB_SEEDS = range(900, 920)          # boundaries are tuned here
VALID_SEEDS = range(3000, 3200)        # and tested here, on 200 fresh seeds
assert not (set(CALIB_SEEDS) & set(VALID_SEEDS)), 'seed blocks must be disjoint'

cohort = make_cohort(N_PAT, COHORT_SEED, 0.18)


def null_cohort(margin):
    c = make_cohort(N_PAT, COHORT_SEED, 0.18)
    c.tau = np.full(N_PAT, margin)
    return c


def rate_with_se(margin, alpha1, alpha2, seeds, **kw):
    """Overall one-sided type I error and its Monte Carlo standard error.

    Each seed is one independent replication of the entire cohort, so the
    across-seed standard deviation divided by sqrt(number of seeds) is the
    Monte Carlo standard error of the mean rate.
    """
    nc = null_cohort(margin)
    rates = []
    for s in seeds:
        r = run_protocol(nc, cv=CV, measurement_seed=s, theta=margin + 1e-9,
                         null_margin=margin, alpha1=alpha1, alpha2=alpha2,
                         order='randomized', **kw)
        an = ~r['dropped']
        rates.append(((r['cls'] == 'R') & an).sum() / max(an.sum(), 1))
    rates = np.asarray(rates)
    m = float(rates.mean())
    se = float(rates.std(ddof=1) / np.sqrt(len(rates)))
    return m, se, (m - 1.96 * se, m + 1.96 * se)


def oc_with_se(alpha1, alpha2, seeds, **kw):
    """Sensitivity and specificity with Monte Carlo standard errors."""
    out = {}
    runs = [run_protocol(cohort, cv=CV, measurement_seed=s, theta=THETA,
                         null_margin=THETA, alpha1=alpha1, alpha2=alpha2,
                         order='randomized', **kw) for s in seeds]
    for k in ('sens', 'spec', 'stage2_rate'):
        a = np.array([r[k] for r in runs], float)
        out[k] = (float(a.mean()), float(a.std(ddof=1) / np.sqrt(len(a))))
    return out


def banner(t):
    print('\n' + '=' * 82)
    print(t)
    print('=' * 82)


banner('R12. CALIBRATION AND VALIDATION USE DISJOINT SEED BLOCKS')
print(f'  calibration seeds {CALIB_SEEDS.start}-{CALIB_SEEDS.stop - 1} '
      f'({len(CALIB_SEEDS)} replications)')
print(f'  validation  seeds {VALID_SEEDS.start}-{VALID_SEEDS.stop - 1} '
      f'({len(VALID_SEEDS)} replications)')
print(f'  intersection: {len(set(CALIB_SEEDS) & set(VALID_SEEDS))}')

cal = calibrate_boundaries(cohort, cv=CV, null_margin=THETA, alpha_total=0.05,
                           seeds=CALIB_SEEDS, order='randomized')
nom1, nom2 = cal['nominal']
print(f'\n  O\'Brien-Fleming nominal levels: alpha1 {nom1:.4f}, alpha2 {nom2:.4f}')
print(f'  calibrated:                     alpha1 {cal["alpha1"]:.4f}, '
      f'alpha2 {cal["alpha2"]:.4f}  (common scale {cal["scale"]:.3f})')
print(f'  achieved on the CALIBRATION set: {cal["achieved"]*100:.2f}%  '
      '(in-sample, reported for comparison only)')

banner('R12. OVERALL TYPE I ERROR ON THE HELD-OUT VALIDATION SEEDS')
print('  Null cohort: every tau set to the margin, so this is the error rate at')
print('  the boundary of H0. Intervals are 1.96 x Monte Carlo standard error.\n')
print(f'  {"rule":<44}{"error":>9}{"MC SE":>9}   95% interval')
print('  ' + '-' * 80)
rows = [
    ('zero-null, alpha 0.05 per look', 0.0, 0.05, 0.05),
    ('zero-null, calibrated', 0.0, 0.0241, 0.0547),
    ('matched-null, alpha 0.05 per look', THETA, 0.05, 0.05),
    ('matched-null, calibrated (primary)', THETA, cal['alpha1'], cal['alpha2']),
]
for label, margin, a1, a2 in rows:
    m, se, ci = rate_with_se(margin, a1, a2, VALID_SEEDS)
    print(f'  {label:<44}{m*100:>8.2f}%{se*100:>8.2f}%   '
          f'[{ci[0]*100:.2f}%, {ci[1]*100:.2f}%]')
print('\n  The primary rule\'s interval is the one to read against the 5% target.')

banner('R12. TYPE I ERROR UNDER DEPARTURES FROM THE REFERENCE MODEL')
print('  Same calibrated boundaries, held-out seeds, one departure at a time.')
print('  The boundary uses the realized number of measurements in each arm. An')
print('  earlier version used the nominal design count, which made the boundary')
print('  too low whenever measurements went missing and pushed the error rate to')
print('  6.16% at 10% missing and 7.57% at 20% missing.\n')
print(f'  {"departure":<44}{"error":>9}{"MC SE":>9}   95% interval')
print('  ' + '-' * 80)
DEPARTURES = [
    ('reference (independent errors)', {}),
    ('within-period correlation, rho = 0.3', dict(rho=0.3)),
    ('within-period correlation, rho = 0.5', dict(rho=0.5)),
    ('10% of measurements missing', dict(p_missing=0.10)),
    ('20% of measurements missing', dict(p_missing=0.20)),
    ('20% dropout', dict(p_dropout=0.20)),
    ('onset ramp, time constant 2 weeks', dict(onset_weeks=2.0)),
    ('washout 4 weeks', dict(washout_weeks=4)),
]
for label, kw in DEPARTURES:
    m, se, ci = rate_with_se(THETA, cal['alpha1'], cal['alpha2'],
                             VALID_SEEDS, **kw)
    flag = '  <-- exceeds 5%' if ci[0] > 0.05 else ''
    print(f'  {label:<44}{m*100:>8.2f}%{se*100:>8.2f}%   '
          f'[{ci[0]*100:.2f}%, {ci[1]*100:.2f}%]{flag}')

banner('R12. OPERATING CHARACTERISTICS WITH MONTE CARLO UNCERTAINTY')
print('  Reference cohort, held-out seeds, primary rule.\n')
oc = oc_with_se(cal['alpha1'], cal['alpha2'], VALID_SEEDS)
for k, name in (('sens', 'sensitivity'), ('spec', 'specificity'),
                ('stage2_rate', 'Stage 2 entry')):
    m, se = oc[k]
    print(f'  {name:<20}{m*100:>7.2f}%  MC SE {se*100:.2f}%  '
          f'95% interval [{(m-1.96*se)*100:.2f}%, {(m+1.96*se)*100:.2f}%]')

banner('R12. WHY THE CALIBRATED ALPHA2 EXCEEDS 0.05')
n1, n2 = 6, 9
print(f'  Stage 1 nominal {nom1:.4f} -> calibrated {cal["alpha1"]:.4f}')
print(f'  Stage 2 nominal {nom2:.4f} -> calibrated {cal["alpha2"]:.4f}')
print(f'\n  Stage 1 critical value: '
      f'{critical_value(CV, n1, cal["alpha1"], THETA)*100:.1f}%')
print(f'  Stage 2 critical value: '
      f'{critical_value(CV, n2, cal["alpha2"], THETA)*100:.1f}%')
print("""
  alpha2 is not a per-test level applied to an unselected patient. Only the
  patients whose Stage 1 estimate fell between the null margin and the Stage 1
  critical value reach the second look; those above were already called
  responders and those below were already called non-responders. alpha2 is the
  conditional level applied to that selected middle group, and the quantity
  constrained to 5% is the error of the whole procedure, not either look. A
  conditional level above 0.05 is therefore expected rather than anomalous, and
  the number to compare against 5% is the validated overall rate above.""")

print('\n  Reproduce with: python3 validation/sequential_error_validation.py')
