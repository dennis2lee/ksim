"""
Serial correlation as a decision-rule problem (review item R13).

The previous version treated within-period correlation as a limitation. It is
better understood as a defect in the boundary itself: the analytic critical
value is derived under independent errors, and under AR(1) correlation the
standard error of an arm mean is larger than that derivation assumes, so the
boundary sits too low and the procedure over-rejects.

Two quantities were previously reported in a way that invited conflation:

  classification specificity  denominator: simulated non-responders, that is
                              every patient with tau below the margin. Reported
                              as 94.7% at rho = 0.3.

  boundary type I error       denominator: patients whose tau equals the margin
                              exactly. Reported as 7.22% at rho = 0.3.

They answer different questions. Specificity stays high because most simulated
non-responders have tau well below the margin and are easy to classify; the
error rate at the boundary is what the alpha level is supposed to control, and
that is what correlation damages.

This script implements a correlation-adjusted boundary using the AR(1) variance
inflation factor and compares three strategies: the independent boundary, a
boundary matched to the true rho, and a conservative boundary assuming a
plausible worst case when rho is unknown.
"""

import numpy as np
from nof1_core import (make_cohort, run_protocol, calibrate_boundaries,
                       critical_value, sd_inflation)

N_PAT, COHORT_SEED, THETA, CV = 1000, 777, 0.10, 0.15
K = 3                                   # measurements per period, the AR(1) block
SEEDS = range(3000, 3100)               # held out from calibration
RHOS = [0.0, 0.2, 0.3, 0.4, 0.5]
RHO_CONSERVATIVE = 0.5                  # assumed when rho is unknown

cohort = make_cohort(N_PAT, COHORT_SEED, 0.18)
cal = calibrate_boundaries(cohort, cv=CV, null_margin=THETA, alpha_total=0.05,
                           seeds=range(900, 920), order='randomized')
A1, A2 = cal['alpha1'], cal['alpha2']


def null_cohort():
    c = make_cohort(N_PAT, COHORT_SEED, 0.18)
    c.tau = np.full(N_PAT, THETA)
    return c


NULL = null_cohort()


def boundary_error(rho_true, rho_assumed):
    """Type I error at the margin, with MC standard error."""
    rates = []
    for s in SEEDS:
        r = run_protocol(NULL, cv=CV, measurement_seed=s, theta=THETA + 1e-9,
                         null_margin=THETA, alpha1=A1, alpha2=A2,
                         order='randomized', rho=rho_true,
                         rho_assumed=rho_assumed)
        rates.append((r['cls'] == 'R').mean())
    a = np.asarray(rates)
    return float(a.mean()), float(a.std(ddof=1) / np.sqrt(len(a)))


def oc(rho_true, rho_assumed):
    runs = [run_protocol(cohort, cv=CV, measurement_seed=s, theta=THETA,
                         null_margin=THETA, alpha1=A1, alpha2=A2,
                         order='randomized', rho=rho_true,
                         rho_assumed=rho_assumed) for s in SEEDS]
    return (float(np.mean([r['sens'] for r in runs])),
            float(np.mean([r['spec'] for r in runs])))


def banner(t):
    print('\n' + '=' * 86)
    print(t)
    print('=' * 86)


banner('R13. THE VARIANCE INFLATION THE INDEPENDENT BOUNDARY IGNORES')
print(f'  AR(1) within blocks of k = {K}; blocks independent, so the factor does')
print('  not depend on how many blocks an arm contains.\n')
print(f'  {"rho":>6}{"SD inflation":>15}{"crit, independent":>20}'
      f'{"crit, adjusted":>17}')
print('  ' + '-' * 60)
for rho in RHOS:
    print(f'  {rho:>6.1f}{sd_inflation(rho, K):>15.4f}'
          f'{critical_value(CV, 6, 0.05, THETA)*100:>19.1f}%'
          f'{critical_value(CV, 6, 0.05, THETA, rho, K)*100:>16.1f}%')

banner('R13. TWO DIFFERENT QUANTITIES, STATED SEPARATELY')
print('  Independent boundary throughout; only the true rho changes.\n')
print(f'  {"true rho":>9}{"boundary type I error":>26}{"classification specificity":>29}')
print(f'  {"":>9}{"(tau = margin exactly)":>26}{"(all tau below margin)":>29}')
print('  ' + '-' * 66)
for rho in RHOS:
    e, se = boundary_error(rho, 0.0)
    _, sp = oc(rho, 0.0)
    print(f'  {rho:>9.1f}{e*100:>19.2f}% ± {se*100:.2f}{sp*100:>27.1f}%')
print('\n  Specificity moves by a few points while the error rate at the boundary')
print('  more than doubles. Reporting only the first hides the failure.')

banner('R13. THREE BOUNDARY STRATEGIES')
print('  matched      : the boundary uses the true rho (an upper bound on what')
print('                 a well-informed trial could do)')
print('  independent  : the boundary assumes rho = 0 (the previous version)')
print(f'  conservative : the boundary always assumes rho = {RHO_CONSERVATIVE},')
print('                 which is what a trial without a correlation estimate')
print('                 can prespecify\n')
print(f'  {"true rho":>9} | {"independent":^22} | {"matched":^22} | '
      f'{"conservative":^22}')
print(f'  {"":>9} | {"error":>9}{"sens":>13} | {"error":>9}{"sens":>13} | '
      f'{"error":>9}{"sens":>13}')
print('  ' + '-' * 84)
summary = {}
for rho in RHOS:
    cells = []
    for name, assumed in (('independent', 0.0), ('matched', rho),
                          ('conservative', RHO_CONSERVATIVE)):
        e, se = boundary_error(rho, assumed)
        sn, sp = oc(rho, assumed)
        summary[(rho, name)] = (e, se, sn, sp)
        cells.append(f'{e*100:>8.2f}%{sn*100:>12.1f}%')
    print(f'  {rho:>9.1f} | {cells[0]} | {cells[1]} | {cells[2]}')

banner('SUMMARY')
ind_max = max(summary[(r, 'independent')][0] for r in RHOS)
mat_max = max(summary[(r, 'matched')][0] for r in RHOS)
con_max = max(summary[(r, 'conservative')][0] for r in RHOS)
con_sens0 = summary[(0.0, 'conservative')][2]
ind_sens0 = summary[(0.0, 'independent')][2]
print(f'  Worst-case type I error across rho 0 to 0.5:')
print(f'    independent boundary  {ind_max*100:.2f}%')
print(f'    matched boundary      {mat_max*100:.2f}%')
print(f'    conservative boundary {con_max*100:.2f}%')
print(f'\n  Cost of the conservative boundary when errors are in fact')
print(f'  independent: sensitivity {ind_sens0*100:.1f}% -> {con_sens0*100:.1f}% '
      f'({(con_sens0-ind_sens0)*100:+.1f} points).')
print("""
  A trial that cannot estimate the within-period correlation in advance should
  prespecify the conservative boundary and accept that loss, or plan an interim
  recalibration once enough within-period pairs have been observed to estimate
  rho. Reporting a 5% overall type I error without naming the independent-error
  reference model overstates what the procedure guarantees.""")

print('\n  Reproduce with: python3 validation/correlation_boundary_analysis.py')
