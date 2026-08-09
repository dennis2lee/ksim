"""
Design efficiency on equal terms (review item R5).

The manuscript compared "lower the CV from 0.25 to 0.15" against "double the
number of cycles" and concluded the first was more efficient because it saved
24 weeks. That comparison was not symmetric. Extending the protocol was priced
at its true cost in weeks, while the CV reduction was simply granted: no
mechanism, no success probability, and a cost table that counted only the
duplicate assay.

This script prices both levers in the same currency. The reachable CV values
come from the variance-component analysis (R4) rather than being assumed, so a
standardization package buys a specific CV rather than an arbitrary one, and
CV and measurement count are optimized jointly rather than treated as
independent levers.

Burden is reported as three quantities a trialist actually trades off:
protocol duration in weeks, number of venipunctures, and number of assays.
Patient-facing standardization burden is flagged separately because it is not
commensurable with the others.

Outputs feed manuscript Table 6 and Table S13.
"""

import numpy as np
from nof1_core import (make_cohort, run_protocol, critical_value, total_cv,
                       calibrate_boundaries)

N_PAT = 1000
COHORT_SEED = 777
THETA = 0.10
REP_SEEDS = range(700, 720)

cohort = make_cohort(N_PAT, COHORT_SEED, 0.18)
_cal = calibrate_boundaries(cohort, cv=0.15, null_margin=THETA,
                            alpha_total=0.05, seeds=range(900, 920),
                            order='randomized')
PRIMARY = dict(theta=THETA, null_margin=THETA,
               alpha1=_cal['alpha1'], alpha2=_cal['alpha2'], order='randomized')

# Standardization packages, with the CV each one reaches under the intermediate
# variance-component scenario of variance_components_analysis.py
# (cv_bio_irreducible 0.254, cv_bio_standardizable 0.254, cv_pre 0.08,
#  cv_analytical 0.06).
PACKAGES = [
    ('none', total_cv(cv_bio=np.hypot(0.254, 0.254), cv_pre=0.08,
                      cv_analytical=0.06, n_replicates=1), 0, 1),
    ('sampling standardization', total_cv(cv_bio=0.254, cv_pre=0.04,
                                          cv_analytical=0.06, n_replicates=1),
     1, 1),
    ('standardization + duplicate assay',
     total_cv(cv_bio=0.254, cv_pre=0.04, cv_analytical=0.06, n_replicates=2),
     1, 2),
    ('hypothetical CV 0.15 (not reachable above)', 0.15, 2, 2),
]

DESIGNS = [(2, 3), (3, 3), (4, 3), (2, 2), (2, 4), (3, 4)]


def banner(t):
    print('\n' + '=' * 96)
    print(t)
    print('=' * 96)


banner('R5. DECISION THRESHOLD BY DESIGN AND ACHIEVABLE CV')
print('  Critical value for H0: tau <= 0.10 at one-sided alpha 0.05.')
print('  Weeks and venipunctures are Stage 1 only; Stage 2 adds 12 weeks and')
print('  one more cycle of draws for the patients who enter it.\n')
print(f'  {"package":<42}{"CV":>6}{"design":>9}{"n/arm":>7}'
      f'{"crit":>8}{"weeks":>7}{"venip":>7}{"assays":>8}')
print('  ' + '-' * 94)
for pkg, cv, std_burden, reps in PACKAGES:
    for cycles, mpp in DESIGNS:
        n_arm = cycles * mpp
        weeks = cycles * 12
        venip = cycles * 2 * mpp
        crit = critical_value(cv, n_arm, 0.05, THETA)
        print(f'  {pkg:<42}{cv:>6.3f}{f"{cycles}x{mpp}":>9}{n_arm:>7}'
              f'{crit*100:>7.1f}%{weeks:>7}{venip:>7}{venip*reps:>8}')

banner('R5. OPERATING CHARACTERISTICS ON THE SAME GRID (primary rule)')
print(f'  {"package":<42}{"CV":>6}{"design":>9}{"sens":>8}{"spec":>8}'
      f'{"weak":>7}{"weeks":>7}{"venip":>7}')
print('  ' + '-' * 94)
grid = {}
for pkg, cv, std_burden, reps in PACKAGES:
    for cycles, mpp in DESIGNS:
        runs = [run_protocol(cohort, cv=cv, measurement_seed=s,
                             n_cycles=cycles, meas_per_period=mpp, **PRIMARY)
                for s in REP_SEEDS]
        sens = float(np.mean([r['sens'] for r in runs]))
        spec = float(np.mean([r['spec'] for r in runs]))
        weak = float(np.mean([r['detect']['weak'] for r in runs]))
        s2 = float(np.mean([r['stage2_rate'] for r in runs]))
        weeks = cycles * 12 + 12 * s2
        venip = cycles * 2 * mpp + 2 * mpp * s2
        grid[(pkg, cycles, mpp)] = dict(cv=cv, sens=sens, spec=spec, weak=weak,
                                        weeks=weeks, venip=venip,
                                        assays=venip * reps,
                                        std_burden=std_burden)
        print(f'  {pkg:<42}{cv:>6.3f}{f"{cycles}x{mpp}":>9}{sens*100:>7.1f}%'
              f'{spec*100:>7.1f}%{weak*100:>6.0f}%{weeks:>7.1f}{venip:>7.1f}')

banner('R5. EFFICIENT FRONTIER: SENSITIVITY VERSUS PROTOCOL DURATION')
print('  A design is on the frontier if no other design achieves at least the')
print('  same sensitivity in no more weeks and no more venipunctures.\n')
items = [(k, v) for k, v in grid.items()
         if not k[0].startswith('hypothetical')]
frontier = []
for k, v in items:
    dominated = any(
        (w['sens'] >= v['sens'] + 1e-9 and w['weeks'] <= v['weeks']
         and w['venip'] <= v['venip'] and w['std_burden'] <= v['std_burden'])
        for k2, w in items if k2 != k)
    if not dominated:
        frontier.append((k, v))
frontier.sort(key=lambda kv: kv[1]['weeks'])
print(f'  {"package":<42}{"design":>9}{"sens":>8}{"weeks":>8}{"venip":>8}')
print('  ' + '-' * 76)
for (pkg, cycles, mpp), v in frontier:
    print(f'  {pkg:<42}{f"{cycles}x{mpp}":>9}{v["sens"]*100:>7.1f}%'
          f'{v["weeks"]:>8.1f}{v["venip"]:>8.1f}')

banner('R5. THE LEVERS ARE NOT INDEPENDENT')
base = grid[('none', 2, 3)]
print(f'  Starting point: no standardization, 2x3 design, CV {base["cv"]:.3f}, '
      f'sensitivity {base["sens"]*100:.1f}%\n')
for label, key in (('add sampling standardization, keep 2x3',
                    ('sampling standardization', 2, 3)),
                   ('keep CV, go to 4 cycles',
                    ('none', 4, 3)),
                   ('keep CV, 4 measurements per period',
                    ('none', 2, 4)),
                   ('standardization AND 3 cycles',
                    ('sampling standardization', 3, 3)),
                   ('standardization AND 4 measurements per period',
                    ('sampling standardization', 2, 4))):
    v = grid[key]
    print(f'  {label:<48}sens {v["sens"]*100:>5.1f}% '
          f'({(v["sens"]-base["sens"])*100:+5.1f} pp), '
          f'{v["weeks"]:>5.1f} wk, {v["venip"]:>5.1f} venip')
gain_std = grid[('sampling standardization', 2, 3)]['sens'] - base['sens']
gain_cyc = grid[('none', 4, 3)]['sens'] - base['sens']
gain_both = grid[('sampling standardization', 3, 3)]['sens'] - base['sens']
print(f'\n  Standardization alone: {gain_std*100:+.1f} pp. '
      f'Doubling cycles alone: {gain_cyc*100:+.1f} pp.')
print(f'  Both (standardization + one extra cycle): {gain_both*100:+.1f} pp, '
      f'against {(gain_std+gain_cyc)*100:+.1f} pp if the two gains added.')
print('  The levers act on the same standard error, so their effects are')
print('  sub-additive; treating them as independent double counts the benefit.')

print('\n  Reproduce with: python3 validation/efficiency_analysis.py')
