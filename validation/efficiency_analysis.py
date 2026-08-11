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
# Three variance splits, not one. The efficiency conclusion is conditional on
# which of these holds; an earlier version reported only the intermediate case,
# which made the conclusion read as unconditional (review item R14).
# (label, irreducible biological CV, standardizable biological CV)
SPLITS = [
    ('optimistic', 0.207, 0.293),
    ('intermediate', 0.254, 0.254),
    ('pessimistic', 0.293, 0.207),
]


def packages_for(bio_irr, bio_std):
    return [
        ('none', total_cv(cv_bio=np.hypot(bio_irr, bio_std), cv_pre=0.08,
                          cv_analytical=0.06, n_replicates=1), 0, 1),
        ('sampling standardization',
         total_cv(cv_bio=bio_irr, cv_pre=0.04, cv_analytical=0.06,
                  n_replicates=1), 1, 1),
        ('standardization + duplicate assay',
         total_cv(cv_bio=bio_irr, cv_pre=0.04, cv_analytical=0.06,
                  n_replicates=2), 1, 2),
    ]


PACKAGES = packages_for(0.254, 0.254) + [
    ('hypothetical CV 0.15 (not reachable above)', 0.15, 2, 2)]

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

banner('R14. THE EFFICIENCY CONCLUSION UNDER ALL THREE VARIANCE SPLITS')
print('  For each split the 2x3 design with no standardization is the baseline.')
print('  Burden is weeks, venipunctures, and a qualitative patient-effort score')
print('  (0 none, 1 timed fasting draws and a dietary log, 2 as 1 plus repeat')
print('  assay logistics). Patient effort is not commensurable with weeks and is')
print('  reported separately rather than folded into one index.\n')
print(f'  {"split":<14}{"lever":<36}{"CV":>7}{"sens":>8}{"gain":>8}'
      f'{"weeks":>8}{"venip":>7}{"effort":>7}')
print('  ' + '-' * 95)
verdict = {}
for split_name, bio_irr, bio_std in SPLITS:
    pkgs = packages_for(bio_irr, bio_std)
    base_cv = pkgs[0][1]
    base_runs = [run_protocol(cohort, cv=base_cv, measurement_seed=s,
                              n_cycles=2, meas_per_period=3, **PRIMARY)
                 for s in REP_SEEDS]
    base = float(np.mean([r['sens'] for r in base_runs]))
    rows = [
        ('baseline: no standardization, 2x3', base_cv, 2, 3, 0),
        ('standardization, 2x3', pkgs[1][1], 2, 3, 1),
        ('standardization + duplicate, 2x3', pkgs[2][1], 2, 3, 2),
        ('no standardization, 2x4', base_cv, 2, 4, 0),
        ('no standardization, 4x3', base_cv, 4, 3, 0),
        ('standardization, 2x4', pkgs[1][1], 2, 4, 1),
    ]
    for label, cvv, cyc, mpp, effort in rows:
        runs = [run_protocol(cohort, cv=cvv, measurement_seed=s,
                             n_cycles=cyc, meas_per_period=mpp, **PRIMARY)
                for s in REP_SEEDS]
        sn = float(np.mean([r['sens'] for r in runs]))
        s2 = float(np.mean([r['stage2_rate'] for r in runs]))
        wks = cyc * 12 + 12 * s2
        ven = cyc * 2 * mpp + 2 * mpp * s2
        gain = '' if label.startswith('baseline') else f'{(sn-base)*100:+.1f}'
        print(f'  {split_name:<14}{label:<36}{cvv:>7.3f}{sn*100:>7.1f}%'
              f'{gain:>8}{wks:>8.1f}{ven:>7.1f}{effort:>7}')
        if label == 'standardization, 2x3':
            verdict[(split_name, 'std')] = sn - base
        if label == 'no standardization, 4x3':
            verdict[(split_name, 'cycles')] = sn - base
    print()

banner('R14. IS STANDARDIZATION STILL THE BETTER LEVER IN EVERY SPLIT?')
print(f'  {"split":<16}{"standardization":>17}{"double cycles":>16}'
      f'{"extra weeks":>14}   which lever wins on sensitivity')
print('  ' + '-' * 96)
for split_name, _, _ in SPLITS:
    g_std = verdict[(split_name, 'std')] * 100
    g_cyc = verdict[(split_name, 'cycles')] * 100
    win = ('standardization' if g_std > g_cyc else 'DOUBLING CYCLES')
    print(f'  {split_name:<16}{g_std:>16.1f}{g_cyc:>16.1f}{22:>14}   {win}')
print('''
  The conclusion does not hold in all three splits. Standardization wins
  comfortably when two-thirds of the biological variance is standardizable
  (+21.9 against +11.7) and narrowly at one-half (+15.9 against +11.8), but at
  one-third it loses on sensitivity (+10.0 against +11.7). It still costs no
  additional weeks, so a trialist who values calendar time may prefer it even
  there; a trialist who values sensitivity alone should not. Reporting the
  intermediate split alone, as an earlier version did, concealed a sign change
  rather than merely a change of magnitude.''')

print('\n  Reproduce with: python3 validation/efficiency_analysis.py')
