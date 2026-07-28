"""
NON-RESPONDER FRACTION SENSITIVITY (supplementary Table S5).

The reference simulation imposes an 18% non-responder fraction. That value is a
scenario assumption, not an estimate: no published trial reports individual IS
responses in enough detail to estimate it, and the low between-study I2 of the
fiber meta-analysis says nothing about within-study individual heterogeneity.

This script sweeps the imposed fraction through the full adaptive protocol
(Stage 1 -> decision -> Stage 2 for borderline -> dichotomous final call) and
reports how the operating characteristics move.

Within a seed the cohort is generated once and the same uniform draw is
thresholded at each fraction, so the only thing that changes between rows is
which patients are assigned the non-responder effect distribution.

Two conventions are reported, matching the rest of the manuscript:
  - seed 777, the single representative execution behind Table 2
  - mean of 50 single-run replications (seeds 700-749), the Table S1/S3 convention
"""
import numpy as np

N_PAT   = 1000
CV_STD  = 0.15
N_S1    = 6       # measures per arm, stage 1
N_S2    = 3       # additional measures per arm, stage 2
DT_S1   = 1.645 * CV_STD * np.sqrt(2.0 / N_S1)
DT_S2   = 1.645 * CV_STD * np.sqrt(2.0 / (N_S1 + N_S2))

FRACTIONS = [0.10, 0.18, 0.20, 0.30]
REF_FRAC  = 0.18
SEEDS_50  = range(700, 750)


def cohort(seed, frac):
    """Reproduces the reference cohort draw order exactly, including the
    non-responder assignment (see reproduce_manuscript_numbers.py), with the
    imposed non-responder fraction as the only free parameter."""
    rng = np.random.default_rng(seed)
    egfr = np.clip(rng.normal(22, 6, N_PAT), 10, 35)
    gut = np.clip(rng.normal(1.0, 0.35, N_PAT), 0.3, 1.8)
    _slope = np.clip(rng.normal(-2.0, 0.8, N_PAT), -5.0, -0.5)
    _raw = (25.0 / egfr) ** 1.2 * np.clip(rng.normal(1.0, 0.30, N_PAT), 0.3, 2.5)
    bis = np.clip(5.4 + 3.6 * (_raw - _raw.mean()) / _raw.std(), 0.5, None)
    tau = np.clip(rng.normal(0.30 * gut, 0.14), 0, 0.70)
    nr = rng.random(N_PAT) < frac
    tau[nr] = np.clip(rng.normal(0.03, 0.02, nr.sum()), 0, 0.08)
    return bis, tau


def run_once(bis, tau, seed):
    """One pass through the full adaptive protocol for every patient, in the
    same draw order as the reference (A1, B1, A2, B2, then stage 2)."""
    rng = np.random.default_rng(seed)
    cls = np.full(N_PAT, 'N', dtype='U1')
    went_s2 = np.zeros(N_PAT, bool)
    half = N_S1 // 2
    for p in range(N_PAT):
        b, t = bis[p], tau[p]
        A1 = b * (1 + rng.normal(0, CV_STD, half))
        B1 = b * (1 - t) * (1 + rng.normal(0, CV_STD, half))
        A2 = b * (1 + rng.normal(0, CV_STD, half))
        B2 = b * (1 - t) * (1 + rng.normal(0, CV_STD, half))
        A = np.concatenate([A1, A2]); B = np.concatenate([B1, B2])
        obs1 = (A.mean() - B.mean()) / A.mean() if A.mean() > 0 else 0
        if obs1 > DT_S1:
            cls[p] = 'R'
        elif obs1 < 0:
            cls[p] = 'N'
        else:
            went_s2[p] = True
            A3 = b * (1 + rng.normal(0, CV_STD, N_S2))
            B3 = b * (1 - t) * (1 + rng.normal(0, CV_STD, N_S2))
            Ac = np.concatenate([A, A3]); Bc = np.concatenate([B, B3])
            obs2 = (Ac.mean() - Bc.mean()) / Ac.mean() if Ac.mean() > 0 else 0
            cls[p] = 'R' if obs2 > DT_S2 else 'N'
    return cls, went_s2


def metrics(cls, tau, went_s2):
    true_resp = tau >= 0.10
    non_r = ~true_resp
    tp = ((cls == 'R') & true_resp).sum(); fn = ((cls == 'N') & true_resp).sum()
    fp = ((cls == 'R') & non_r).sum();     tn = ((cls == 'N') & non_r).sum()
    weak = true_resp & (tau < 0.20)
    return dict(
        realized_nr=non_r.sum() / N_PAT,
        sens=tp / (tp + fn) if tp + fn else float('nan'),
        spec=tn / (tn + fp) if tn + fp else float('nan'),
        fp_rate=fp / (fp + tn) if fp + tn else float('nan'),
        weak_det=((cls == 'R') & weak).sum() / weak.sum() if weak.sum() else float('nan'),
        s2_rate=went_s2.mean(),
        tp=tp, fn=fn, fp=fp, tn=tn, n_resp=true_resp.sum(), n_nr=non_r.sum())


def evaluate(frac, seed):
    bis, tau = cohort(seed, frac)
    cls, went_s2 = run_once(bis, tau, seed)
    return metrics(cls, tau, went_s2)


print("=" * 86)
print("NON-RESPONDER FRACTION SENSITIVITY (full adaptive protocol, CV = 0.15)")
print("=" * 86)
print(f"  DT6 = {DT_S1*100:.1f}%   DT9 = {DT_S2*100:.1f}%   N = {N_PAT}")

print(f"\n  (A) Single representative execution, seed 777")
print(f"  {'imposed':>8}{'realized':>10}{'resp':>7}{'NR':>6}{'sens':>8}{'spec':>8}"
      f"{'FP':>7}{'weak':>7}{'stage2':>8}")
print("  " + "-" * 74)
rows_777 = {}
for f in FRACTIONS:
    m = evaluate(f, 777); rows_777[f] = m
    tag = "  <- reference" if f == REF_FRAC else ""
    print(f"  {f*100:>7.0f}%{m['realized_nr']*100:>9.1f}%{m['n_resp']:>7}{m['n_nr']:>6}"
          f"{m['sens']*100:>7.1f}%{m['spec']*100:>7.1f}%{m['fp_rate']*100:>6.1f}%"
          f"{m['weak_det']*100:>6.0f}%{m['s2_rate']*100:>7.0f}%{tag}")

print(f"\n  (B) Mean of 50 single-run replications, seeds 700-749")
print(f"  {'imposed':>8}{'realized':>10}{'sens':>8}{'spec':>8}{'FP':>7}{'weak':>7}{'stage2':>8}")
print("  " + "-" * 56)
means = {}
for f in FRACTIONS:
    acc = [evaluate(f, s) for s in SEEDS_50]
    m = {k: float(np.mean([a[k] for a in acc])) for k in
         ('realized_nr', 'sens', 'spec', 'fp_rate', 'weak_det', 's2_rate')}
    means[f] = m
    print(f"  {f*100:>7.0f}%{m['realized_nr']*100:>9.1f}%{m['sens']*100:>7.1f}%"
          f"{m['spec']*100:>7.1f}%{m['fp_rate']*100:>6.1f}%{m['weak_det']*100:>6.0f}%"
          f"{m['s2_rate']*100:>7.0f}%")

ref = means[REF_FRAC]
print(f"\n  (C) Change from the {REF_FRAC*100:.0f}% reference (50-seed means, percentage points)")
print(f"  {'imposed':>8}{'d sens':>9}{'d spec':>9}{'d FP':>8}{'d weak':>9}")
print("  " + "-" * 43)
for f in FRACTIONS:
    m = means[f]
    print(f"  {f*100:>7.0f}%{(m['sens']-ref['sens'])*100:>+8.1f}{(m['spec']-ref['spec'])*100:>+9.1f}"
          f"{(m['fp_rate']-ref['fp_rate'])*100:>+8.1f}{(m['weak_det']-ref['weak_det'])*100:>+9.1f}")

sens_range = (min(m['sens'] for m in means.values()), max(m['sens'] for m in means.values()))
spec_range = (min(m['spec'] for m in means.values()), max(m['spec'] for m in means.values()))
print(f"""
  INTERPRETATION
  Imposing 10% to 30% non-responders moves the realized non-responder share
  (tau < 10%, which also collects low-effect draws from the responder
  distribution) across {min(m['realized_nr'] for m in means.values())*100:.0f}% to {max(m['realized_nr'] for m in means.values())*100:.0f}% of the cohort. Over that range
  sensitivity stays within {sens_range[0]*100:.1f}-{sens_range[1]*100:.1f}% and specificity within
  {spec_range[0]*100:.1f}-{spec_range[1]*100:.1f}%. Classification performance is therefore driven by where
  each patient's effect sits relative to the decision threshold, not by how
  many non-responders the generating model contains.
""")
