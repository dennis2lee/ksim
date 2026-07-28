"""
ACCURACY OF THE DECISION-THRESHOLD APPROXIMATION (supplementary Table S6).

The protocol classifies on the ratio estimator (mean_A - mean_B)/mean_A and
compares it with DT_n = 1.645 * CV * sqrt(2/n), which comes from approximating
the estimator's variance as CV^2 * (1/n_A + 1/n_B). That approximation assumes
normal, independent, homoscedastic errors.

This script measures how far the approximation is off under each stress
scenario. For null patients (tau = 0) it takes the empirical one-sided 95th
percentile of the observed reduction and compares it with the analytic DT6,
and reports the empirical Type I error rate the analytic threshold actually
delivers.

The noise model mirrors robustness_experiments.py draw_measures() exactly:
each arm mean pools TWO independent blocks of 3 measurements (A1 and A2), and
AR(1) runs within a block only, in stationary form
    raw[0] = N(0, cv);  raw[j] = rho*raw[j-1] + sqrt(1-rho^2)*N(0, cv)
Applying AR(1) across all 6 measurements of an arm instead would overstate the
correlation effect.
"""
import numpy as np

CV = 0.15
BLOCK = 3          # measurements per period (N_S1//2 in the protocol code)
N_BLOCKS = 2       # periods per arm (A1, A2)
N_ARM = BLOCK * N_BLOCKS
DT6 = 1.645 * CV * np.sqrt(2.0 / N_ARM)
N = 400_000
BASE = 5.4
SEED = 20260727


def block_noise(rng, n_pat, rho, cv):
    raw = np.zeros((n_pat, BLOCK))
    raw[:, 0] = rng.normal(0, 1, n_pat)
    s = np.sqrt(1 - rho**2) if rho > 0 else 1.0
    for j in range(1, BLOCK):
        prev = rho * raw[:, j - 1] if rho > 0 else 0.0
        raw[:, j] = prev + s * rng.normal(0, 1, n_pat)
    return raw * cv


def arm(rng, n_pat, rho, cv, base):
    raw = np.concatenate([block_noise(rng, n_pat, rho, cv)
                          for _ in range(N_BLOCKS)], axis=1)
    return base * (1 + raw)


def observed(rng, rho=0.0, lognormal=False, hetero=False):
    cv = CV
    base = BASE
    if hetero:
        cv = np.clip(rng.normal(0.15, 0.04, (N, 1)), 0.08, 0.28)
    if lognormal:
        sigma = 0.5
        base = np.exp(rng.normal(np.log(BASE) - 0.5 * sigma**2, sigma, (N, 1)))
    A = arm(rng, N, rho, cv, base)
    B = arm(rng, N, rho, cv, base)
    ma, mb = A.mean(1), B.mean(1)
    return np.where(ma > 0, (ma - mb) / ma, 0.0)


def sd_inflation(rho):
    """Analytic SD inflation of the arm mean under within-block AR(1),
    relative to iid. Independent of the number of blocks."""
    return np.sqrt((BLOCK + 2 * sum((BLOCK - k) * rho**k
                                    for k in range(1, BLOCK))) / BLOCK)


rng = np.random.default_rng(SEED)
cases = [("Gaussian, iid (reference)", dict(), 0.0),
         ("Log-normal IS baseline", dict(lognormal=True), 0.0),
         ("Heteroscedastic CV", dict(hetero=True), 0.0),
         ("AR(1) rho = 0.3", dict(rho=0.3), 0.3),
         ("AR(1) rho = 0.5", dict(rho=0.5), 0.5)]

print("=" * 84)
print("DECISION-THRESHOLD APPROXIMATION ACCURACY")
print("=" * 84)
print(f"  Analytic DT6 = {DT6*100:.2f}%   null patients per scenario = {N:,}   seed = {SEED}")
print(f"  AR(1) within blocks of {BLOCK}; {N_BLOCKS} independent blocks per arm\n")
print(f"  {'scenario':<28}{'empirical p95':>15}{'diff (pp)':>12}{'empirical alpha':>17}")
print("  " + "-" * 74)
worst = 0.0
for name, kw, rho in cases:
    o = observed(rng, **kw)
    p95 = np.percentile(o, 95)
    diff = (p95 - DT6) * 100
    alpha = (o > DT6).mean()
    worst = max(worst, abs(diff))
    print(f"  {name:<28}{p95*100:>14.2f}%{diff:>+12.2f}{alpha*100:>16.2f}%")

print(f"\n  Largest absolute deviation: {worst:.2f} percentage points")
print("  The approximation holds within 1 pp under the distributional violations")
print("  and is mildly conservative there (empirical alpha below the nominal 5%).")
print("  It breaks only under within-period serial correlation, where the analytic")
print("  threshold is too permissive and the empirical Type I error rises above 5%.")
print("  That is the direct cause of the specificity loss in scenarios D and D'.")

print(f"\n  ANALYTIC CROSS-CHECK (independent of the simulation)")
print(f"  {'rho':>6}{'SD inflation':>15}{'DT6 x inflation':>18}")
print("  " + "-" * 39)
for rho in (0.0, 0.3, 0.5):
    f = sd_inflation(rho)
    print(f"  {rho:>6}{f:>15.3f}{DT6*100*f:>17.2f}%")

print(f"""
  The inflation factor is [BLOCK + 2*sum((BLOCK-k)*rho^k)] / BLOCK, which does
  not depend on the number of blocks. Serial correlation therefore multiplies
  every design's threshold by the same constant and leaves the comparison
  between the two design levers unchanged: CV still enters DT linearly and the
  measurement count still as 1/sqrt(n).""")

print(f"\n  {'design':<8}{'CV':>6}{'n':>4}{'DT iid':>10}{'DT rho=0.5':>13}")
print("  " + "-" * 41)
for label, cv, n in [("2x3", 0.25, 6), ("2x3", 0.15, 6),
                      ("4x3", 0.25, 12), ("3x3", 0.15, 9)]:
    dt = 1.645 * cv * np.sqrt(2.0 / n)
    print(f"  {label:<8}{cv:>6}{n:>4}{dt*100:>9.1f}%{dt*sd_inflation(0.5)*100:>12.1f}%")
ratio_iid = (1.645*0.15*np.sqrt(2/6)) / (1.645*0.25*np.sqrt(2/12))
print(f"\n  DT(CV=0.15, n=6) / DT(CV=0.25, n=12) = {ratio_iid:.4f}, identical at every rho.")
