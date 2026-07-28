# ksim — In-Silico N-of-1 Trial Design for Uremic Toxin Clearance in CKD

Simulation framework for designing and validating personalized n-of-1 crossover
trials that detect individual indoxyl sulfate (IS) reduction from gut-clearance
interventions in chronic kidney disease (CKD stage 3b–4).

**Paper**: Lee P, Lee TJ. "Measurement Noise Optimization as a Design Lever for N-of-1
Trials: In-Silico Validation of Personalized Gut-Clearance Protocols in Chronic
Kidney Disease." (2026, targeting CPT: Pharmacometrics & Systems Pharmacology)

## The Problem

Gut-targeted interventions (dietary fiber, probiotics, oral sorbents) reduce
indoxyl sulfate in CKD — a 2025 meta-analysis of 11 RCTs confirms a pooled
effect (SMD −0.34). But these trials report only population averages. A
clinician cannot tell whether a specific patient responds or not. The
within-person biological variability of IS (CV 35.9%, Pretorius 2013) makes single pre-post
comparisons unreliable.

## What This Code Does

1. **Generates virtual CKD cohorts** from published clinical data (IS levels,
   eGFR distributions, treatment effect sizes, biological variability)
2. **Simulates an n-of-1 crossover protocol** (2-cycle AB design with adaptive
   second stage for borderline patients)
3. **Quantifies the classification-efficiency frontier**: the relationship between
   measurement noise (CV), the one-sided decision threshold (DT), and the number of patients
   correctly classified as responders or non-responders
4. **Tests robustness** under 8 distributional violations (log-normal IS, bimodal
   treatment effects, correlated noise, partial carryover, combined worst case)

The central finding: reducing measurement CV from 0.25 to 0.15 through
standardized sampling (fasting, timed draws, duplicate assays) is more efficient
than doubling the number of crossover cycles — saving 24 weeks of protocol time
at a cost of one extra blood tube per visit.

## Repository Structure

```
ksim/
├── simulations/          14 scripts — toxin models and intervention engineering
│   ├── gut_clearance_model.py          Baseline IS/PCS/urea generation & clearance
│   ├── drug_dev_simulation.py          Monte Carlo (N=20K), adherence, dose-response
│   ├── iter1_safety_feasibility.py     Side effects, feasibility gates, dialysis delay
│   ├── iter2_creative_bypass.py        Replace infeasible arms with deployable alternatives
│   ├── iter3_competing_risk.py         74-year-old competing risk model
│   ├── iter4_what_moves_needle.py      Toxin route vs AKI/CV prevention
│   ├── iter5_toxin_mortality_coupling.py  f_tox parameter sweep
│   ├── devA_sachet_formulation.py      Optimal fiber blend (RS/inulin/acacia)
│   ├── devC_collision_timing.py        Sorbent–probiotic timing separation
│   ├── dev_ketoacid_protein.py         Protein target + ketoacid optimization
│   ├── dev_integration_validation.py   Naive vs engineered stack comparison
│   ├── final_outcome_and_redteam.py    AST-120 drug interaction analysis
│   ├── redteam_loop.py                 15-round adversarial refinement (R: 48→82)
│   └── why_experts_failed_revalidation.py  10-pass expert failure escape analysis
│
├── protocol/             3 scripts — n-of-1 trial design and clinical decision engine
│   ├── nof1_virtual_cohort.py          N=100 virtual cohort, crossover power analysis
│   ├── nof1_weak_rescue.py             Weak responder strategies (adaptive + CV reduction)
│   └── clinical_protocol.py            Step-by-step protocol + classify_patient() function
│
├── validation/           8 scripts — robustness, calibration, and reproducibility
│   ├── robustness_experiments.py       8 distributional stress tests (N=1000 each)
│   ├── evsi_analysis.py                Classification-efficiency frontier (NCC vs CV, 50-seed mean)
│   ├── large_scale_validation.py       N=1000 end-to-end + published RCT reproduction
│   ├── reproduce_manuscript_numbers.py Single script regenerating all manuscript figures
│   ├── threshold_sensitivity.py        Responder-definition threshold sweep (5/10/15/20%)
│   ├── sensitivity_analysis.py         Tornado sweep across 8 literature-calibrated params
│   ├── literature_recalibration.py     Before/after comparison with published data
│   └── novelty_and_gap_analysis.py     RCT vs n-of-1 gap quantification
│
├── LIMITATIONS.md        Known limitations and anticipated reviewer responses
├── FUTURE_WORK.md        Research roadmap (14 directions, prioritized)
└── README.md             This file
```

## Quick Start

**Requirements**: Python 3.8+, NumPy. Matplotlib needed only for figure generation.

```bash
pip install numpy matplotlib

# Run any script independently — all are self-contained
python simulations/gut_clearance_model.py
python protocol/nof1_virtual_cohort.py
python validation/robustness_experiments.py

# Reproduce every number in the manuscript
python validation/reproduce_manuscript_numbers.py

# Run all simulations
for f in simulations/*.py protocol/*.py validation/*.py; do
  echo "=== $f ===" && python "$f"
done
```

## Key Results

| Metric | Value | Source script |
|--------|-------|--------------|
| Single-run sensitivity | 95.1% | `validation/reproduce_manuscript_numbers.py` |
| Single-run specificity | 85.7% | same |
| Worst-case sensitivity (combined violations) | 92.0% | `validation/robustness_experiments.py` |
| Worst-case specificity | 81.1% | same |
| Weak responder detection (τ 10–20%) | 76% | same |
| Peak NCC (classification-efficiency frontier) | CV ≈ 0.12–0.15 | `validation/evsi_analysis.py` |
| Average protocol duration | 26.7 weeks | `validation/reproduce_manuscript_numbers.py` |
| Patients needing Stage 2 | 23% | same |
| Specificity across responder thresholds (θ 5→20%) | 89.3→65.8% | `validation/threshold_sensitivity.py` |

## Parameter Sources

| Parameter | Value | Literature source |
|-----------|-------|-------------------|
| IS biological CV | 0.25 (conservative modeling baseline) | Pretorius et al. 2013, *Clin Chim Acta*, reports 35.9% for total serum IS in healthy volunteers. We simulate from 0.25, which understates rather than inflates the gain from standardization. |
| IS baseline CKD4 | 5.4 ± 3.6 µg/mL | Lin et al. 2011 |
| Fiber IS reduction (pooled) | SMD −0.34 | Wathanavasin et al. 2025, *Toxins*, 11 RCTs, N=398 |
| CKD4 eGFR slope | −2.0 mL/min/1.73 m²/yr | CRIC Study / MDRD |
| Non-responder fraction | 18% (conservative) | Estimated from meta I² + crossover data |
| Gut microbiome variability | SD 0.35 | Wu et al. 2011, *Science* |

## Clinical Protocol Summary

The protocol classifies individual CKD patients as IS responders or non-responders:

1. **Eligibility**: CKD 3b–4, stable eGFR, IS assay available
2. **Standardize**: Fasting AM draws, duplicate assays → CV ≈ 0.15
3. **Stage 1**: 2 × 3 AB crossover (24 weeks, 12 visits)
4. **Classify**: obs_red > 14% → Responder; obs_red < 0% → Non-responder; 0–14% → Borderline
5. **Stage 2**: Borderline patients get 1 extra cycle (wk 24–36), threshold 12%
6. **Follow-up**: Responders monitored quarterly; non-responders deprescribed

The decision engine is implemented in `protocol/clinical_protocol.py` as the
`classify_patient()` function — input IS measurements, output classification.

## License

Released under the [MIT License](LICENSE). Copyright (c) 2026 Pyeongwoo Lee and Timothy Juheon Lee.

This project is for research and educational purposes. Simulation parameters
are illustrative and calibrated to published data — not for clinical dosing
decisions without empirical validation.

## Citation

```
Lee P, Lee TJ. Measurement Noise Optimization as a Design Lever for N-of-1 Trials:
In-Silico Validation of Personalized Gut-Clearance Protocols in Chronic Kidney
Disease. 2026. https://github.com/dennis2lee/ksim
```

## Author

Pyeongwoo Lee — Independent Researcher, Sunnyvale, CA, USA
Timothy Juheon Lee — Independent Researcher, Sunnyvale, CA, USA
Contact: dennis2.lee@gmail.com
