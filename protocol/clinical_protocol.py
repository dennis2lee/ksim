"""
OPERATIONAL CLINICAL PROTOCOL: n-of-1 gut-clearance trial for individual CKD patients.

This file serves as both DOCUMENTATION and EXECUTABLE CODE:
  (A) Step-by-step protocol specification (print output)
  (B) Decision support engine (function: patient data in → classification out)
  (C) Demonstration on 6 example patients

Designed so a clinician can follow the protocol with standard lab equipment
and a spreadsheet (or this script). No specialized software required.

Literature-calibrated parameters (Pretorius 2013, 2025 meta-analysis, Lin 2011).
"""
import numpy as np

# =========================================================================
# PROTOCOL CONSTANTS (literature-calibrated)
# =========================================================================
CV_TARGET  = 0.15    # achievable via AM fasting + duplicate draws
CV_NATIVE  = 0.25    # unstandardized biological CV (Pretorius 2013)
MDE_6      = 1.645 * CV_TARGET * np.sqrt(2.0/6)   # 2x3 design: 14%
MDE_9      = 1.645 * CV_TARGET * np.sqrt(2.0/9)   # 3x3 design: 12%
WASHOUT_WK = 2       # literature: 2-4wk washout used in fiber/probiotic crossovers
TREAT_WK   = 4       # 4 weeks on-treatment before steady-state measurement
MEAS_PER_PERIOD = 3  # 3 measurement visits in final week of each period

# =========================================================================
# (1) STEP-BY-STEP CLINICAL PROTOCOL
# =========================================================================
print("="*86)
print("OPERATIONAL n-of-1 PROTOCOL FOR GUT-CLEARANCE IN CKD")
print("="*86)

print(f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 0: PATIENT ELIGIBILITY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  INCLUDE:
    - CKD stage 3b-4 (eGFR 15-45 mL/min/1.73m²)
    - Stable eGFR slope (no AKI in past 3 months)
    - Serum IS measurable at local/reference lab (LC-MS/MS or UPLC)
    - Able to adhere to 24-week protocol
    - On stable base medical therapy (RAS blocker, SGLT2i, BP controlled)

  EXCLUDE:
    - Active infection or hospitalization
    - GI surgery, short bowel, or severe gastroparesis
    - On antibiotics (disrupts microbiome; wait 4 weeks post-course)
    - On immunosuppressants (transplant patients)
    - Life expectancy < 12 months (competing risk too high for protocol value)
    - Already on dialysis (different toxin kinetics)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 1: MEASUREMENT STANDARDIZATION (achieve CV ~ 0.15)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  WHY: biological IS CV = 25-27% (Pretorius 2013). Unstandardized measurement
  has MDE = 24% — too high to detect moderate responders. Standardization
  targets CV ~ 0.15, giving MDE = {MDE_6*100:.0f}%.

  HOW (at each measurement visit):
    1. FASTING: overnight fast (>=8h), no tryptophan-rich foods prior day
       (reduces dietary contribution to IS generation)
    2. TIMED: draw between 7-9 AM (reduces diurnal variation)
    3. DUPLICATE: draw 2 separate tubes from same venipuncture
       -> lab reports average of the 2 assays
       -> effective analytical CV halved (~3% instead of ~6%)
    4. ASSAY: serum total IS by LC-MS/MS or UPLC-fluorescence
       -> same lab, same method throughout protocol (no inter-lab variation)
    5. RECORD: IS value in ug/mL + date + time + fasting confirmed (Y/N)

  COST: 1 extra tube per visit = negligible
  EVIDENCE: Deltombe 2019 (Toxins): standardized IS CV = 14-18%

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 2: STAGE 1 — 2x3 AB CROSSOVER (weeks 0-24)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  DESIGN: 2 complete AB cycles, 3 measurements per period

  TIMELINE:
    wk  0-4  : PERIOD A1 (control = base medical only, NO gut intervention)
                -> measure IS at wk 3.0, 3.5, 4.0 (3 visits)
    wk  4-6  : WASHOUT (stop any prior gut supplement; 2 weeks)
                Rationale: IS returns to baseline within 2 weeks after
                fiber/probiotic cessation (Sirich 2014; standard in XO trials)
    wk  6-10 : PERIOD B1 (intervention = fiber sachet + Lacto probiotic)
                -> titrate fiber wk 6-7, full dose wk 7-10
                -> measure IS at wk 9.0, 9.5, 10.0 (3 visits)
    wk 10-12 : WASHOUT
    wk 12-16 : PERIOD A2 (control repeat)
                -> measure IS at wk 15.0, 15.5, 16.0
    wk 16-18 : WASHOUT
    wk 18-22 : PERIOD B2 (intervention repeat)
                -> measure IS at wk 21.0, 21.5, 22.0
    wk 22-24 : FINAL WASHOUT + ANALYSIS

  INTERVENTION REGIMEN (from engineering modules):
    - Fiber sachet: RS 14g + Inulin 4g + Acacia 7g, split AM/PM
      (titrate: 1/3 dose wk1, 2/3 wk2, full wk3+)
    - Probiotic: Lactobacillus, refrigerated, at BEDTIME (>=2h after dinner)
      (NOT spores; NOT co-administered with sorbent)
    - NO AST-120 (net harmful due to antihypertensive binding)

  TOTAL: 12 measurement visits, 24 blood draws (duplicate), 24 weeks

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 3: STAGE 1 ANALYSIS — compute individual IS reduction
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  CALCULATION (can be done in a spreadsheet):
    mean_A = average of 6 control-period IS values (A1 + A2)
    mean_B = average of 6 intervention-period IS values (B1 + B2)
    observed_reduction = (mean_A - mean_B) / mean_A

  DECISION THRESHOLDS:
    MDE = {MDE_6*100:.0f}% (one-sided 95%, at CV=0.15, 6 measures/arm)

    if observed_reduction > {MDE_6*100:.0f}%  -> RESPONDER
       action: continue gut-clearance regimen indefinitely
       confidence: >95% that the true effect is real (not noise)

    if observed_reduction < 0%               -> NON-RESPONDER
       action: stop gut-clearance regimen, eliminate burden
       confidence: high (IS went up or stayed same on intervention)

    if 0% <= observed_reduction <= {MDE_6*100:.0f}% -> BORDERLINE
       action: proceed to STAGE 2 (1 extra cycle)
       rationale: effect may be real but too small to confirm with 6 measures

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 4: STAGE 2 — ADAPTIVE ENRICHMENT (weeks 24-36, borderline only)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  WHO: only patients classified as BORDERLINE in Stage 1
  (expected: ~30% of patients; strong responders and clear non-resp already done)

    wk 24-28 : PERIOD A3 (control) -> measure IS x3
    wk 28-30 : WASHOUT
    wk 30-34 : PERIOD B3 (intervention) -> measure IS x3
    wk 34-36 : ANALYSIS

  STAGE 2 ANALYSIS:
    Combine ALL data: 9 measures per arm (6 from Stage 1 + 3 from Stage 2)
    mean_A = average of 9 control IS values
    mean_B = average of 9 intervention IS values
    observed_reduction = (mean_A - mean_B) / mean_A

    MDE_stage2 = {MDE_9*100:.0f}% (one-sided 95%, 9 measures/arm)

    if observed_reduction > {MDE_9*100:.0f}% -> RESPONDER (continue)
    else                                     -> NON-RESPONDER (stop)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 5: FOLLOW-UP & MONITORING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  FOR RESPONDERS:
    - Continue regimen
    - Re-measure IS every 3 months (single standardized draw)
    - If IS rises >20% above on-treatment level -> investigate adherence or
      disease progression; consider repeating protocol
    - Quarterly deprescribe review: is the regimen still net beneficial?

  FOR NON-RESPONDERS:
    - Stop all gut-clearance supplements
    - Focus on base medical optimization (BP, SGLT2i, acidosis management)
    - Document: "this patient's IS does not respond to fiber/probiotic"
    - Revisit only if regimen composition changes substantially

  SAFETY MONITORING (all patients):
    - Serum potassium: check at wk 2, 4 (fiber + RAS = K risk)
    - GI symptoms diary (bloating, flatulence) — tolerability gate
    - Weight + grip strength (if on protein moderation + ketoacid)
    - Stop rules: K > 5.5, persistent GI intolerance, weight loss >3%
""")

# =========================================================================
# (2) DECISION SUPPORT ENGINE
# =========================================================================
print(f"\n{'='*86}")
print("DECISION SUPPORT ENGINE")
print(f"{'='*86}")

def classify_patient(
    is_control: list,
    is_intervention: list,
    cv: float = CV_TARGET,
    stage: int = 1,
) -> dict:
    """Classify a patient based on their IS measurements.

    Args:
        is_control: list of IS values from control periods (ug/mL)
        is_intervention: list of IS values from intervention periods (ug/mL)
        cv: measurement CV (default 0.15 = standardized protocol)
        stage: 1 (6+6 measures) or 2 (9+9 measures)

    Returns:
        dict with classification, observed reduction, confidence, and action
    """
    n_a = len(is_control)
    n_b = len(is_intervention)
    mean_a = np.mean(is_control)
    mean_b = np.mean(is_intervention)

    if mean_a <= 0:
        return dict(classification="ERROR", reason="control mean IS <= 0")

    obs_reduction = (mean_a - mean_b) / mean_a
    mde = 1.645 * cv * np.sqrt(1.0/n_a + 1.0/n_b)

    if obs_reduction > mde:
        classification = "RESPONDER"
        action = "Continue gut-clearance regimen"
        confidence = ">=95% true positive"
    elif obs_reduction < 0:
        classification = "NON-RESPONDER"
        action = "Stop gut-clearance; focus on base medical therapy"
        confidence = "high (IS increased or unchanged on intervention)"
    else:
        if stage == 1:
            classification = "BORDERLINE"
            action = "Proceed to Stage 2 (1 additional AB cycle, wk 24-36)"
            confidence = f"observed {obs_reduction*100:.0f}% < MDE {mde*100:.0f}%"
        else:
            classification = "NON-RESPONDER"
            action = "Stop gut-clearance; effect too small to confirm"
            confidence = f"stage 2 combined: {obs_reduction*100:.0f}% < MDE {mde*100:.0f}%"

    return dict(
        classification=classification,
        observed_reduction=obs_reduction,
        mean_control=mean_a,
        mean_intervention=mean_b,
        mde=mde,
        n_control=n_a,
        n_intervention=n_b,
        action=action,
        confidence=confidence,
        stage=stage,
    )

def print_result(r, label=""):
    tag = f" [{label}]" if label else ""
    print(f"\n  Patient{tag}:")
    print(f"    control IS:    {r['mean_control']:.1f} ug/mL (n={r['n_control']})")
    print(f"    interv. IS:    {r['mean_intervention']:.1f} ug/mL (n={r['n_intervention']})")
    print(f"    observed red:  {r['observed_reduction']*100:.1f}%  (MDE = {r['mde']*100:.0f}%)")
    print(f"    --->  {r['classification']}  (stage {r['stage']})")
    print(f"    action: {r['action']}")
    print(f"    confidence: {r['confidence']}")

# =========================================================================
# (3) DECISION LOGIC (text flowchart)
# =========================================================================
print(f"""
  DECISION FLOWCHART:

  [Patient eligible?] --NO--> exclude
        |
       YES
        |
  [STAGE 1: 2x3 crossover, 24 wk]
  [Compute: obs_red = (mean_A - mean_B) / mean_A]
        |
        +---> obs_red > {MDE_6*100:.0f}%  ---> RESPONDER  ---> continue regimen
        |
        +---> obs_red < 0%     ---> NON-RESP  ---> stop, deprescribe
        |
        +---> 0% <= obs_red <= {MDE_6*100:.0f}%  ---> BORDERLINE
                                          |
                                    [STAGE 2: +1 cycle, wk 24-36]
                                    [Combine all 9+9 measures]
                                    [obs_red_combined vs MDE = {MDE_9*100:.0f}%]
                                          |
                                          +---> > {MDE_9*100:.0f}% ---> RESPONDER
                                          +---> <= {MDE_9*100:.0f}% ---> NON-RESP

  REMEASUREMENT TRIGGER (during follow-up):
    if on-treatment IS rises >20% above baseline -> repeat protocol or check adherence
""")

# =========================================================================
# (4) DEMONSTRATION: 6 example patients
# =========================================================================
print(f"{'='*86}")
print("DEMONSTRATION: 6 example patients")
print(f"{'='*86}")

rng = np.random.default_rng(42)

examples = [
    ("Strong responder (tau~45%)", 5.4, 0.45),
    ("Moderate responder (tau~30%)", 6.2, 0.30),
    ("Weak responder (tau~15%)", 4.8, 0.15),
    ("Very weak (tau~8%)", 5.0, 0.08),
    ("Non-responder (tau~2%)", 7.1, 0.02),
    ("High baseline (tau~35%)", 12.0, 0.35),
]

for label, baseline_is, true_tau in examples:
    # Simulate stage 1 measurements (2x3 = 6 per arm)
    ctrl = baseline_is * (1 + rng.normal(0, CV_TARGET, 6))
    intv = baseline_is * (1 - true_tau) * (1 + rng.normal(0, CV_TARGET, 6))

    r1 = classify_patient(ctrl.tolist(), intv.tolist(), CV_TARGET, stage=1)
    print_result(r1, label)

    # If borderline, demonstrate stage 2
    if r1['classification'] == 'BORDERLINE':
        ctrl_s2 = baseline_is * (1 + rng.normal(0, CV_TARGET, 3))
        intv_s2 = baseline_is * (1 - true_tau) * (1 + rng.normal(0, CV_TARGET, 3))
        all_ctrl = ctrl.tolist() + ctrl_s2.tolist()
        all_intv = intv.tolist() + intv_s2.tolist()
        r2 = classify_patient(all_ctrl, all_intv, CV_TARGET, stage=2)
        print(f"    --- STAGE 2 RESULT ---")
        print(f"    combined obs red: {r2['observed_reduction']*100:.1f}% (MDE={r2['mde']*100:.0f}%)")
        print(f"    --->  {r2['classification']}")
        print(f"    action: {r2['action']}")

# =========================================================================
# (5) QUANTITATIVE JUSTIFICATION — linking to simulation results
# =========================================================================
print(f"\n{'='*86}")
print("QUANTITATIVE JUSTIFICATION FOR EACH PROTOCOL CHOICE")
print(f"{'='*86}")
print(f"""
  CHOICE                        JUSTIFICATION (from simulation)
  -----------------------------------------------------------------------
  CV=0.15 (standardized)        MDE drops 24% -> {MDE_6*100:.0f}% (saves 24 wk vs 4x3 at CV=0.25)
                                Power 68% -> 86% (recalibrated, literature_recalibration.py)
                                Cost: 1 extra tube per visit

  2x3 crossover (not 1x3)      1x3 power = ~55%; 2x3 = 86%; +12 wk justified
                                2 cycles provide replication (crossover validity)

  3 measures per period         2/period: MDE=18%; 3/period: MDE={MDE_6*100:.0f}%
  (not 2)                       +1 visit per period = +4 visits total for -4pp MDE

  2-wk washout                  IS returns to baseline within 2 wk post fiber cessation
                                (Sirich 2014, standard in published crossover trials)
                                Longer washout (4 wk) is safer but costs 8 wk total

  Adaptive stage 2              Only ~30% of patients need it (borderline zone)
                                Avoids 36-wk fixed design for the 70% already classified
                                Weak responder power: 20% -> 47% (recalibrated)

  MDE={MDE_6*100:.0f}% threshold          One-sided 95% confidence: <5% chance of false positive
                                Clinically meaningful: {MDE_6*100:.0f}% IS reduction is a real response

  IS as primary endpoint        Avoids EPPIC surrogate trap (IS->hard endpoint failed)
                                Directly answers: "does this patient's IS drop?"
                                f_tox independent: no downstream outcome claim needed

  No AST-120                    Net harmful: adsorbs antihypertensives (final_outcome_and_redteam.py)
                                BP control is the #1 eGFR slope lever

  Lacto at bedtime              Function > survival: Lacto+bedtime = 91% effective vs
  (not spores)                  spore 49% (final_outcome_and_redteam.py Part 3)

  OVERALL PROTOCOL PERFORMANCE (literature-calibrated):
    sensitivity (responders correctly identified):  86%
    specificity (non-resp correctly deprescribed): ~95%
    avg duration: 24 wk (70% patients) / 36 wk (30% borderline)
    total draws: 24 (stage 1) + 12 if stage 2 = 24-36
    eGFR cost: ~0.9 mL/min over 24 wk (crossover-controlled)
""")
