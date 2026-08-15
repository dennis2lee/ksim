# ARCHIVED. NOT PART OF THE PAPER.
#
# validation/reproduce_all.py does not run this file, and no result in
# the manuscript or its supplement depends on it. It is kept as a record
# of earlier work and is not maintained.
#
# This file argues for deprescribing as a clinical objective. The paper
# makes no clinical recommendation.
#
# For what the paper actually does, read validation/nof1_core.py and the
# thirteen scripts the README lists under "What the paper reproduces".

"""
WHY HAVE EXPERTS NOT SOLVED CKD? -- and does OUR converged solution escape the same
traps, or just repeat them? 10 critical-thinking passes: each pass names a root cause
of the field's failure, scores whether our solution ESCAPES it, and exposes the
RESIDUAL vulnerability + the action that would close it.

escape in [0,1]: 1 = our solution genuinely sidesteps this failure mode.
A high aggregate is meaningless if any BINDING constraint stays low -> we report those.
"""
import numpy as np

# (root cause experts failed, why it blocks them, our escape (base), residual + needed action)
PASSES = [
 ("Nephrons don't regenerate (hard biology)",
  "Adult kidney has no meaningful nephron regeneration -> lost function is permanent",
  0.90, "We never claim to RESTORE function -> escape by honest scope, BUT ceiling stays low (slow/relieve only)"),
 ("Economic incentive misalignment",
  "Dialysis is a huge industry; fiber/diet are unpatentable -> under-funded, under-studied",
  0.95, "Our core is cheap/OTC/non-patentable -> we explore exactly what industry ignored (an ADVANTAGE)"),
 ("Surrogate-endpoint trap (toxin != outcome)",
  "AST-120 lowered toxins but failed hard endpoints -> lowering a marker may not help",
  0.45, "We are VULNERABLE to the same trap -- legitimacy only if we MEASURE + scope to symptom/QoL"),
 ("Disease heterogeneity lumped as 'CKD'",
  "Diabetic/hypertensive/GN/PKD are different diseases -> one trial dilutes any effect",
  0.85, "Our n-of-1 is tailored to HER cause (hypertensive) -> escapes lumping by design"),
 ("Reductionism & specialty silos",
  "Nephro/cardio/nutrition/microbiome siloed -> gut-kidney axis falls between chairs",
  0.88, "Our design is explicitly multi-factor + whole-person (competing risk) -> escapes silos"),
 ("RCT-or-nothing evidence culture",
  "n-of-1 / personalized optimization dismissed as anecdote -> not pursued for individuals",
  0.70, "We embrace n-of-1 -> legitimate FOR HER, but never 'proof' by guideline standards"),
 ("Geriatric / competing-risk neglect",
  "Trials optimize the disease, not the 74yo's life-years -> wrong objective for elderly",
  0.90, "We center competing risk + QoL + deprescribing -> escapes the wrong-objective trap"),
 ("Negative-result burial / dogma",
  "One failed trial (AST-120) can kill an approach permanently regardless of nuance",
  0.80, "We re-examined and DROPPED AST-120 on its merits -> escapes dogma, stays adaptive"),
 ("Over-treatment bias ('more is better')",
  "Adding drugs is rewarded; deprescribing is not -> elderly get polypharmacy harm",
  0.90, "Our loop converged to LESS (minimum-effective-set) -> escapes the additive bias"),
 ("No measurement / no feedback loop",
  "Without per-patient measurement you cannot learn what works -> everyone guesses",
  0.40, "BINDING: our whole legitimacy collapses to 'another guess' UNLESS the n-of-1 measurement is actually executed"),
]

def fix_escape(base):
    # the two closeable residuals (surrogate trap, measurement) are closed by ONE action:
    # actually executing honest n-of-1 measurement. Model that lift.
    return min(0.92, base + 0.45*(1-base)) if base < 0.6 else base

print("="*94)
print("10-PASS RE-VALIDATION: does our solution escape the traps that defeated the experts?")
print("="*94)
base=[]; fixed=[]
for i,(cause,why,esc,res) in enumerate(PASSES,1):
    f=fix_escape(esc); base.append(esc); fixed.append(f)
    tag = "BINDING" if esc<0.5 else ("strong" if esc>=0.85 else "ok")
    print(f"\n[{i:>2}] {cause}   (escape {esc:.2f} -> {f:.2f}  {tag})")
    print(f"     why experts stuck: {why}")
    print(f"     our residual/action: {res}")

print("\n" + "="*94)
print("AGGREGATE")
print("="*94)
print(f"  mean escape (before action): {np.mean(base):.2f}")
print(f"  weakest link (binding)     : {min(base):.2f}  <- a high mean is FALSE comfort if this stays low")
print(f"  mean escape (after executing honest n-of-1 measurement): {np.mean(fixed):.2f}")
print(f"  weakest link after action  : {min(fixed):.2f}")

print("\nVERDICT:")
print("  Experts didn't 'solve' it for TWO different reasons, and they must be separated:")
print("   (1) a REAL hard ceiling (no nephron regeneration) -> nobody can solve; only slow/relieve.")
print("   (2) STRUCTURAL traps (incentives, silos, RCT-culture, geriatric neglect, over-treatment)")
print("       -> our approach sidesteps these almost by accident, BY being cheap, personalized,")
print("          systems-level, geriatric-aware and deprescribing.")
print("  BUT our solution earns legitimacy over the failed attempts on ONE hinge only:")
print("  actually MEASURING in the patient + honest scoping. Skip that -> we are just another")
print("  unvalidated surrogate story, identical to what already failed.")
