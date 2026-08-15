# ARCHIVED. NOT PART OF THE PAPER.
#
# validation/reproduce_all.py does not run this file, and no result in
# the manuscript or its supplement depends on it. It is kept as a record
# of earlier work and is not maintained.
#
# This file is exploratory work on toxin models and intervention
# engineering. It informed no number in the paper.
#
# For what the paper actually does, read validation/nof1_core.py and the
# thirteen scripts the README lists under "What the paper reproduces".

"""
DEV-A: optimize the fiber-sachet recipe.
Decision: grams/day of resistant starch (RS), inulin/FOS (INU), acacia gum (ACA).
Maximize combined indoxyl-sulfate (IS) + p-cresyl-sulfate (PCS) reduction, SUBJECT TO
  - GI gas/tolerability budget (rapid fermenters bloat -> wreck adherence)
  - low potassium / phosphate (stage-4 ceilings)
  - total dose <= 25 g/day
Each fiber's per-toxin reduction saturates: r = Emax*(1-exp(-dose/k)).
"""
import numpy as np
import itertools

# fiber : Emax_IS, k_IS, Emax_PCS, k_PCS, gas_coeff, K_mg_per_g, P_mg_per_g
F = {
 "RS":  dict(eIS=0.45,kIS=12, ePCS=0.20,kPCS=12, gas=0.30, K=0.5, P=0.7),  # slow ferment, low gas, ~0 K/P; best for IS
 "INU": dict(eIS=0.25,kIS=8,  ePCS=0.45,kPCS=8,  gas=1.00, K=3.0, P=2.0),  # fast ferment, high gas; best for PCS
 "ACA": dict(eIS=0.15,kIS=10, ePCS=0.25,kPCS=10, gas=0.25, K=2.0, P=1.5),  # gentlest, padding/tolerability
}
GAS_BUDGET = 11.0     # tolerability cap (units); above this -> bloating -> drop-out
K_CAP, P_CAP = 90.0, 70.0   # mg/day from the sachet (keep well under daily budgets)
TOTAL_CAP = 25.0
TARGET_FERMENTABLE = 15.0    # want >=15 g fermentable for prebiotic effect

def red(dose, e, k):
    return e*(1-np.exp(-dose/k)) if dose>0 else 0.0

def evaluate(dRS,dINU,dACA):
    # combined reduction = 1 - product(1 - per-fiber reduction)
    IS = 1 - (1-red(dRS,F["RS"]["eIS"],F["RS"]["kIS"]))*(1-red(dINU,F["INU"]["eIS"],F["INU"]["kIS"]))*(1-red(dACA,F["ACA"]["eIS"],F["ACA"]["kIS"]))
    PCS= 1 - (1-red(dRS,F["RS"]["ePCS"],F["RS"]["kPCS"]))*(1-red(dINU,F["INU"]["ePCS"],F["INU"]["kPCS"]))*(1-red(dACA,F["ACA"]["ePCS"],F["ACA"]["kPCS"]))
    gas = dRS*F["RS"]["gas"]+dINU*F["INU"]["gas"]+dACA*F["ACA"]["gas"]
    K = dRS*F["RS"]["K"]+dINU*F["INU"]["K"]+dACA*F["ACA"]["K"]
    P = dRS*F["RS"]["P"]+dINU*F["INU"]["P"]+dACA*F["ACA"]["P"]
    ferment = dINU+dACA+0.4*dRS    # RS partially fermented
    return IS,PCS,gas,K,P,ferment

best=None
for dRS,dINU,dACA in itertools.product(range(0,21), range(0,21), range(0,21)):
    if dRS+dINU+dACA>TOTAL_CAP: continue
    IS,PCS,gas,K,P,ferment = evaluate(dRS,dINU,dACA)
    if gas>GAS_BUDGET or K>K_CAP or P>P_CAP or ferment<TARGET_FERMENTABLE: continue
    score = 0.5*IS+0.5*PCS
    if best is None or score>best[0]:
        best=(score,dRS,dINU,dACA,IS,PCS,gas,K,P,ferment)

print("="*70)
print("DEV-A  OPTIMAL FIBER-SACHET RECIPE (per day)")
print("="*70)
s,dRS,dINU,dACA,IS,PCS,gas,K,P,ferment = best
print(f"  Resistant starch (Hi-Maize/감자전분) : {dRS:>4} g")
print(f"  Inulin / FOS (치커리 이눌린)         : {dINU:>4} g")
print(f"  Acacia gum (아라비아검)              : {dACA:>4} g")
print(f"  --------------------------------------------")
print(f"  total                                : {dRS+dINU+dACA:>4} g   (fermentable {ferment:.0f} g)")
print(f"  expected IS reduction  : {IS*100:>4.0f} %")
print(f"  expected PCS reduction : {PCS*100:>4.0f} %")
print(f"  gas/tolerability       : {gas:>4.1f} / {GAS_BUDGET} budget")
print(f"  potassium load         : {K:>4.0f} mg/day")
print(f"  phosphate load         : {P:>4.0f} mg/day")

print("\n  Split: 2 sachets/day (breakfast + dinner) -> per sachet:")
print(f"    RS {dRS/2:.1f} g + INU {dINU/2:.1f} g + ACA {dACA/2:.1f} g, in a full glass of water")

print("\n  TITRATION (avoid bloating): start 1/3 dose wk1, 2/3 wk2, full wk3.")
for wk,frac in [("wk1",1/3),("wk2",2/3),("wk3+",1.0)]:
    print(f"    {wk:<5}: RS {dRS*frac:>4.1f} g  INU {dINU*frac:>4.1f} g  ACA {dACA*frac:>4.1f} g  (gas {gas*frac:.1f})")

# sensitivity: what does relaxing the gas budget buy?
print("\n  GAS-BUDGET SENSITIVITY (if she tolerates more fermentation):")
for gb in [8,11,14,18]:
    b=None
    for dRS,dINU,dACA in itertools.product(range(0,21),range(0,21),range(0,21)):
        if dRS+dINU+dACA>TOTAL_CAP: continue
        IS,PCS,gas,K,P,ferment=evaluate(dRS,dINU,dACA)
        if gas>gb or K>K_CAP or P>P_CAP or ferment<TARGET_FERMENTABLE: continue
        sc=0.5*IS+0.5*PCS
        if b is None or sc>b[0]: b=(sc,dRS,dINU,dACA,IS,PCS)
    if b: print(f"    gas<= {gb:>2}: RS{b[1]:>3} INU{b[2]:>3} ACA{b[3]:>3} -> IS {b[4]*100:.0f}% PCS {b[5]*100:.0f}%")
