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
DEV-C: collision simulation. AST-120 (activated carbon) adsorbs/inactivates the
probiotic when both occupy the same gut zone at the same time. Question: how many
hours must separate the probiotic dose from the sorbent dose to keep most bacteria alive?

Model: each oral dose creates a 'presence pulse' in the colonic interaction zone,
described by a gamma residence curve (mean transit ~6 h, spread). Probiotic taken at
t=0, sorbent at t=Δ. Killing is proportional to the temporal OVERLAP of the two pulses.
Surviving fraction = 1 - maxkill * overlap(Δ)/overlap(0).
"""
import numpy as np
from math import gamma as gammafn

def gpdf(x, a, scale):
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    m = x > 0
    out[m] = x[m]**(a-1)*np.exp(-x[m]/scale)/(scale**a*gammafn(a))
    return out

# gamma residence pulse: mean = shape*scale, choose mean~6h, sd~3h
shape, scale = 4.0, 1.5     # mean 6h, sd = sqrt(4)*1.5 = 3h
t = np.linspace(0, 24, 2400)
dt = t[1]-t[0]

def pulse(shift):
    return gpdf(t - shift, shape, scale)

P0 = pulse(0.0)                 # probiotic taken at t=0
MAXKILL = 0.70                  # fraction killed if taken together (Δ=0)
overlap0 = np.sum(P0*pulse(0.0))*dt

def survival(delta):
    S = pulse(delta)            # sorbent taken at t=delta
    ov = np.sum(P0*S)*dt
    return 1 - MAXKILL*(ov/overlap0)

print("="*64)
print("DEV-C  PROBIOTIC SURVIVAL vs SORBENT SEPARATION INTERVAL")
print("="*64)
print(f"{'separation Δ (h)':>18}{'probiotic survival':>22}")
target=0.80; rec=None
for d in [0,0.5,1,1.5,2,2.5,3,3.5,4,5,6,8,10]:
    s=survival(d)
    if s>=target and rec is None: rec=d
    flag = "  <-- first >=80%" if d==rec else ""
    bar="#"*int(s*30)
    print(f"{d:>18.1f}{s*100:>19.0f}%  {bar}{flag}")

rectxt = f"~{rec:.1f} h" if rec is not None else ">10 h (impractical)"
print(f"\nMinimum separation for >=80% probiotic survival (broad/worst-case): {rectxt}")

# sensitivity to the INTERACTION-WINDOW WIDTH (the real driver) and kill strength
print("\nSENSITIVITY: interaction-window width (sd) -- narrower zone = timing works better")
for sd in [1.0, 2.0, 3.0]:
    sh=(6.0/sd)**2; sc=6.0/sh
    Pm=gpdf(t,sh,sc); ov0=np.sum(Pm*Pm)*dt
    rec2=None
    for d in np.arange(0,10.01,0.25):
        s=1-MAXKILL*(np.sum(Pm*gpdf(t-d,sh,sc))*dt/ov0)
        if s>=0.80 and rec2 is None: rec2=d
    print(f"  sd={sd:.0f} h  -> need {('~%.1f h'%rec2) if rec2 is not None else '>10 h'} separation")

print("\nSENSITIVITY: max-kill strength (more selective sorbent = less killing)")
for mk in [0.7,0.5,0.3]:
    rec3=None
    for d in np.arange(0,10.01,0.25):
        s=1-mk*(np.sum(P0*pulse(d))*dt/overlap0)
        if s>=0.80 and rec3 is None: rec3=d
    print(f"  maxkill={mk:.1f} -> need {('~%.1f h'%rec3) if rec3 is not None else '>10 h'} separation")
