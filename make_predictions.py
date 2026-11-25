"""Freeze the blind predictions for radii 10^4 to 10^5 (writes predictions/predictions.json).

Model: predict.Model (singular-series Poisson model, densities only, no prime computation).
  * Bias patterns use the plain model (kappa as in the prime ideal theorem).
  * The record radius r* is predicted with kappa_eff = C_RSTAR * kappa, where C_RSTAR = 1.25 is the single
    parameter fitted on the 27 configurations of the paper (radii 500 to 4600): it gives the measured r*
    exactly in 20 of 27 cases.  Since r* may still be off by one lattice shell, the bias pattern is also
    given conditionally on each lattice distance near the predicted r*; the verification uses the pattern
    at the measured r*.
"""
import json
import math
import os

import numpy as np

from holes import FIELDS, QF
from predict import Model

C_RSTAR = 1.25
TARGET_POINTS = 4_000_000          # approximate number of lattice points per test annulus
TEST_RADII = (10_000, 31_623, 100_000)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'predictions')
os.makedirs(OUT, exist_ok=True)


def width(ring, R):
    covol = math.sqrt(-ring.D) / 2
    return int(round(TARGET_POINTS * covol / (ring.theta * R)))


def rstar_pred(M, R1, R2):
    k0 = M.ring.kappa
    M.ring.kappa = k0 * C_RSTAR
    try:
        return M.annulus(R1, R2)['rstar']
    finally:
        M.ring.kappa = k0


def transitions(M, Rmin=10_000, Rmax=100_000, steps=300):
    Rs = np.exp(np.linspace(math.log(Rmin), math.log(Rmax), steps))
    vals = [rstar_pred(M, R, R + width(M.ring, R)) for R in Rs]
    out = []
    for i in range(1, len(Rs)):
        if vals[i] != vals[i - 1]:
            lo, hi = Rs[i - 1], Rs[i]
            for _ in range(25):
                mid = math.sqrt(lo * hi)
                if rstar_pred(M, mid, mid + width(M.ring, mid)) == vals[i - 1]:
                    lo = mid
                else:
                    hi = mid
            out.append({'R': round(float(hi)), 'rstar_from': vals[i - 1], 'rstar_to': vals[i]})
    return vals[0], out


def main():
    pred = {'model': 'singular-series Poisson model, y = 31, quantile 0.999', 'C_RSTAR': C_RSTAR,
            'fields': {}}
    for D in FIELDS:
        M = Model(D)
        r0, tr = transitions(M)
        cfgs = []
        for R in TEST_RADII:
            W = width(M.ring, R)
            rs = rstar_pred(M, R, R + W)
            k = int(np.argmin(np.abs(M.shells - rs)))
            cond = {}
            for kk in range(max(0, k - 2), min(len(M.shells), k + 3)):
                p = M.annulus(R, R + W, rstar=float(M.shells[kk]))
                cond[f'{M.shells[kk]:.6f}'] = {'by_omega': {str(j): round(v, 4) for j, v in p['by_omega'].items()},
                                               'class_share': {str(j): v for j, v in p['class_share'].items()},
                                               'by_ideal': [{'norm': d['norm'], 'gen': d['gen'],
                                                             'ratio': round(d['ratio'], 4)} for d in p['by_ideal']]}
            cfgs.append({'R1': R, 'R2': R + W, 'rstar_pred': rs, 'conditional': cond})
        pred['fields'][str(D)] = {'rstar_at_1e4': r0, 'rstar_transitions': tr, 'tests': cfgs}
        print(f'D={D:5d}: r*(1e4)={r0:.3f}, transitions {[(t["R"], round(t["rstar_to"], 3)) for t in tr]}')
        for c in cfgs:
            pat = c['conditional'][f'{c["rstar_pred"]:.6f}']['by_omega']
            print(f'    [{c["R1"]}, {c["R2"]}] r*={c["rstar_pred"]:.3f}  omega pattern {pat}')
    with open(os.path.join(OUT, 'predictions.json'), 'w') as f:
        json.dump(pred, f, indent=1)


if __name__ == '__main__':
    main()
