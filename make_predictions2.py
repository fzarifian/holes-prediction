"""Freeze the second set of blind predictions (radii 2x10^5 to 10^6) -> predictions_v2/predictions_v2.json.

Model: predict2.Model2 (Poisson model + Hardy-Littlewood pair correlations), with no fitted parameter.
It was checked, before this freeze, on the 54 configurations already measured (radii 500 to 10^5):
record radius exact in 41 of 54 and within one lattice shell in 53 of 54; direction of the bias at the
measured record radius correct in 104 of 114 clear predictions.

Test annuli: R <= |z| <= R + 100 for R in TEST_RADII, inside the angular window 0 <= arg z <= phi with
phi chosen to give about 4 million lattice points (phi <= the fundamental sector).
"""
import json
import math
import os

import numpy as np

from holes import FIELDS
from predict2 import Model2

TEST_RADII = (200_000, 450_000, 1_000_000)
WIDTH = 100
TARGET_POINTS = 4_000_000
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'predictions_v2')
os.makedirs(OUT, exist_ok=True)


def rmax_for(D):
    return {-3: 7.2, -4: 7.6}.get(D, 10.6)


class Interp:
    """P(h >= shell | m) on a grid of log N, interpolated linearly in 1/log N."""

    def __init__(self, M, Rlo, Rhi, n=12):
        self.M = M
        self.L = np.log(np.exp(np.linspace(math.log(Rlo), math.log(Rhi), n)) ** 2)
        self.P = np.stack([np.exp(M.log_tail(L)[0]) for L in self.L])          # (n, #m, #shells)
        self.x = 1 / self.L

    def at(self, R):
        x = 1 / math.log(R ** 2)
        i = int(np.clip(np.searchsorted(-self.x, -x) - 1, 0, len(self.x) - 2))
        t = (x - self.x[i]) / (self.x[i + 1] - self.x[i])
        return (1 - t) * self.P[i] + t * self.P[i + 1]

    def predict(self, R, q=0.999, rstar=None):
        M = self.M
        P = self.at(R + WIDTH / 2)
        logN = math.log((R + WIDTH / 2) ** 2)
        dens = M.dens.copy()
        dens[0] *= 1 - M.ring.kappa / logN * M.prime_factor
        dens /= dens.sum()
        tail = dens @ P
        k = int(np.max(np.where(tail >= 1 - q)[0])) if rstar is None else int(np.argmin(np.abs(M.shells - rstar)))
        w = P[:, k]; tot = dens @ w
        by_omega, share = {}, {}
        for j in range(int(M.omega.max()) + 1):
            s = M.omega == j
            if dens[s].sum() > 0:
                by_omega[str(j)] = round(float((dens[s] @ w[s]) / dens[s].sum() / tot), 4)
                share[str(j)] = float(dens[s].sum())
        by_ideal = [{'norm': n, 'gen': [int(pa), int(pb)],
                     'ratio': round(float((dens[M.bits[:, j]] @ w[M.bits[:, j]]) / dens[M.bits[:, j]].sum() / tot), 4)}
                    for j, (pa, pb, n) in enumerate(M.sp)]
        return k, {'by_omega': by_omega, 'class_share': share, 'by_ideal': by_ideal}


def main():
    pred = {'model': 'Poisson + Hardy-Littlewood pair correlations (predict2.py), y = 31, quantile 0.999, '
                     'no fitted parameter', 'width': WIDTH, 'fields': {}}
    for D in FIELDS:
        M = Model2(D, rmax_for(D))
        I = Interp(M, 1.5e5, 1.05e6)
        Rs = np.exp(np.linspace(math.log(1e5), math.log(1e6), 600))
        ks = [I.predict(R)[0] for R in Rs]
        trans = []
        for i in range(1, len(Rs)):
            if ks[i] != ks[i - 1]:
                lo, hi = Rs[i - 1], Rs[i]
                for _ in range(30):
                    mid = math.sqrt(lo * hi)
                    if I.predict(mid)[0] == ks[i - 1]:
                        lo = mid
                    else:
                        hi = mid
                trans.append({'R': round(float(hi)), 'rstar_from': float(M.shells[ks[i - 1]]),
                              'rstar_to': float(M.shells[ks[i]])})
        tests = []
        covol = math.sqrt(-D) / 2
        for R in TEST_RADII:
            phi = min(M.ring.theta, TARGET_POINTS * covol / (R * WIDTH))
            k, _ = I.predict(R)
            cond = {}
            for kk in range(max(0, k - 2), min(len(M.shells), k + 3)):
                cond[f'{M.shells[kk]:.6f}'] = I.predict(R, rstar=float(M.shells[kk]))[1]
            tests.append({'R1': R, 'R2': R + WIDTH, 'phi': phi, 'rstar_pred': float(M.shells[k]), 'conditional': cond})
        pred['fields'][str(D)] = {'rstar_at_1e5': float(M.shells[ks[0]]), 'rstar_transitions': trans, 'tests': tests}
        print(f'D={D:5d}: r*(1e5)={M.shells[ks[0]]:.3f}, transitions {[(t["R"], round(t["rstar_to"], 3)) for t in trans]}')
        for c in tests:
            pat = c['conditional'][f'{c["rstar_pred"]:.6f}']
            show = {k: v for k, v in pat['by_omega'].items() if pat['class_share'][k] >= 0.025}
            print(f'    [{c["R1"]}, {c["R2"]}] phi={c["phi"]:.4f} r*={c["rstar_pred"]:.3f} {show}', flush=True)
    with open(os.path.join(OUT, 'predictions_v2.json'), 'w') as f:
        json.dump(pred, f, indent=1)


if __name__ == '__main__':
    main()
