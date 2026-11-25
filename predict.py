"""Blind predictions from the singular-series Poisson model, using only densities (no prime computation).

For a field D and an annulus R1 <= |z| <= R2, the model gives, for each signature m (set of small prime
ideals dividing z, density from the Chinese remainder theorem),
    P(h(z) >= d | m) = E_N[ exp(-kappa * Lambda_m(d-) / log N) ],     Lambda_m(d-) = sum_{0<|u|<d} S_m(u),
where E_N averages over the norms N = |z|^2 of the annulus (area weighting).  From this:
  * r*  = largest lattice distance d with P(h >= d) >= 1 - q   (q = 0.999, as np.quantile on the data),
  * the predicted over-representation of each omega class, and of 'p | z' for each small ideal p,
    among the points with h >= r*.
Non-prime points only: the class m = 0 loses the primes (probability kappa/log N * prod N/(N-1)).
"""
import json
import math
import os

import numpy as np

from holes import FIELDS, QF, offsets

Y = 31


class Model:
    def __init__(self, D, y=Y, rmax=13.0):
        self.ring = ring = QF(D)
        self.sp = sp = ring.small_primes(y)
        s = len(sp)
        self.norms = np.array([n for _, _, n in sp], float)
        ua, ub, d = offsets(ring, rmax)
        o = np.argsort(d, kind='stable'); ua, ub, d = ua[o], ub[o], d[o]
        self.shells = np.unique(np.round(d, 9))
        # index of the first offset with |u| >= shell, so that offsets [:idx] are |u| < shell
        self.first = np.searchsorted(np.round(d, 9), self.shells, side='left')
        udiv = np.stack([ring.divides(pa, pb, n, ua, ub) for pa, pb, n in sp], 1)       # (#u, s)
        masks = np.arange(1 << s)
        bits = ((masks[:, None] >> np.arange(s)[None, :]) & 1).astype(bool)                 # (#m, s)
        self.bits = bits
        self.omega = bits.sum(1)
        N = self.norms
        self.dens = np.prod(np.where(bits, 1 / N, 1 - 1 / N), axis=1)
        one, none = N / (N - 1), N * (N - 2) / (N - 1) ** 2
        # Lambda_m(shell-) for every mask and shell, built by blocks of masks
        cum = np.zeros((len(masks), len(self.shells)))
        for a in range(0, len(masks), 512):
            zb = bits[a:a + 512][:, None, :]                       # (b, 1, s)
            ud = udiv[None, :, :]                                   # (1, #u, s)
            f = np.where(zb & ud, 0.0, np.where(zb ^ ud, one, none))
            S = f.prod(axis=2)                                      # (b, #u)
            c = np.concatenate([np.zeros((S.shape[0], 1)), np.cumsum(S, axis=1)], axis=1)
            cum[a:a + 512] = c[:, self.first]
        self.lam = cum
        self.prime_factor = float(np.prod(N / (N - 1)))

    def annulus(self, R1, R2, q=0.999, nr=64, rstar=None):
        """Predictions for the annulus R1 <= |z| <= R2.

        If rstar is given (a lattice distance), the bias pattern is computed at that radius instead of
        the model's own r* (conditional prediction)."""
        kappa = self.ring.kappa
        t = (np.arange(nr) + 0.5) / nr
        r = np.sqrt(R1 ** 2 + t * (R2 ** 2 - R1 ** 2))             # area-uniform radii
        logN = np.log(r ** 2)
        # P(h >= shell | m), averaged over N
        P = np.exp(-kappa * self.lam[:, :, None] / logN[None, None, :]).mean(axis=2)   # (#m, #shells)
        dens = self.dens.copy()
        dens[0] *= 1 - float(np.mean(kappa / logN)) * self.prime_factor
        dens /= dens.sum()
        tail = dens @ P                                             # P(h >= shell)
        if rstar is None:
            k = int(np.max(np.where(tail >= 1 - q)[0]))
        else:
            k = int(np.argmin(np.abs(self.shells - rstar)))
        rstar = float(self.shells[k])
        w = P[:, k]
        tot = dens @ w
        by_omega = {}
        for j in range(int(self.omega.max()) + 1):
            s = self.omega == j
            if dens[s].sum() > 0:
                by_omega[j] = float((dens[s] @ w[s]) / dens[s].sum() / tot)
        by_ideal = []
        for j, (pa, pb, n) in enumerate(self.sp):
            s = self.bits[:, j]
            by_ideal.append({'norm': n, 'gen': [int(pa), int(pb)], 'share': float(dens[s].sum()),
                             'ratio': float((dens[s] @ w[s]) / dens[s].sum() / tot)})
        return {'rstar': rstar, 'tail_at_rstar': float(tail[k]), 'tail_next': float(tail[k + 1]),
                'class_share': {j: float(dens[self.omega == j].sum()) for j in by_omega},
                'by_omega': by_omega, 'by_ideal': by_ideal}

    def rstar_transitions(self, Rmin, Rmax, width, q=0.999, steps=400):
        """Radii R in [Rmin, Rmax] where the predicted r* of the annulus [R, R + width] changes."""
        Rs = np.exp(np.linspace(math.log(Rmin), math.log(Rmax), steps))
        vals = [self.annulus(R, R + width, q)['rstar'] for R in Rs]
        out = []
        for i in range(1, len(Rs)):
            if vals[i] != vals[i - 1]:
                lo, hi = Rs[i - 1], Rs[i]
                for _ in range(30):
                    mid = math.sqrt(lo * hi)
                    if self.annulus(mid, mid + width, q)['rstar'] == vals[i - 1]:
                        lo = mid
                    else:
                        hi = mid
                out.append({'R': float(hi), 'from': vals[i - 1], 'to': vals[i]})
        return vals[0], out


if __name__ == '__main__':
    for D in FIELDS:
        M = Model(D)
        print(D, len(M.sp), 'ideals', [M.annulus(R1, R2)['rstar'] for R1, R2 in ((500, 1100), (2000, 2600), (4000, 4600))])
