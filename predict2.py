"""Second model: Poisson model + pair correlations from the Hardy-Littlewood pair singular series.

For a signature m (small prime ideals of norm <= y dividing z) and the offsets u with 0 < |u| < r:
    p_u = kappa * S_m(u) / log N
    log P(h(z) >= r | m) = sum_u log(1 - p_u) + sum_{u<v} (p_uv - p_u p_v) / ((1 - p_u)(1 - p_v)),
    p_uv = p_u p_v R_m(u, v),
where R_m(u, v) = P(z+u, z+v both prime) / (P(z+u prime) P(z+v prime)) is the ratio of local factors:
  * small ideals (norm <= y), conditioned on whether they divide z;
  * all larger prime ideals (unconditioned, z uniform mod p):  (1 - nu/N) / (1 - 1/N)^2 with nu = 1 if
    p | u - v and nu = 2 otherwise (product taken up to norm X = 10^6; the tail is negligible).
This is the second-order (pair) term of the expansion of P(no prime in the disc); it has no free parameter.
"""
import math

import numpy as np

from holes import QF, kron, offsets, is_prime

Y = 31
X_LARGE = 10 ** 6


def large_ideal_constant(D, y, X):
    """L0 = sum over prime ideals y < N(p) <= X of log((1 - 2/N)/(1 - 1/N)^2)."""
    sv = np.ones(X + 1, bool); sv[:2] = False
    for i in range(2, int(X ** .5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    L0 = 0.0
    f = lambda N: math.log1p(-2 / N) - 2 * math.log1p(-1 / N)
    for p in np.nonzero(sv)[0]:
        p = int(p); c = kron(D, p)
        if c == 1:
            if p > y:
                L0 += 2 * f(p)
        elif c == 0:
            if p > y:
                L0 += f(p)
        else:
            if y < p * p <= X:
                L0 += f(p * p)
    return L0


class Model2:
    def __init__(self, D, rmax, y=Y):
        self.ring = ring = QF(D)
        self.sp = sp = ring.small_primes(y)
        s = len(sp)
        ua, ub, d = offsets(ring, rmax)
        o = np.argsort(d, kind='stable'); ua, ub, d = ua[o], ub[o], d[o]
        self.ua, self.ub, self.d = ua, ub, d
        rd = np.round(d, 9)
        self.shells = np.unique(rd)
        self.first = np.searchsorted(rd, self.shells, side='left')     # offsets [:first[k]] have |u| < shell k
        n = len(ua)
        N = np.array([nn for _, _, nn in sp], float)
        self.norms = N
        udiv = np.stack([ring.divides(pa, pb, nn, ua, ub) for pa, pb, nn in sp], 1)       # (n, s)
        masks = np.arange(1 << s)
        bits = ((masks[:, None] >> np.arange(s)[None, :]) & 1).astype(bool)
        self.bits, self.omega = bits, bits.sum(1)
        self.dens = np.prod(np.where(bits, 1 / N, 1 - 1 / N), axis=1)
        one, none = N / (N - 1), N * (N - 2) / (N - 1) ** 2
        f = np.where(bits[:, None, :] & udiv[None], 0.0,
                     np.where(bits[:, None, :] ^ udiv[None], one, none))
        self.S = f.prod(axis=2)                                                             # (#m, n)
        # pairs i < j
        I, J = np.triu_indices(n, 1)
        self.I, self.J = I, J
        wa, wb = ua[I] - ua[J], ub[I] - ub[J]
        # small ideals: log A_j(pair), used when the ideal does NOT divide z
        logA = np.zeros((s, len(I)))
        for j, (pa, pb, nn) in enumerate(sp):
            du, dv = udiv[I, j], udiv[J, j]
            dw = ring.divides(pa, pb, nn, wa, wb)
            nu = np.where(du & dv, 0, np.where(du ^ dv, 1, np.where(dw, 1, 2)))
            fu = nn - 1 - (~du); fv = nn - 1 - (~dv)
            ok = (fu > 0) & (fv > 0)
            A = np.where(ok, (nn - 1 - nu) * (nn - 1) / np.where(ok, fu * fv, 1), 1.0)
            logA[j] = np.log(np.maximum(A, 1e-300))
        self.logA = logA
        # large ideals: log G(w) = L0 + sum over p | w with N(p) > y of log((1 - 1/N)/(1 - 2/N))
        L0 = large_ideal_constant(D, y, X_LARGE)
        Nw = ring.norm(wa, wb)
        big = [x for x in ring.small_primes(int(Nw.max()) + 1) if x[2] > y]
        logG = np.full(len(I), L0)
        for pa, pb, nn in big:
            dw = ring.divides(pa, pb, nn, wa, wb)
            logG += np.where(dw, math.log1p(-1 / nn) - math.log1p(-2 / nn), 0.0)
        self.logG = logG
        self.prime_factor = float(np.prod(N / (N - 1)))

    def log_tail(self, logN, block=256):
        """log P(h >= shell | m) for every mask and shell, at a single value of log N."""
        kappa = self.ring.kappa
        p = np.minimum(kappa * self.S / logN, 0.99)                                  # (#m, n)
        lp = np.log1p(-p)
        first_order = np.concatenate([np.zeros((len(p), 1)), np.cumsum(lp, axis=1)], axis=1)[:, self.first]
        nm, n = p.shape
        corr = np.zeros((nm, len(self.shells)))
        notbits = (~self.bits).astype(float)
        for a in range(0, nm, block):
            R = np.exp(notbits[a:a + block] @ self.logA + self.logG[None, :])          # (b, #pairs)
            pi, pj = p[a:a + block][:, self.I], p[a:a + block][:, self.J]
            T = pi * pj * (R - 1) / ((1 - pi) * (1 - pj))
            # pair (i, j) with i < j counts for shells k with first[k] > j
            byj = np.stack([np.bincount(self.J, weights=t, minlength=n) for t in T])   # (b, n)
            cum = np.concatenate([np.zeros((len(byj), 1)), np.cumsum(byj, axis=1)], axis=1)
            corr[a:a + block] = cum[:, self.first]
        return first_order + corr, first_order

    def annulus(self, R1, R2, q=0.999, rstar=None, pairs=True, nr=3):
        kappa = self.ring.kappa
        t = (np.arange(nr) + 0.5) / nr
        rr = np.sqrt(R1 ** 2 + t * (R2 ** 2 - R1 ** 2))
        logNs = np.log(rr ** 2)
        Ps = []
        for L in logNs:
            full, first = self.log_tail(L)
            Ps.append(np.exp(full if pairs else first))
        P = np.mean(Ps, axis=0)                                                          # (#m, #shells)
        dens = self.dens.copy()
        dens[0] *= 1 - float(np.mean(kappa / logNs)) * self.prime_factor
        dens /= dens.sum()
        tail = dens @ P
        if rstar is None:
            k = int(np.max(np.where(tail >= 1 - q)[0]))
        else:
            k = int(np.argmin(np.abs(self.shells - rstar)))
        w = P[:, k]; tot = dens @ w
        by_omega = {}
        for j in range(int(self.omega.max()) + 1):
            s = self.omega == j
            if dens[s].sum() > 0:
                by_omega[j] = float((dens[s] @ w[s]) / dens[s].sum() / tot)
        by_ideal = []
        for j, (pa, pb, n) in enumerate(self.sp):
            s = self.bits[:, j]
            by_ideal.append({'norm': n, 'gen': [int(pa), int(pb)],
                             'ratio': float((dens[s] @ w[s]) / dens[s].sum() / tot)})
        return {'rstar': float(self.shells[k]), 'tail': tail, 'shells': self.shells,
                'class_share': {j: float(dens[self.omega == j].sum()) for j in by_omega},
                'by_omega': by_omega, 'by_ideal': by_ideal}
