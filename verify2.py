"""Measurement code for the second pre-registered test (radii 2x10^5 to 10^6).  Frozen before any measurement there.

    python3 verify2.py check         # reproduces configurations already measured (radii 500 to 10^5)
    python3 verify2.py run           # measures the test annuli listed in predictions_v2/predictions_v2.json

Same statistics as verify.py (version 1).  New: an angular window 0 <= arg z <= phi (phi <= theta, the
fundamental sector), so that thin, long annuli at radii up to 10^6 stay small in memory; prime elements
are equidistributed in angle (Hecke), so the window is representative.

Differences with holes.Config: lattice points are generated in the annulus only (not the whole disc),
and primality uses a segmented sieve on the norms of the annulus, so that radii up to 10^5 are feasible.
The statistics are those of holes.model_by_omega: r* = 99.9% quantile of h over the non-prime points of
the fundamental sector, and for each class the share among points with h >= r* divided by its overall
share.  Uncertainty: points with h >= r* are grouped into holes (connected components at distance <= 2)
and holes are resampled (bootstrap, 2000 draws, fixed seed).
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

from holes import QF, kron

HERE = os.path.dirname(os.path.abspath(__file__))
Y = 31
MARGIN = 12


def annulus_points(ring, R1, R2, phi, margin=MARGIN):
    Ro, Ri = R2 + margin, max(0, R1 - margin)
    ymax = Ro * math.sin(phi) + margin + 2
    bs = np.arange(int(math.floor(-(margin + 2) / ring.s)), int(math.ceil(ymax / ring.s)) + 1)
    A, B = [], []
    for b in bs:
        Yv = b * ring.s
        if abs(Yv) > Ro:
            continue
        xo = math.sqrt(Ro * Ro - Yv * Yv)
        xi = math.sqrt(Ri * Ri - Yv * Yv) if abs(Yv) < Ri else 0.0
        for lo, hi in ((max(-margin, xi - 1), xo), (-xo, min(-xi + 1, margin)) if xi > 0 else (1, 0)):
            if hi < lo:
                continue
            aa = np.arange(math.ceil(lo - b * ring.t), math.floor(hi - b * ring.t) + 1)
            if len(aa):
                A.append(aa); B.append(np.full(len(aa), b))
    A = np.concatenate(A).astype(np.int64); B = np.concatenate(B).astype(np.int64)
    key = A * 4_000_003 + B                                    # drop duplicates from overlapping intervals
    _, idx = np.unique(key, return_index=True)
    A, B = A[idx], B[idx]
    X, Yc = ring.xy(A, B)
    r = np.hypot(X, Yc); ang = np.arctan2(Yc, X)
    ok = (r >= Ri) & (r <= Ro)
    dist_lo = np.where(ang < 0, r * np.sin(np.minimum(-ang, math.pi / 2)), 0)
    dist_hi = np.where(ang > phi, r * np.sin(np.minimum(ang - phi, math.pi / 2)), 0)
    ok &= (dist_lo <= margin) & (dist_hi <= margin)
    return A[ok], B[ok], X[ok], Yc[ok], r[ok], ang[ok]


def small_sieve(n):
    sv = np.ones(n + 1, bool); sv[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return sv


def segmented_is_prime(N, lo, hi):
    """Boolean array: N[i] is a rational prime (all N in [lo, hi])."""
    base = np.nonzero(small_sieve(int(math.isqrt(hi)) + 1))[0]
    seg = np.ones(hi - lo + 1, bool)
    for p in base:
        p = int(p)
        start = max(p * p, ((lo + p - 1) // p) * p)
        seg[start - lo::p] = False
    if lo <= 1:
        seg[:2 - lo] = False
    return seg[N - lo]


def primality(ring, A, B):
    N = ring.norm(A, B)
    lo, hi = int(N.min()), int(N.max())
    P = segmented_is_prime(N, lo, hi)
    q = np.rint(np.sqrt(N.astype(float))).astype(np.int64)
    sq = (q * q == N) & (q >= 2)
    if sq.any():
        qs = np.unique(q[sq])
        sv = small_sieve(int(qs.max()))
        inert = [int(x) for x in qs if sv[x] and kron(ring.D, int(x)) == -1]
        P |= sq & np.isin(q, inert)
    return P, N


def measure(D, R1, R2, phi=None, boot=2000, seed=12345):
    ring = QF(D)
    phi = ring.theta if phi is None else min(phi, ring.theta)
    A, B, X, Yc, r, ang = annulus_points(ring, R1, R2, phi)
    P, N = primality(ring, A, B)
    tree = cKDTree(np.c_[X[P], Yc[P]])
    core = (r >= R1) & (r <= R2) & (ang >= 0) & (ang <= phi) & (~P)
    Ac, Bc, Xc, Ycc, Nc = A[core], B[core], X[core], Yc[core], N[core]
    h, _ = tree.query(np.c_[Xc, Ycc])
    rstar = float(np.quantile(h, 0.999))
    rec = h >= rstar
    sp = ring.small_primes(Y)
    divs = np.stack([ring.divides(pa, pb, n, Ac, Bc) for pa, pb, n in sp], 1)
    omega = divs.sum(1)
    # holes: connected components of the points with h >= r* (link distance <= 2)
    idx = np.nonzero(rec)[0]
    pts = np.c_[Xc[idx], Ycc[idx]]
    parent = np.arange(len(idx))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i, j in cKDTree(pts).query_pairs(2.0):
        a, b = find(i), find(j)
        if a != b:
            parent[a] = b
    comp = np.array([find(i) for i in range(len(idx))])
    _, comp = np.unique(comp, return_inverse=True)
    ncomp = int(comp.max()) + 1
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, ncomp, size=(boot, ncomp))
    weights = np.stack([np.bincount(d, minlength=ncomp) for d in draws])          # (boot, ncomp)

    def stat(member):
        """observed ratio and bootstrap 95% interval for a boolean class over all core points."""
        share = member.mean()
        obs = member[idx].mean() / share
        cnt = np.bincount(comp, weights=member[idx].astype(float), minlength=ncomp)
        tot = np.bincount(comp, minlength=ncomp).astype(float)
        b = (weights @ cnt) / (weights @ tot) / share
        return {'share': float(share), 'n_points_class': int(member.sum()), 'n_rec_class': int(member[idx].sum()),
                'expected_rec': float(share * len(idx)), 'obs': float(obs),
                'ci95': [float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))],
                'frac_boot_above_1': float((b > 1).mean())}

    out = {'D': D, 'R1': R1, 'R2': R2, 'phi': phi, 'n_points': int(core.sum()), 'rstar': rstar,
           'n_rec_points': int(rec.sum()), 'n_holes': ncomp, 'by_omega': {}, 'by_ideal': []}
    for k in range(int(omega.max()) + 1):
        s = omega == k
        if s.sum() >= 500:
            out['by_omega'][str(k)] = stat(s)
    seen = set()
    for j, (pa, pb, n) in enumerate(sp):
        s = divs[:, j]
        key = (n, round(float(s.mean()), 4))
        if key in seen:
            continue
        seen.add(key)
        d = stat(s); d.update({'norm': n, 'gen': [int(pa), int(pb)]})
        out['by_ideal'].append(d)
    return out


def check():
    """Reproduce measured configurations: two of the paper (full sector) and two of test 1 (full sector)."""
    F = json.load(open(os.path.join(HERE, 'results', 'fields.json')))
    V = {(m['D'], m['R1']): m for m in json.load(open(os.path.join(HERE, 'verification', 'measurements.json')))}
    for D, R1, R2 in ((-3, 500, 1100), (-163, 2000, 2600)):
        m = measure(D, R1, R2, boot=200); ref = F[f'{D}_{R1}']
        print(f'D={D} [{R1},{R2}]: points equal {m["n_points"] == ref["n_points"]}, '
              f'r* equal {abs(m["rstar"] - ref["model_y31"]["rstar"]) < 1e-9}')
    for D, R1 in ((-7, 100000), (-67, 31623)):
        ref = V[(D, R1)]; m = measure(D, R1, ref['R2'], boot=2000)
        same = all(abs(m['by_omega'][k]['obs'] - v['obs']) < 1e-12 for k, v in ref['by_omega'].items())
        print(f'D={D} [{R1},{ref["R2"]}]: points equal {m["n_points"] == ref["n_points"]}, '
              f'r* equal {m["rstar"] == ref["rstar"]}, omega ratios equal {same}')


def run():
    pred = json.load(open(os.path.join(HERE, 'predictions_v2', 'predictions_v2.json')))
    os.makedirs(os.path.join(HERE, 'verification_v2'), exist_ok=True)
    res = []
    for D, f in pred['fields'].items():
        for c in f['tests']:
            t = time.time(); m = measure(int(D), c['R1'], c['R2'], c['phi'])
            res.append(m)
            print(f'D={D} [{c["R1"]},{c["R2"]}]: {m["n_points"]} points, r*={m["rstar"]:.3f}, '
                  f'{m["n_holes"]} holes ({time.time() - t:.0f} s)', flush=True)
    json.dump(res, open(os.path.join(HERE, 'verification_v2', 'measurements_v2.json'), 'w'), indent=1)


if __name__ == '__main__':
    {'check': check, 'run': run}[sys.argv[1]]()
