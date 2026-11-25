"""Hole statistics for prime elements in the nine imaginary quadratic fields of class number one.

Ring of integers O = Z[alpha]:
    alpha = (1 + sqrt(D)) / 2  if D = 1 mod 4   (norm a^2 + a b + m b^2,  m = (1 - D) / 4)
    alpha = sqrt(D / 4)        if D = 0 mod 4   (norm a^2 + k b^2,        k = -D / 4)
An element a + b*alpha is stored as the integer pair (a, b).

All analyses are restricted to a fundamental sector of the unit/conjugation symmetry
(30 deg for D = -3, 45 deg for D = -4, 90 deg otherwise), so that each hole is counted once.
"""
import math
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import chi2

FIELDS = (-3, -4, -7, -8, -11, -19, -43, -67, -163)
ANNULI = ((500, 1100), (2000, 2600), (4000, 4600))


# ----------------------------------------------------------------------------- arithmetic

def is_prime(n):
    if n < 2:
        return False
    for d in range(2, int(n ** .5) + 1):
        if n % d == 0:
            return False
    return True


def kron(D, p):
    """Kronecker symbol (D / p) for a rational prime p."""
    if p == 2:
        if D % 2 == 0:
            return 0
        return 1 if D % 8 in (1, 7) else -1
    if D % p == 0:
        return 0
    return 1 if pow(D % p, (p - 1) // 2, p) == 1 else -1


class QF:
    """Ring of integers of Q(sqrt(D)), D one of the nine class-number-one discriminants."""

    def __init__(self, D):
        self.D = D
        if D % 4 == 0:
            self.k = -D // 4; self.t = 0.0; self.s = math.sqrt(self.k); self.half = False
        else:
            self.m = (1 - D) // 4; self.t = 0.5; self.s = math.sqrt(-D) / 2; self.half = True
        if D == -3:
            self.units = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, 1), (1, -1)]; self.theta = math.pi / 6
        elif D == -4:
            self.units = [(1, 0), (-1, 0), (0, 1), (0, -1)]; self.theta = math.pi / 4
        else:
            self.units = [(1, 0), (-1, 0)]; self.theta = math.pi / 2
        # density of prime elements near z is ~ kappa / log N(z)  (prime ideal theorem, units included)
        self.kappa = len(self.units) * math.sqrt(-D) / (2 * math.pi)

    def xy(self, a, b):
        return a + b * self.t, b * self.s

    def norm(self, a, b):
        return a * a + a * b + self.m * b * b if self.half else a * a + self.k * b * b

    def mul(self, a, b, c, d):
        if self.half:
            return a * c - self.m * b * d, a * d + b * c + b * d
        return a * c - self.k * b * d, a * d + b * c

    def conj(self, a, b):
        return (a + b, -b) if self.half else (a, -b)

    def divides(self, pa, pb, n, za, zb):
        """pi | z  <=>  z * conj(pi) = 0 mod N(pi)   (pi = pa + pb*alpha, N(pi) = n)."""
        ca, cb = self.conj(pa, pb)
        ra, rb = self.mul(za, zb, ca, cb)
        return (ra % n == 0) & (rb % n == 0)

    def small_primes(self, y):
        """One generator (a, b, N) for each prime ideal of norm <= y."""
        out = []
        for p in range(2, y + 1):
            if not is_prime(p):
                continue
            c = kron(self.D, p)
            if c == -1:                                   # inert: ideal (p) of norm p^2
                if p * p <= y:
                    out.append((p, 0, p * p))
                continue
            found = None
            B = int(math.sqrt(4 * p / (-self.D))) + 2
            for b in range(0, B + 1):
                for a in range(-p - b, p + 1):
                    if self.norm(a, b) == p:
                        found = (a, b); break
                if found:
                    break
            a, b = found
            out.append((a, b, p))
            if c == 1:                                    # split: add the conjugate ideal
                ca, cb = self.conj(a, b)
                out.append((ca, cb, p))
        return out


# ----------------------------------------------------------------------------- lattice and primes

def sector_points(ring, R1, R2, margin=12):
    """Lattice points within distance `margin` of the sector {R1 <= |z| <= R2, 0 <= arg z <= theta}."""
    Ro, Ri = R2 + margin, max(0, R1 - margin)
    ymax = Ro * math.sin(ring.theta) + margin + 2
    bs = np.arange(int(math.floor(-(margin + 2) / ring.s)), int(math.ceil(ymax / ring.s)) + 1)
    A, B = [], []
    for b in bs:
        Y = b * ring.s
        if abs(Y) > Ro:
            continue
        xmax = math.sqrt(Ro * Ro - Y * Y); xmin = -margin
        aa = np.arange(math.ceil(xmin - b * ring.t), math.floor(xmax - b * ring.t) + 1)
        A.append(aa); B.append(np.full(len(aa), b))
    A = np.concatenate(A).astype(np.int64); B = np.concatenate(B).astype(np.int64)
    X, Y = ring.xy(A, B)
    r = np.hypot(X, Y); ang = np.arctan2(Y, X)
    ok = (r >= Ri) & (r <= Ro)
    dist_lo = np.where(ang < 0, r * np.sin(np.minimum(-ang, math.pi / 2)), 0)
    dist_hi = np.where(ang > ring.theta, r * np.sin(np.minimum(ang - ring.theta, math.pi / 2)), 0)
    ok &= (dist_lo <= margin) & (dist_hi <= margin)
    return A[ok], B[ok], X[ok], Y[ok], r[ok], ang[ok]


def primality(ring, A, B, Nmax):
    """z is prime iff N(z) is a rational prime, or N(z) = q^2 with q an inert rational prime."""
    sv = np.ones(Nmax + 1, bool); sv[:2] = False
    for i in range(2, int(Nmax ** .5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    N = ring.norm(A, B)
    q = np.rint(np.sqrt(N)).astype(np.int64)
    P = np.zeros(len(N), bool)
    ok = (N >= 2) & (N <= Nmax); P[ok] = sv[N[ok]]
    sq = (q * q == N) & (q >= 2) & (q <= Nmax)
    sq &= sv[np.clip(q, 0, Nmax)]
    inert = {int(x) for x in np.unique(q[sq]) if kron(ring.D, int(x)) == -1}
    P |= sq & np.isin(q, list(inert))
    return P, N


class Config:
    """All non-prime lattice points of one annulus (fundamental sector) with their hole radius h(z)."""

    def __init__(self, ring, R1, R2):
        self.ring, self.R1, self.R2 = ring, R1, R2
        A, B, X, Y, r, ang = sector_points(ring, R1, R2)
        P, N = primality(ring, A, B, int((R2 + 14) ** 2) + 10)
        self.A, self.B, self.P = A, B, P
        self.tree = cKDTree(np.c_[X[P], Y[P]])
        self.prime_share = float(P[(r >= R1) & (r <= R2)].mean())
        core = (r >= R1) & (r <= R2) & (ang >= 0) & (ang <= ring.theta) & (~P)
        self.Ac, self.Bc, self.Xc, self.Yc, self.Nc = A[core], B[core], X[core], Y[core], N[core]
        self.h, _ = self.tree.query(np.c_[self.Xc, self.Yc])
        self.n = int(core.sum())

    def records(self, K=150, sep=10.0):
        """Indices of the K largest h, greedily keeping points more than `sep` apart (one per hole)."""
        order = np.argsort(-self.h); chosen = []; pts = np.c_[self.Xc, self.Yc]
        for i in order[:300000]:
            if len(chosen) >= K:
                break
            if all((pts[i, 0] - pts[j, 0]) ** 2 + (pts[i, 1] - pts[j, 1]) ** 2 > sep * sep for j in chosen):
                chosen.append(i)
        return np.array(chosen)

    def signature(self, sp):
        """Bit mask of the small prime ideals (list sp) dividing each point."""
        mask = np.zeros(len(self.Ac), np.int64)
        for j, (pa, pb, n) in enumerate(sp):
            mask |= self.ring.divides(pa, pb, n, self.Ac, self.Bc).astype(np.int64) << j
        return mask


def omega_of(mask, nbits):
    return np.array([bin(int(m)).count('1') for m in range(1 << nbits)])[mask] if nbits <= 16 else \
        np.array([bin(int(m)).count('1') for m in mask])


# ----------------------------------------------------------------------------- statistics

def ratio_table(omega, chosen):
    """Over-representation of each omega class among the records, with a chi-square test."""
    base = np.bincount(omega, minlength=8)[:8] / len(omega)
    obs = np.bincount(omega[chosen], minlength=8)[:8].astype(float)
    exp = base * obs.sum()
    m = exp > 5
    x2 = ((obs[m] - exp[m]) ** 2 / exp[m]).sum()
    return {'obs': obs.astype(int).tolist(),
            'ratio': [round(float(o / e), 2) if e > 0.5 else None for o, e in zip(obs, exp)],
            'p': float(chi2.sf(x2, max(1, m.sum() - 1)))}


def neighbour_rate(cfg, omega, step=15):
    """Share of unit neighbours z + u that are prime, by omega class of z (every `step`-th point)."""
    key = set(zip(cfg.A[cfg.P].tolist(), cfg.B[cfg.P].tolist()))
    out = {}
    idx = np.arange(0, len(cfg.Ac), step)
    for k in range(0, 6):
        s = idx[omega[idx] == k]
        if len(s) < 200:
            continue
        c = t = 0
        for i in s:
            for du, dv in cfg.ring.units:
                t += 1; c += ((int(cfg.Ac[i]) + du, int(cfg.Bc[i]) + dv) in key)
        out[k] = round(c / t, 4)
    return out


# ----------------------------------------------------------------------------- singular series model

def offsets(ring, rmax):
    """Lattice offsets u with 0 < |u| < rmax, and their lengths."""
    R = int(rmax * 2) + 4
    ua, ub = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1)); ua = ua.ravel(); ub = ub.ravel()
    ux, uy = ring.xy(ua, ub); d = np.hypot(ux, uy)
    s = (d > 0) & (d < rmax)
    return ua[s], ub[s], d[s]


def lam(ring, sp, zdiv, ua, ub, d, r):
    """Lambda_D(r) = sum over 0 < |u| < r of the conditional singular series S_D(u).

    zdiv[j] is True when the j-th small prime ideal divides z.  Local factors:
        0                 if pi | z and pi | u
        N/(N-1)           if pi divides exactly one of z, u
        N(N-2)/(N-1)^2    if pi divides neither.
    """
    sel = d < r
    S = np.ones(sel.sum())
    for j, (pa, pb, n) in enumerate(sp):
        ud = ring.divides(pa, pb, n, ua[sel], ub[sel])
        both = zdiv[j] & ud; one = zdiv[j] ^ ud
        S *= np.where(both, 0.0, np.where(one, n / (n - 1), n * (n - 2) / (n - 1) ** 2))
    return S.sum(), int(sel.sum())


def poisson_weights(cfg, sp, mask, r, ua, ub, d):
    """w(z) = exp(-kappa Lambda_{D(z)}(r) / log N(z)): model probability that h(z) >= r (up to a constant)."""
    Lm = {}
    for m in np.unique(mask):
        zd = np.array([(int(m) >> j) & 1 for j in range(len(sp))], bool)
        Lm[int(m)] = lam(cfg.ring, sp, zd, ua, ub, d, r)[0]
    L = np.array([Lm[int(m)] for m in mask])
    return np.exp(-cfg.ring.kappa * L / np.log(cfg.Nc))


def model_by_omega(cfg, sp, mask, omega, q=0.999):
    """Observed vs predicted over-representation of each omega class among points with h >= r* (quantile q)."""
    rstar = float(np.quantile(cfg.h, q))
    ua, ub, d = offsets(cfg.ring, rstar + 1.0)
    w = poisson_weights(cfg, sp, mask, rstar, ua, ub, d)
    rec = cfg.h >= rstar
    out = {'rstar': rstar}
    for k in range(0, 6):
        s = omega == k
        if s.sum() < 500:
            continue
        out[k] = (round(float(rec[s].mean() / rec.mean()), 2), round(float(w[s].mean() / w.mean()), 2))
    return out


# ----------------------------------------------------------------------------- one configuration

def analyse(ring, R1, R2, ys=(13, 31, 100), y_model=31, K=150, sep=10.0):
    cfg = Config(ring, R1, R2)
    chosen = cfg.records(K, sep)
    res = {'D': ring.D, 'R': [R1, R2], 'n_points': cfg.n, 'prime_share': cfg.prime_share,
           'record_h': [float(cfg.h[chosen].min()), float(cfg.h[chosen].max())]}
    for y in ys:
        sp = ring.small_primes(y)
        omega = np.zeros(len(cfg.Ac), int)
        for pa, pb, n in sp:
            omega += ring.divides(pa, pb, n, cfg.Ac, cfg.Bc)
        res[f'records_y{y}'] = ratio_table(omega, chosen)
        if y in (13, 31):
            res[f'neigh_y{y}'] = neighbour_rate(cfg, omega)
    sp = ring.small_primes(y_model)
    mask = cfg.signature(sp)
    res['model_y31'] = model_by_omega(cfg, sp, mask, omega_of(mask, len(sp)))
    res['small_prime_norms'] = [n for _, _, n in sp]
    return res, cfg, chosen
