"""Post hoc diagnostic (not pre-registered): model vs measured tail probability P(h >= d) near r*,
for the 27 configurations of the second test.  Writes verification_v2/tail_diagnostics.json."""
import json, math, numpy as np, sys
sys.path.insert(0, '.')
from scipy.spatial import cKDTree
from holes import QF
from verify2 import annulus_points, primality
from predict2 import Model2
from make_predictions2 import Interp, rmax_for, WIDTH
p = json.load(open('predictions_v2/predictions_v2.json'))
out = []
for D, f in p['fields'].items():
    D = int(D); M = Model2(D, rmax_for(D)); I = Interp(M, 1.5e5, 1.05e6)
    for c in f['tests']:
        ring = QF(D)
        A, B, X, Y, r, ang = annulus_points(ring, c['R1'], c['R2'], c['phi'])
        P, N = primality(ring, A, B)
        tree = cKDTree(np.c_[X[P], Y[P]])
        core = (r >= c['R1']) & (r <= c['R2']) & (ang >= 0) & (ang <= c['phi']) & (~P)
        h, _ = tree.query(np.c_[X[core], Y[core]])
        R = c['R1']; Pm = I.at(R + WIDTH / 2); logN = math.log((R + WIDTH / 2) ** 2)
        dens = M.dens.copy(); dens[0] *= 1 - M.ring.kappa / logN * M.prime_factor; dens /= dens.sum()
        tail = dens @ Pm
        k = int(np.argmin(abs(M.shells - c['rstar_pred'])))
        rows = []
        for kk in range(k - 1, k + 3):
            s = M.shells[kk]
            emp = float((h >= s - 1e-9).mean())
            rows.append((round(float(s), 3), float(tail[kk]), emp))
        meas = float(np.quantile(h, 0.999))
        out.append((D, R, c['rstar_pred'], meas, rows))
        print(D, R, 'pred %.3f meas %.3f' % (c['rstar_pred'], meas),
              ' '.join(f'[{s}: model {t*1e3:.2f} emp {e*1e3:.2f}]' for s, t, e in rows), flush=True)
json.dump(out, open('verification_v2/tail_diagnostics.json', 'w'))
