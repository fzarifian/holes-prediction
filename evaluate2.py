"""Apply the pre-registered criteria of PREREGISTRATION_v2.md to verification_v2/measurements_v2.json."""
import json
import math
import os

import numpy as np
from scipy.stats import binomtest

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_EXPECTED = 100      # expected number of points of the class among the points with h >= r*
CLEAR = 0.2             # |log predicted ratio| above which a prediction counts as clear
MIN_PAIRS = 10


def lookup(cond, rstar):
    for k, v in cond.items():
        if abs(float(k) - rstar) < 1e-6:
            return v
    return None


def verdict(ok, n):
    if n < MIN_PAIRS:
        return f'inconclusive ({n} eligible pairs < {MIN_PAIRS})'
    p = binomtest(ok, n, 0.5, alternative='greater').pvalue
    good = ok / n >= 0.75 and p < 0.01
    return f'{"SUCCESS" if good else "FAILURE"}: {ok}/{n} = {ok / n:.0%} correct side, one-sided binomial p = {p:.1e}'


def main():
    pred = json.load(open(os.path.join(HERE, 'predictions_v2', 'predictions_v2.json')))
    meas = json.load(open(os.path.join(HERE, 'verification_v2', 'measurements_v2.json')))
    M = {(m['D'], m['R1']): m for m in meas}
    p1 = [0, 0, 0, 0]; p3 = [0, 0]; exact = within = ntot = 0; flips = []
    for D, f in pred['fields'].items():
        D = int(D); uncond = {}
        for c in f['tests']:
            m = M[(D, c['R1'])]; ntot += 1
            shells = sorted(float(k) for k in c['conditional'])
            kp = int(np.argmin([abs(s - c['rstar_pred']) for s in shells]))
            km = int(np.argmin([abs(s - m['rstar']) for s in shells]))
            hit = abs(shells[km] - m['rstar']) < 1e-6
            exact += hit and km == kp; within += hit and abs(km - kp) <= 1
            pat = lookup(c['conditional'], m['rstar'])
            if pat is None:
                print(f'D={D} R1={c["R1"]}: measured r* = {m["rstar"]:.4f} outside the prediction table')
                continue
            for k, s in m['by_omega'].items():
                pr = pat['by_omega'].get(k)
                if pr is None or s['expected_rec'] < MIN_EXPECTED or abs(math.log(pr)) <= CLEAR:
                    continue
                ok = (s['obs'] - 1) * (pr - 1) > 0
                p1[0] += ok; p1[1] += 1
                if not (s['ci95'][0] <= 1 <= s['ci95'][1]):
                    p1[2] += ok; p1[3] += 1
            for s in m['by_ideal']:
                pr = next((d['ratio'] for d in pat['by_ideal'] if d['norm'] == s['norm'] and d['gen'] == s['gen']), None)
                if pr is None or s['expected_rec'] < MIN_EXPECTED or abs(math.log(pr)) <= CLEAR:
                    continue
                p3[0] += (s['obs'] - 1) * (pr - 1) > 0; p3[1] += 1
            up = lookup(c['conditional'], c['rstar_pred'])
            uncond[c['R1']] = (up, m)
        radii = sorted(uncond)
        for a, b in zip(radii, radii[1:]):
            (pa, ma), (pb, mb) = uncond[a], uncond[b]
            for k in pa['by_omega']:
                ra, rb = pa['by_omega'][k], pb['by_omega'].get(k)
                if rb is None or (ra - 1) * (rb - 1) >= 0 or min(abs(math.log(ra)), abs(math.log(rb))) <= 0.1:
                    continue
                sa, sb = ma['by_omega'].get(k), mb['by_omega'].get(k)
                if not sa or not sb or min(sa['expected_rec'], sb['expected_rec']) < MIN_EXPECTED:
                    continue
                flips.append((D, k, a, b, ra, rb, sa['obs'], sb['obs'],
                              (sa['obs'] - 1) * (sb['obs'] - 1) < 0 and (sa['obs'] - 1) * (ra - 1) > 0))
    print('P1 (co-primary, omega classes, conditional on measured r*):', verdict(p1[0], p1[1]))
    print(f'    among pairs whose bootstrap interval excludes 1: {p1[2]}/{p1[3]}')
    print(f'P2 (co-primary, record radius): exact {exact}/{ntot}, within one shell {within}/{ntot} ->',
          'SUCCESS' if exact >= ntot / 2 and within >= 0.85 * ntot else 'FAILURE')
    print('P3 (ideals, conditional on measured r*):', verdict(p3[0], p3[1]))
    nobs = sum(f[-1] for f in flips)
    if len(flips) < 3:
        v4 = f'inconclusive ({len(flips)} eligible reversals < 3)'
    else:
        v4 = 'SUCCESS' if nobs >= 2 * len(flips) / 3 else 'FAILURE'
    print(f'P4 (predicted reversals between test radii, at predicted r*): {nobs}/{len(flips)} observed -> {v4}')
    for f in flips:
        print(f'    D={f[0]} omega={f[1]}: R={f[2]} pred {f[4]:.2f} obs {f[6]:.2f} | R={f[3]} pred {f[5]:.2f} obs {f[7]:.2f}'
              f' -> {"observed" if f[-1] else "not observed"}')


if __name__ == '__main__':
    main()
