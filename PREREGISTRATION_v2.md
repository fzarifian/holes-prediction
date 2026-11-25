# Second pre-registered prediction: holes between prime elements at radii 2×10^5 to 10^6

Author: Fabien Zarifian. Date of registration: 25 November 2026.

This is the second test of the singular-series model for the holes between prime elements of the nine
imaginary quadratic fields of class number one. The first test (release `v1.0-prediction`, results in
`RESULTS.md`) confirmed the predicted direction of the bias at the measured record radius, but failed to
predict the record radius r* itself, and with it the scales where the bias reverses.

As before, this document fixes the predictions, the measurement code and the success criteria before any
measurement. Nothing about holes or primes at radii above 10^5 had been measured when these files were
frozen; the only computation run at those radii was a timing test of point generation and sieving at
radius 10^6 (point counts only, no statistics).

## 1. What changed in the model

The first model treated the events "z + u is prime" as independent. The second model adds the pair term
of the Hardy–Littlewood conjectures (`predict2.py`):

    log P(h(z) >= r | signature) = sum_u log(1 - p_u) + sum_(u<v) (p_uv - p_u p_v) / ((1 - p_u)(1 - p_v)),

with p_u = kappa S(u) / log N and p_uv = p_u p_v R(u, v), where R(u, v) is the ratio of the local factors
of the pair singular series (small prime ideals conditioned on whether they divide z, all larger prime
ideals up to norm 10^6 unconditioned).

**No parameter is fitted.** The calibration constant of the first test (`C_RSTAR = 1.25`) is dropped.

## 2. Status of the earlier data

The form of the correction was chosen after the failure of the first test, so the 54 configurations
already measured (radii 500 to 10^5) are development data, not evidence. On them, with no adjustment:

* record radius r* exact in 41 of 54 configurations and within one lattice shell in 53 of 54
  (first model with its fitted constant: 25 of 54 exact);
* direction of the bias at the measured r* correct in 104 of 114 clear predictions.

The code was run once in this form; no option was tuned to improve these numbers.

## 3. Test configurations (fixed)

For each of the 9 fields, three annuli R <= |z| <= R + 100 with R in {200000, 450000, 1000000}, inside the
angular window 0 <= arg z <= φ, with φ chosen to give about 4 million lattice points (27 configurations;
values of φ below and in `predictions_v2/predictions_v2.json`). The window is a sub-sector of the
fundamental sector; it keeps the computation small at radius 10^6.

## 4. Measured quantities (fixed; code `verify2.py`)

As in the first test: non-prime lattice points of the window; h(z); r* = 99.9% quantile of h; for each
omega class with at least 500 points and for each prime ideal of norm <= 31 (conjugates merged), the share
among points with h >= r* divided by the overall share; bootstrap over holes (2000 draws, seed 12345).
`verify2.py check` reproduces exactly four configurations measured earlier.

## 5. Criteria (fixed; code `evaluate2.py`)

Eligibility as in the first test: expected number of class points among the points with h >= r* at
least 100, and clear prediction (|log predicted ratio| > 0.2).

* **P1 (co-primary): direction of the bias by omega class, at the measured r*.** Success if at least 75% of
  eligible pairs are on the correct side of 1 and the one-sided binomial test against 50% gives p < 0.01
  (inconclusive below 10 eligible pairs).
* **P2 (co-primary): record radius.** Success if the measured r* equals the predicted r* in at least half of
  the 27 configurations and is within one lattice shell in at least 85%.
* **P3: direction of the bias by prime ideal**, at the measured r*: same rule as P1.
* **P4: reversals of the bias between consecutive test radii, at the predicted r*.** A reversal is eligible
  if the predicted ratios lie on opposite sides of 1 with |log ratio| > 0.1 at both radii and the class
  has an expected count of at least 100 at both. Success if at least two thirds of the eligible reversals
  are observed (inconclusive below 3).

Any deviation from this protocol will be reported together with the result.

## 6. Predictions

Predicted record radius and its changes between 10^5 and 10^6 (annuli of width 100):

| D | predicted r* at R = 10^5 | predicted changes of r* (R: new r*) |
|---|---|---|
| -3 | 4.583 | 110872: 5.000, 219411: 5.196, 397593: 5.292 |
| -4 | 5.831 | 103982: 6.000, 122448: 6.083, 309568: 6.325, 553428: 6.403 |
| -7 | 8.185 | 105169: 8.426, 124515: 8.485, 196253: 8.602, 344164: 8.775, 444196: 8.888, 526884: 9.000, 556897: 9.274, 949935: 9.381 |
| -8 | 8.246 | 113229: 8.485, 142497: 8.544, 188201: 8.660, 259046: 8.718, 337059: 9.000, 671770: 9.055, 905402: 9.110 |
| -11 | 8.307 | 107998: 8.426, 136335: 8.660, 337512: 8.944, 459384: 9.000 |
| -19 | 8.544 | 108090: 8.718, 125238: 8.775, 275107: 8.944, 426273: 9.000, 531951: 9.110, 794806: 9.220 |
| -43 | 8.246 | 121514: 8.888, 207837: 9.000, 275997: 9.110, 542498: 9.592, 937030: 9.849 |
| -67 | 8.544 | 173692: 8.718, 367580: 9.000, 558740: 9.110 |
| -163 | 8.000 | 143734: 8.426, 407415: 9.000, 763085: 9.110 |

Predicted over-representation among points with h >= r*, at the predicted r* (full conditional tables and
all prime ideals in `predictions_v2/predictions_v2.json`):

| D | annulus (width 100) | window φ (rad) | predicted r* | predicted over-representation among points with h >= r* (ω classes with share >= 2.5%) |
|---|---|---|---|---|
| -3 | [200000, 200100] | 0.1732 | 5.000 | ω=0: 1.07, ω=1: 1.00, ω=2: 0.93, ω=3: 0.98 |
| -3 | [450000, 450100] | 0.0770 | 5.292 | ω=0: 1.05, ω=1: 1.03, ω=2: 0.93, ω=3: 0.92 |
| -3 | [1000000, 1000100] | 0.0346 | 5.292 | ω=0: 1.05, ω=1: 1.03, ω=2: 0.93, ω=3: 0.92 |
| -4 | [200000, 200100] | 0.2000 | 6.083 | ω=0: 0.75, ω=1: 0.88, ω=2: 1.10, ω=3: 1.35 |
| -4 | [450000, 450100] | 0.0889 | 6.325 | ω=0: 1.36, ω=1: 0.98, ω=2: 0.88, ω=3: 0.86 |
| -4 | [1000000, 1000100] | 0.0400 | 6.403 | ω=0: 0.72, ω=1: 0.99, ω=2: 1.09, ω=3: 1.14 |
| -7 | [200000, 200100] | 0.2646 | 8.602 | ω=0: 0.75, ω=1: 1.03, ω=2: 1.05, ω=3: 1.00, ω=4: 0.93 |
| -7 | [450000, 450100] | 0.1176 | 8.888 | ω=0: 1.26, ω=1: 1.05, ω=2: 0.90, ω=3: 0.93, ω=4: 0.97 |
| -7 | [1000000, 1000100] | 0.0529 | 9.381 | ω=0: 1.92, ω=1: 1.14, ω=2: 0.77, ω=3: 0.68, ω=4: 0.64 |
| -8 | [200000, 200100] | 0.2828 | 8.660 | ω=0: 0.71, ω=1: 1.12, ω=2: 1.01, ω=3: 0.94, ω=4: 0.89 |
| -8 | [450000, 450100] | 0.1257 | 9.000 | ω=0: 0.84, ω=1: 1.11, ω=2: 0.96, ω=3: 0.96, ω=4: 0.97 |
| -8 | [1000000, 1000100] | 0.0566 | 9.110 | ω=0: 1.26, ω=1: 0.92, ω=2: 0.91, ω=3: 1.13, ω=4: 1.26 |
| -11 | [200000, 200100] | 0.3317 | 8.660 | ω=0: 1.48, ω=1: 1.08, ω=2: 0.90, ω=3: 0.74, ω=4: 0.59 |
| -11 | [450000, 450100] | 0.1474 | 8.944 | ω=0: 1.21, ω=1: 0.96, ω=2: 0.99, ω=3: 0.98, ω=4: 0.88 |
| -11 | [1000000, 1000100] | 0.0663 | 9.000 | ω=0: 1.22, ω=1: 0.97, ω=2: 0.95, ω=3: 0.98, ω=4: 0.95 |
| -19 | [200000, 200100] | 0.4359 | 8.775 | ω=0: 1.35, ω=1: 1.04, ω=2: 0.88, ω=3: 0.82, ω=4: 0.81 |
| -19 | [450000, 450100] | 0.1937 | 9.000 | ω=0: 1.21, ω=1: 1.00, ω=2: 0.92, ω=3: 0.92, ω=4: 1.01 |
| -19 | [1000000, 1000100] | 0.0872 | 9.220 | ω=0: 1.31, ω=1: 1.03, ω=2: 0.89, ω=3: 0.83, ω=4: 0.84 |
| -43 | [200000, 200100] | 0.6557 | 8.888 | ω=0: 0.98, ω=1: 1.06, ω=2: 0.97, ω=3: 0.84 |
| -43 | [450000, 450100] | 0.2914 | 9.110 | ω=0: 1.10, ω=1: 1.03, ω=2: 0.89, ω=3: 0.75 |
| -43 | [1000000, 1000100] | 0.1311 | 9.849 | ω=0: 1.01, ω=1: 1.06, ω=2: 0.95, ω=3: 0.80 |
| -67 | [200000, 200100] | 0.8185 | 8.718 | ω=0: 1.19, ω=1: 0.95, ω=2: 0.78, ω=3: 0.66 |
| -67 | [450000, 450100] | 0.3638 | 9.000 | ω=0: 1.01, ω=1: 1.02, ω=2: 0.95, ω=3: 0.87 |
| -67 | [1000000, 1000100] | 0.1637 | 9.110 | ω=0: 1.01, ω=1: 1.01, ω=2: 0.96, ω=3: 0.91 |
| -163 | [200000, 200100] | 1.2767 | 8.426 | ω=0: 1.10, ω=1: 0.89, ω=2: 0.71 |
| -163 | [450000, 450100] | 0.5674 | 9.000 | ω=0: 1.16, ω=1: 0.82, ω=2: 0.56 |
| -163 | [1000000, 1000100] | 0.2553 | 9.110 | ω=0: 1.10, ω=1: 0.89, ω=2: 0.72 |

Most striking predicted reversals: for D = −4, centres without small factor go from ×0.75 (R = 2×10^5)
to ×1.36 (R = 4.5×10^5) and back to ×0.72 (R = 10^6); for D = −7, from ×0.75 to ×1.26 and then ×1.92;
for D = −8, from ×0.84 (R = 4.5×10^5) to ×1.26 (R = 10^6).

## 7. Files and SHA-256 fingerprints

    a58bd09ac6fbbcced38ac099b4722e1336945904cec9837d296685ade81500cc  holes.py
    d00ad996dbf4215fc6052b7274097afe1986dcedf29e2aaaf5ba1476fb1a6d18  predict2.py
    11247d56fc92f8c09b3f54decee8d45e5414597596966d8a6cc282a67430a3d2  make_predictions2.py
    5cb1261c70b6609df4dfb33bb77b8d13da498472c443ecf4256ac87d30c2c3d3  verify2.py
    328a3773af540e19e4af22821279a896b637712b3548178083c83a3c772a664e  evaluate2.py
    4476c87b50992ef86c49e9ce4431dd817d5a8d25259b5305ba1257130801b28d  predictions_v2/predictions_v2.json
