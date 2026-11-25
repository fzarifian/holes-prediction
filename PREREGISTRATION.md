# Pre-registered prediction: holes between prime elements at radii 10^4 to 10^5

Author: Fabien Zarifian. Date of registration: 25 November 2026.

This document fixes, before any measurement, (i) the predictions of a singular-series model for the
holes between prime elements of the nine imaginary quadratic fields of class number one at radii
10^4 to 10^5, (ii) the code that will measure them, and (iii) the criteria that will decide whether the
predictions succeeded. Nothing about holes or primes at radii above 4600 had been measured when these
files were frozen; the only computation run at those radii was a timing test of point generation and
sieving (point counts only, no statistics).

## 1. Background

For a lattice point z of the ring of integers O of Q(sqrt D), D in {-3, -4, -7, -8, -11, -19, -43, -67, -163},
let h(z) be the distance to the nearest prime element and omega(z) the number of prime ideals of norm
<= 31 dividing z. A conditional singular series S_D(u) gives, for each signature D(z), the relative
probability that z + u is prime. A Poisson model then gives, for each signature, the probability that
h(z) >= r as exp(-kappa * Lambda_D(r) / log N(z)), with Lambda_D(r) the sum of S_D(u) over 0 < |u| < r.
The over-representation of each class among the largest holes is predicted from this, and its sign can
change with the scale because Lambda_D(r) - mean changes sign as r crosses successive lattice shells.
See the accompanying preprint and `predict.py`.

## 2. What was fitted on earlier data (radii 500 to 4600, 27 configurations)

* The bias pattern has no free parameter. At the measured record radius, the model puts the observed
  over-representation on the correct side of 1 in 118 of 146 omega classes, and in 53 of the 56 classes
  where it predicts a clear effect (|log ratio| > 0.2).
* The plain model overestimates the record radius r* by one or two lattice shells in every one of the 27
  configurations. One parameter corrects this: r* is predicted with kappa_eff = 1.25 kappa
  (`C_RSTAR`), which gives the measured r* exactly in 20 of 27 configurations and within one shell in 25.
  This parameter is used only to predict r*, never to predict the bias pattern.

## 3. Test configurations (fixed)

For each of the 9 fields, three annuli R1 <= |z| <= R2 with R1 in {10000, 31623, 100000}, of width
chosen to contain about 4 million lattice points of the fundamental sector (27 configurations, listed
below and in `predictions/predictions.json`).

## 4. Measured quantities (fixed; code `verify.py`, statistics as in the preprint)

For each configuration: the non-prime lattice points of the fundamental sector; h(z); r* = 99.9%
quantile of h (numpy linear quantile); for each omega class with at least 500 points and for each prime
ideal of norm <= 31 (conjugates merged), the share among the points with h >= r* divided by the overall
share. Uncertainty by bootstrap over holes (connected components of the points with h >= r* at link
distance 2; 2000 draws, seed 12345).

## 5. Criteria (fixed; code `evaluate.py`)

A (configuration, class) pair is eligible if the expected number of points of the class among the points
with h >= r* is at least 100 and the prediction is clear (|log predicted ratio| > 0.2). Predictions are
taken at the measured r* from the conditional tables of `predictions.json` (the model gives a pattern for
each lattice distance near the predicted r*).

* **P1 (primary): omega classes.** Success if at least 75% of eligible pairs have the observed ratio on
  the same side of 1 as the prediction and the one-sided binomial test against 50% gives p < 0.01.
  Inconclusive if fewer than 10 pairs are eligible (about 24 are expected).
* **P2: record radius.** Success if the measured r* equals the predicted r* in at least half of the 27
  configurations and is within one lattice shell in at least 85%.
* **P3: individual prime ideals.** Same rule as P1, for ideals.
* **P4: predicted changes of sign** between consecutive test radii (unconditional predictions, at the
  predicted r*): reported descriptively, with the number observed.

Any deviation from this protocol will be reported together with the result.

## 6. Predictions

Predicted record radius and its changes between 10^4 and 10^5 (annuli of about 4 million points):

| D | predicted r* at R = 10^4 | predicted changes of r* (R: new r*) |
|---|---|---|
| -3 | 4.359 | 37446: 4.583 |
| -4 | 5.099 | 22688: 5.385, 52041: 5.657, 95864: 5.831 |
| -7 | 7.071 | 10206: 7.280, 13260: 7.483, 26523: 7.616, 42572: 7.937, 50846: 8.000 |
| -8 | 7.141 | 11935: 7.348, 22409: 7.550, 42405: 7.681, 56501: 8.000, 70808: 8.124 |
| -11 | 7.000 | 11099: 7.280, 16763: 7.416, 24494: 7.681, 34977: 7.746, 74598: 8.000, 88315: 8.185 |
| -19 | 7.000 | 13678: 7.416, 36815: 7.810, 61087: 7.937, 98678: 8.000 |
| -43 | 7.211 | 17722: 7.280, 39000: 7.681, 83555: 8.000 |
| -67 | 7.681 | 17457: 8.000, 30755: 8.185, 48558: 8.246 |
| -163 | 7.000 | 12983: 7.280, 60038: 7.810 |

Predicted over-representation among points with h >= r*, at the predicted r* (the full tables, for each
possible r* near the prediction and for every prime ideal, are in `predictions/predictions.json`):

| D | annulus | predicted r* | predicted over-representation among points with h >= r* (omega_31 classes with share >= 2.5%) |
|---|---|---|---|
| -3 | [10000, 10662] | 4.359 | ω=0: 1.18, ω=1: 0.99, ω=2: 0.88, ω=3: 0.94 |
| -3 | [31623, 31832] | 4.359 | ω=0: 1.16, ω=1: 0.99, ω=2: 0.89, ω=3: 0.94 |
| -3 | [100000, 100066] | 4.583 | ω=0: 1.48, ω=1: 1.02, ω=2: 0.76, ω=3: 0.65 |
| -4 | [10000, 10509] | 5.099 | ω=0: 1.62, ω=1: 0.92, ω=2: 0.83, ω=3: 0.92 |
| -4 | [31623, 31784] | 5.385 | ω=0: 1.12, ω=1: 0.83, ω=2: 0.98, ω=3: 1.29 |
| -4 | [100000, 100051] | 5.831 | ω=0: 1.50, ω=1: 0.93, ω=2: 0.85, ω=3: 0.91 |
| -7 | [10000, 10337] | 7.071 | ω=0: 1.55, ω=1: 0.96, ω=2: 0.94, ω=3: 0.94, ω=4: 0.96 |
| -7 | [31623, 31730] | 7.616 | ω=0: 1.36, ω=1: 0.96, ω=2: 0.97, ω=3: 0.95, ω=4: 0.94 |
| -7 | [100000, 100034] | 8.000 | ω=0: 1.82, ω=1: 0.94, ω=2: 0.87, ω=3: 0.91, ω=4: 0.97 |
| -8 | [10000, 10360] | 7.141 | ω=0: 1.30, ω=1: 1.00, ω=2: 0.96, ω=3: 0.92, ω=4: 0.90 |
| -8 | [31623, 31737] | 7.550 | ω=0: 0.83, ω=1: 0.98, ω=2: 0.98, ω=3: 1.12, ω=4: 1.29 |
| -8 | [100000, 100036] | 8.124 | ω=0: 1.35, ω=1: 1.09, ω=2: 0.84, ω=3: 0.89, ω=4: 1.00 |
| -11 | [10000, 10422] | 7.000 | ω=0: 1.00, ω=1: 0.93, ω=2: 0.97, ω=3: 1.12, ω=4: 1.33 |
| -11 | [31623, 31757] | 7.681 | ω=0: 1.30, ω=1: 1.05, ω=2: 0.91, ω=3: 0.87, ω=4: 0.90 |
| -11 | [100000, 100042] | 8.185 | ω=0: 1.09, ω=1: 0.97, ω=2: 0.99, ω=3: 1.01, ω=4: 1.00 |
| -19 | [10000, 10555] | 7.000 | ω=0: 1.18, ω=1: 1.02, ω=2: 0.94, ω=3: 0.92, ω=4: 0.91 |
| -19 | [31623, 31799] | 7.416 | ω=0: 1.21, ω=1: 1.03, ω=2: 0.93, ω=3: 0.89, ω=4: 0.90 |
| -19 | [100000, 100055] | 8.000 | ω=0: 1.24, ω=1: 1.03, ω=2: 0.92, ω=3: 0.88, ω=4: 0.90 |
| -43 | [10000, 10835] | 7.211 | ω=0: 1.18, ω=1: 1.00, ω=2: 0.87, ω=3: 0.76 |
| -43 | [31623, 31887] | 7.280 | ω=0: 1.00, ω=1: 1.02, ω=2: 0.99, ω=3: 0.94 |
| -43 | [100000, 100083] | 8.000 | ω=0: 1.13, ω=1: 1.02, ω=2: 0.88, ω=3: 0.77 |
| -67 | [10000, 11042] | 7.681 | ω=0: 1.19, ω=1: 0.94, ω=2: 0.82, ω=3: 0.76 |
| -67 | [31623, 31953] | 8.185 | ω=0: 1.19, ω=1: 0.94, ω=2: 0.81, ω=3: 0.73 |
| -67 | [100000, 100104] | 8.246 | ω=0: 1.20, ω=1: 0.94, ω=2: 0.79, ω=3: 0.70 |
| -163 | [10000, 11626] | 7.000 | ω=0: 0.91, ω=1: 1.10, ω=2: 1.14 |
| -163 | [31623, 32137] | 7.280 | ω=0: 0.98, ω=1: 1.03, ω=2: 0.98 |
| -163 | [100000, 100163] | 7.810 | ω=0: 1.08, ω=1: 0.92, ω=2: 0.75 |

Predicted changes of sign between consecutive test radii, for classes with enough points to be tested:
D = -8, omega = 0: 0.83 at R = 31623, then 1.35 at R = 100000; D = -8, omega = 3: 1.12, then 0.89;
D = -11, omega = 3: 1.12 at R = 10000, then 0.87 at R = 31623.

## 7. Files and SHA-256 fingerprints

The fingerprints below identify the exact files frozen today.

    a58bd09ac6fbbcced38ac099b4722e1336945904cec9837d296685ade81500cc  holes.py
    70b07f3bad59196e097ac38749b9439374b7f678fd69d35986641c3208a25235  predict.py
    9ef462b3a3669708379c946035533be2c899a373da7882391d79ee1541ce6f7c  make_predictions.py
    1131f761cb4016529cef2a2800aeba5e5627fa71b2417e69cf37d5aabf921658  verify.py
    37625f4ce7e435d05304f28f18165e305c7c92a662a25f200a6583a49ebbc7c4  evaluate.py
    b4404e5b49002d4f335bd30ccde6710da60513809ec43542d2b6820543384603  predictions/predictions.json
