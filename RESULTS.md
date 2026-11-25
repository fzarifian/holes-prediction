# Results of the pre-registered test

Measurement run on 25 November 2026 (00:30 UTC), after the release `v1.0-prediction`, with the code of
commit `102f3877b5308e6c325122bc4c713d21f8362cef` unchanged (`python3 verify.py run`, then
`python3 evaluate.py`). No deviation from the protocol of `PREREGISTRATION.md`.

Files: `verification/measurements.json` (SHA-256 `52c04e228b4573ef9780ece770631c0b710be340f957def7ae582a2037bd9fac`)
and `verification/evaluation.txt` (output of `evaluate.py`).

## Summary

| Criterion | Result | Verdict |
|---|---|---|
| **P1 (primary):** direction of the bias by omega class, at the measured r* | 36 of 36 clear predictions correct (all nine fields contribute); 35 of the 35 with a bootstrap interval excluding 1 | **Success** (p = 1.5 × 10⁻¹¹) |
| **P2:** record radius r* | exact in 5 of 27 configurations, within one lattice shell in 20 of 27 | **Failure** |
| **P3:** direction of the bias by prime ideal, at the measured r* | 33 of 34 clear predictions correct | **Success** (p = 2.0 × 10⁻⁹) |
| **P4:** predicted reversals between test radii, at the predicted r* | 2 of 6 observed | descriptive |

## Reading

**The bias pattern is predicted correctly.** Once the record radius r* is known, the singular-series
model, which has no free parameter for this purpose, gave the direction of the over- or
under-representation of every class tested at radii 10⁴ to 10⁵, in all nine fields.

**The record radius is not.** The one-parameter calibration fitted at radii 500 to 4600 does not carry
over: at radii 10⁴ to 10⁵ the measured r* is never smaller than predicted, and is larger by one shell in
15 configurations and by two shells in 7. The largest holes are therefore closer to the plain Poisson
model at these radii than at the smaller ones.

**The reversals depend on r*.** Since the radius at which r* jumps from one lattice shell to the next was
mispredicted, so were most of the scales where the bias reverses (2 of 6). The four failures all come
from configurations where r* was mispredicted; with the measured r* (criterion P1) the direction is
correct in every tested case.

## Exploratory remarks (not pre-registered)

The direction of the r* misses is systematic. A natural hypothesis is that the deviation of the largest
holes from the Poisson model, large at radii below 5000, decreases as log N grows. Testing this, and
predicting r* itself, would need a model of the correlations between nearby primes, not a constant
correction.
