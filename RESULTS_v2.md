# Results of the second pre-registered test

Measurement run on 25 November 2026 (01:20 UTC), after the release `v2.0-prediction`, with the code of
that release unchanged (`python3 verify2.py run`, then `python3 evaluate2.py`). No deviation from the
protocol of `PREREGISTRATION_v2.md`.

Files: `verification_v2/measurements_v2.json` (SHA-256 `2c850497b09f15d0134433915221248a71c40587dbea3d18b57eb0320aeda47c`)
and `verification_v2/evaluation_v2.txt` (output of `evaluate2.py`).

## Summary

| Criterion | Result | Verdict |
|---|---|---|
| **P1 (co-primary):** direction of the bias by omega class, at the measured r* | 27 of 28 clear predictions correct; 27 of the 27 with a bootstrap interval excluding 1 | **Success** (p = 1.1 × 10⁻⁷) |
| **P2 (co-primary):** record radius r* | exact in 22 of 27 configurations, within one lattice shell in 27 of 27 | **Success** |
| **P3:** direction of the bias by prime ideal, at the measured r* | 43 of 45 correct | **Success** (p = 2.9 × 10⁻¹¹) |
| **P4:** predicted reversals between test radii, at the predicted r* | 6 of 6 observed | **Success** |

All four pre-registered criteria are met, at radii up to 10⁶ never measured before the predictions were
published, by a model with no fitted parameter.

## Details and limits

* **Record radius.** The five misses are each by one lattice shell: for D = −7 at all three radii and for
  D = −19 at R = 4.5×10⁵ the measured r* is one shell smaller than predicted; for D = −11 at R = 10⁶ it is
  one shell larger. The misses of the first model all went in the same direction; these do not.
* **Reversals.** The six predicted reversals are all observed, including the double reversal for D = −4
  (centres without small factor: ×0.71, then ×1.33, then ×0.68 for predicted ×0.75, ×1.36, ×0.72). In two
  of them one of the observed values is close to 1, with a bootstrap interval that contains 1
  (D = −4, ω = 3 at R = 4.5×10⁵: 0.97, interval 0.89–1.04; D = −8, ω = 0 at R = 4.5×10⁵: 0.95,
  interval 0.87–1.03). The direction is right in both, but these two are not individually significant.
* **The one miss in P1** is D = −8, ω = 4 at R = 10⁶: predicted ×1.26, observed ×0.95 (interval 0.82–1.08).

## Comparison with the first test

| | First model (independent primes + 1 fitted constant), radii 10⁴–10⁵ | Second model (pair correlations, no fitted constant), radii 2×10⁵–10⁶ |
|---|---|---|
| Bias direction at measured r* | 36 / 36 | 27 / 28 |
| Record radius exact / within one shell | 5 / 27, 20 / 27 | 22 / 27, 27 / 27 |
| Predicted reversals observed | 2 / 6 | 6 / 6 |
