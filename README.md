# Holes between prime elements: a pre-registered prediction [![DOI](https://zenodo.org/badge/1471016843.svg)](https://doi.org/10.5281/zenodo.24518337)

**Fabien Zarifian** (independent researcher), fabien.zarifian@nuevolia.fr

This repository records, **before any measurement**, a set of numerical predictions about the empty regions
("holes") between prime elements in the nine imaginary quadratic fields of class number one
(Gaussian integers, Eisenstein integers, and the rings of integers of Q(√D) for
D = −7, −8, −11, −19, −43, −67, −163).

The point of publishing it now is simple: anyone can check that the predictions were written down
before the data they will be compared with. The measurement will be run after this release, with the
code included here, and the result will be added later, whether the predictions succeed or fail.

## The question

Take a lattice point z and let h(z) be its distance to the nearest prime element. Does the size of the
hole around z depend on which small primes divide z?

At the radii already studied (500 to 4600), it does, and in a curious way: centres with many small
factors are over-represented among the largest holes at some scales and under-represented at others.
A model built from the Hardy–Littlewood singular series explains when and why the direction of this
bias changes. The model has been fitted on those radii only.

## The prediction

For radii 10⁴ to 10⁵, which have never been measured, the model predicts:

* the typical radius of the largest holes in each field, and the scales where it jumps;
* for each class of centres (by number of small prime factors, and for each small prime ideal),
  whether it is over- or under-represented among the largest holes;
* specific reversals of the bias between scales, for example for D = −8 the centres without small
  factors go from under-represented (×0.83) at radius 31 623 to over-represented (×1.35) at radius 100 000.

The full protocol, the success criteria and the tables are in
[`PREREGISTRATION.md`](PREREGISTRATION.md).

## Contents

| File | Role |
|---|---|
| `PREREGISTRATION.md` | Protocol, success criteria, prediction tables, SHA-256 fingerprints |
| `predictions/predictions.json` | All predictions, in machine-readable form |
| `holes.py` | Library: rings, primes, hole radii, singular series |
| `predict.py`, `make_predictions.py` | The model and the script that produced the predictions |
| `verify.py` | The measurement code that will be run on the new radii (frozen) |
| `evaluate.py` | Applies the success criteria to the measurements (frozen) |
| `results/fields.json` | Earlier measurements (radii 500 to 4600) used to fit the model |
| `verification/` | Measurements at radii 10⁴ to 10⁵ and evaluation output (added after the prediction release) |
| `RESULTS.md` | Outcome of the test |
| `PREREGISTRATION_v2.md` | Second test: protocol, criteria, predictions, fingerprints |
| `predict2.py`, `make_predictions2.py` | Second model (pair correlations) and the script that produced its predictions |
| `verify2.py`, `evaluate2.py` | Measurement and evaluation code of the second test (frozen) |
| `predictions_v2/predictions_v2.json` | All predictions of the second test |
| `verification_v2/`, `RESULTS_v2.md` | Measurements and outcome of the second test |

To check that the files have not changed since this release:

    sha256sum -c <(grep -E '^    [0-9a-f]{64}' PREREGISTRATION.md | sed 's/^    //')

Requirements: Python 3 with numpy and scipy.

## Status

- [x] Predictions frozen (25 November 2026, release `v1.0-prediction`)
- [x] Measurement at radii 10⁴ to 10⁵ (`verification/measurements.json`)
- [x] Evaluation against the pre-registered criteria: see [`RESULTS.md`](RESULTS.md)

**Outcome:** the direction of the bias was predicted correctly in 36 of 36 tested classes (primary
criterion, success) and in 33 of 34 prime ideals; the record radius itself was mispredicted by one or
two lattice shells in most configurations (failure), and with it most of the predicted scales of reversal.

### Second test (radii 2×10⁵ to 10⁶)

The first model treated nearby primes as independent. The second model adds the pair correlations of the
Hardy–Littlewood conjectures and drops the fitted constant: it has no free parameter. On the 54
configurations already measured it gives the record radius exactly in 41 cases (25 for the first model)
and within one lattice shell in 53. Since those data were used to develop it, it is tested again, blind,
at radii never measured.

- [x] Predictions frozen (25 November 2026, release `v2.0-prediction`): see [`PREREGISTRATION_v2.md`](PREREGISTRATION_v2.md)
- [x] Measurement at radii 2×10⁵ to 10⁶ (`verification_v2/measurements_v2.json`)
- [x] Evaluation against the pre-registered criteria: see [`RESULTS_v2.md`](RESULTS_v2.md)

**Outcome:** all four pre-registered criteria are met. Record radius exact in 22 of 27 configurations and
within one lattice shell in all 27; direction of the bias correct in 27 of 28 classes and 43 of 45 prime
ideals; the 6 predicted reversals of the bias between scales were all observed (two of them only
marginally).

## License

All files are released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

## Acknowledgement

The code, the computations and a first draft of the texts were produced with the help of an AI assistant
(Claude, Anthropic). The author is responsible for the content.
