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

To check that the files have not changed since this release:

    sha256sum -c <(grep -E '^    [0-9a-f]{64}' PREREGISTRATION.md | sed 's/^    //')

Requirements: Python 3 with numpy and scipy.

## Status

- [x] Predictions frozen (25 November 2026)
- [ ] Measurement at radii 10⁴ to 10⁵
- [ ] Evaluation against the pre-registered criteria

## License

All files are released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

## Acknowledgement

The code, the computations and a first draft of the texts were produced with the help of an AI assistant
(Claude, Anthropic). The author is responsible for the content.
