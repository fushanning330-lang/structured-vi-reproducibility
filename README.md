# Structured VI reproducibility archive

Reproduction materials for the manuscript:

**From factor-space instability to auditable structured variational inference: replicate diagnostics, covariance reproducibility, and computational scaling with a pharmacovigilance application**

Manuscript-facing repository revision: **V8.9.1 (28 August 2026)**.

## Files in the repository root

- `structured-vi-reproducibility_REPOSITORY_READY.zip` — complete frozen reproduction bundle corresponding to the V8.9.1 submission materials. It contains derived numerical data, author-generated analysis/simulation code, source-identity manifests, the E1–E4 computational extension, and the publication-output reproduction route.
- `reproduce_plos_outputs.py` — publication-output reproduction entry point. It regenerates Tables 1–5, main Figures 3–5, and Supporting Figures S1–S6 from frozen derived inputs in the extracted archive.
- `reproduce_matched_four_sensitivity.py` — recomputes the post hoc descriptive matched-four sensitivity using the four optimizer seeds common to P1 and P4 (1001, 2002, 3003, 5005). It does not refit any model or rerun the 239 simulation fits.
- `POSTHOC_MATCHED_FOUR_SENSITIVITY_V8_9_1.csv` — frozen matched-four descriptive sensitivity output used in manuscript V8.9.1.
- `S2_CONTENT_SHA256_MANIFEST_V8_9_1.csv` — SHA256 manifest for the contents of the complete V8.9.1 reproduction archive.
- `requirements.txt` — minimal dependencies for the publication-output reproduction route.
- `LICENSE` — MIT License covering author-generated code in this repository. It does not alter rights or terms applying to third-party materials or raw FDA AEMS/FAERS source data.

## Reproduction

1. Download `structured-vi-reproducibility_REPOSITORY_READY.zip` and extract it to a local directory.
2. Follow `README.txt` inside the extracted archive for the main reproduction route.
3. From the extracted archive, run:

```bash
python reproduce_plos_outputs.py
```

The main route performs internal checks and reports `REPRODUCTION CHECK: PASS` when the frozen publication outputs are reproduced successfully.

For the V8.9.1 matched-four sensitivity, run from the repository root (or point `--data-dir` at the extracted archive):

```bash
python reproduce_matched_four_sensitivity.py --data-dir /path/to/extracted/archive
```

The matched-four analysis is **post hoc and descriptive**. Pairwise rows reuse fitted identities and are not treated as independent experimental replicates. No p-values or inferential confidence intervals are generated.

## Scientific scope

The pharmacovigilance data serve as a high-dimensional computational application. No raw FDA regulatory extract is redistributed in this repository. The manuscript and repository do not claim universal variational-inference convergence, universal particle-number superiority, universal end-to-end runtime superiority, or clinical signal-detection validity.

The V3 architecture is a post-diagnostic adaptive redesign. Representation reproducibility, covariance reproducibility, and known-target accuracy are treated as distinct evidence layers.

## Data source

Raw FDA AEMS/FAERS quarterly public-use files are available from the U.S. Food and Drug Administration. This repository distributes author-generated code and derived/frozen reproduction materials rather than raw regulatory extracts.

## License

Author-generated code is released under the MIT License. See `LICENSE`.
