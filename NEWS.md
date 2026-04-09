# NEWS: SmoothPLS

## SmoothPLS 0.1.3 (2026-04-09)

### Improvements 
* **Documentation:**
  * Launched the official `pkgdown` website (hosted on GitHub Pages) including detailed vignettes and function references.
  * Added the complete compiled PDF manual to `inst/doc/`.
* **README & Branding:**
  * Added a comprehensive quick-start example (One-State Categorical PLS) with visual outputs.
  * Fixed LaTeX equations rendering for cross-compatibility between GitHub and Pandoc/pkgdown.
  * Added institutional links (Decathlon SportsLab, Inria).
* **Continuous Integration:**
  * Configured GitHub Actions workflows for automated `R CMD check` and `pkgdown` site deployment.

## SmoothPLS 0.1.2 (2026-04-08)

### Core Improvements & Stability
* **Numerical Precision**: Optimized categorical integration in `evaluate_id_func_integral` with stricter relative tolerance (`rel.tol`) and increased subdivisions (1000) for high-order B-splines.
* **Analytic Prediction**: Implemented analytic L2 inner product for Scalar Functional Data (SFD) using `fda::inprod`, replacing discrete trapezoidal integration for near-perfect precision.
* **Safety Checks**: Added time-range assertions in `smoothPLS_predict` to prevent silent errors when predicting on data outside the basis domain.

### Bug Fixes & Refactoring
* **Tidyselect Compatibility**: Fixed deprecation warnings by implementing `all_of()` in data pivoting functions.
* **Multivariate Support**: Corrected logical assertions in `smoothPLS` to properly handle mixed lists of categorical and numerical predictors.
* **Integration Robustness**: Added `stop.on.error = FALSE` in segment integration to handle micro-intervals without crashing the full model.

### Testing
* **Core Test Suite**: Added 70 unit tests covering Theorems (univariate equivalence), score orthogonality, and prediction consistency.
* **Edge Cases**: Added tests for time-mismatch handling and multi-state categorical transitions.

---

## SmoothPLS 0.1.1 (2026-03-20)

### Improvements
* **Code Refactoring**: Modularization of internal functions for Lambda matrix evaluation.
* **Synthetic Data**: Improved `generate_X_df` and `generate_Y_df` for more realistic categorical state transitions.
* **S3 Structure Prep**: Initial work on internal objects to support future S3 methods (print, plot, predict).

---

## SmoothPLS 0.1.0 (2025-12-15)

### Initial Release
* **Thesis Milestone**: First functional version used for the initial examples in the doctoral thesis.
* **Core Algorithms**: Implementation of Smooth PLS for Hybrid Functional Data (CFD and SFD).
* **Basis Expansion**: Support for B-spline basis representation of functional predictors.
* **Categorical Handling**: Implementation of the "active area" integration concept for state-based predictors.
