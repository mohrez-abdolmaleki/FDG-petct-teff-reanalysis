# FDG PET/CT External Dose Rate and Effective Half-Life — Reanalysis

Reanalysis of external dose-rate measurements from 36 patients undergoing
¹⁸F-FDG PET/CT, estimating the effective half-life (Teff) of the dose rate
using a linear mixed-effects model that accounts for the repeated-measures
structure of the data (3 measurements per patient).

## Background

An earlier analysis of this dataset reported Teff ≈ 43.5 min. That value
could not be reproduced from the underlying patient-level data by any
aggregation method tried (pooled regression, time-binned averaging,
positional averaging matching the original spreadsheet's own summary
rows) — all methods converge on Teff ≈ 85–95 min instead. This repository
contains the corrected, from-scratch reanalysis.

## Repository structure

```text
├── CSV.csv                                    # raw patient data (as provided)
├── lmm.py                                     # linear mixed model engine
│                                               # (REML/ML, hand-implemented —
│                                               # see note below). This is a
│                                               # MODULE, imported by the
│                                               # notebooks below -- it stays
│                                               # .py, it is not run on its own.
├── 00_validate_lmm.ipynb                      # validates lmm.py against
│                                               # simulated data with known
│                                               # ground-truth parameters --
│                                               # run this first
├── 01_data_cleaning_and_teff_reconciliation.ipynb
│                                               # cleans raw data, resolves the
│                                               # duplicate-row issue (currently
│                                               # patient P08), cross-checks
│                                               # pooled Teff across several
│                                               # naive methods (diagnostic
│                                               # only -- see Methodology)
├── 02_final_model_and_results.ipynb           # descriptive stats, model
│                                               # selection, final Teff/Tbio
│                                               # with bootstrap CI, covariate
│                                               # analysis, figures
├── 03_sex_comparison.ipynb                    # tests sex as a covariate on
│                                               # intercept and slope; null
│                                               # result (see Key results)
├── 04_statsmodels_reimplementation.ipynb      # independent cross-check of
│                                               # the main mixed-model results
│                                               # using statsmodels.MixedLM;
│                                               # not the primary analysis
└── patient_data_clean.csv                     # output of notebook 01
```

Note on `.py` vs `.ipynb`: `lmm.py` is a plain module (functions only, no
top-level analysis) and must stay `.py` so it can be imported. Jupyter
notebooks can't be imported as modules the same way. The analysis notebooks
keep each step's output visible inline, which is generally the more
readable/reviewable format for a GitHub research repository.

Notebook `04_statsmodels_reimplementation.ipynb` is intentionally separate
from the primary analysis. It independently fits the same linear mixed-effects
models using `statsmodels.MixedLM` and is included as a reproducibility and
implementation cross-check rather than as a replacement for `lmm.py`.

## Requirements

Python 3.9+, with:

```text
numpy
pandas
scipy
matplotlib
jupyter
```

The primary analysis does not require `statsmodels`: `lmm.py` is a
self-contained mixed-model implementation (see "Why a hand-written mixed
model" below).

The independent cross-check notebook additionally requires:

```text
statsmodels
```

The `statsmodels` implementation uses `statsmodels.formula.api.mixedlm`
with the BFGS optimizer.

## How to reproduce

Keep `lmm.py` in the same folder as the notebooks, then run in order,
each with **Restart Kernel and Run All** (don't run cells out of order —
later cells depend on earlier ones):

1. `00_validate_lmm.ipynb` — confirms the model implementation is correct
   before trusting it on real data
2. `01_data_cleaning_and_teff_reconciliation.ipynb` — produces
   `patient_data_clean.csv` (requires the raw CSV to include a `sex`
   column; the cleaning step carries it through automatically if present)
3. `02_final_model_and_results.ipynb` — final results and figures
4. `03_sex_comparison.ipynb` — tests sex as a covariate (optional, only
   needed if you want to reproduce the sex-comparison result)
5. `04_statsmodels_reimplementation.ipynb` — optional independent
   cross-check using `statsmodels.MixedLM`

The `statsmodels` notebook should be treated as a validation/cross-check
step, not as a prerequisite for reproducing the primary analysis.

If you'd rather run the notebooks as plain scripts (e.g. in CI), each
notebook can be executed non-interactively with:

```bash
jupyter nbconvert --to notebook --execute 00_validate_lmm.ipynb
```

## Methodology summary

* **Outcome:** log(dose rate), modeled as a linear function of time since
  injection.
* **Data structure:** 36 patients × 3 repeated measurements each (one
  patient has 2 valid measurements after removing a duplicate data-entry
  row — see below). This repeated-measures structure means observations
  within a patient are correlated, so it's modeled as a linear
  mixed-effects model rather than pooling all points as independent or
  fitting each patient separately.
* **Random-effects structure:** compared a random-intercept-only model
  against a random-intercept-plus-slope model. The random-intercept-only
  model was selected based on AIC, BIC, and a likelihood-ratio test — the
  slope-variance model's fitted random effects were degenerate
  (intercept–slope correlation ≈ 1.0), an expected consequence of having
  only 2–3 observations per patient (not enough within-patient
  information to separate an individual decay rate from an individual
  starting level).
* **Uncertainty:** non-parametric bootstrap over patients (not
  observations), refitting the model on each resample, used for the
  final Teff/Tbio confidence intervals rather than relying on asymptotic
  normal-theory intervals given the modest sample size (n=36).
* **Bootstrap grouping:** each bootstrap draw is assigned a unique
  patient-group identifier. This is necessary because the same original
  patient can be sampled more than once in a bootstrap replicate; retaining
  the original patient ID would incorrectly merge repeated bootstrap draws
  into a single mixed-model group.
* **Covariates:** weight, BMI, and injected activity were tested as
  predictors of both the intercept (starting dose rate) and the slope
  (decay rate). Weight and injected activity are correlated (r≈0.76,
  consistent with roughly weight-based dosing); fitting both together
  shows injected activity, not weight, is the variable associated with
  decay rate in this dataset.
* **Independent implementation check:** the main model was independently
  reimplemented with `statsmodels.MixedLM`. The cross-check reproduces the
  main Teff result and the non-significant random-intercept/random-slope
  comparison, while the joint covariate analysis gives very similar
  estimates and p-values. Small numerical differences are expected from
  differences between independent implementations and optimizers.

### Bootstrap implementation note

The bootstrap resamples patients with replacement. Within each replicate,
each sampled patient occurrence receives a unique group ID of the form
`boot_<replicate>_<draw>`. This preserves the intended bootstrap structure
when a patient is selected multiple times.

The bootstrap fits use BFGS in the `statsmodels` cross-check. Convergence
warnings are monitored separately from actual fit failure: a fit is not
classified as failed solely because `statsmodels` emits a
`ConvergenceWarning`; the fitted model must also be checked for
`fit.converged` and valid finite fixed-effect estimates.

## Key results

| Quantity        | Estimate  | 95% CI          |
| --------------- | --------- | --------------- |
| Population Teff | 87.8 min  | [74.7, 105.0]   |
| Population Tbio | 438.6 min | [233.9, 2292.1] |

No statistically significant difference in decay rate or starting dose
rate was found between male and female patients (LRT p=0.77 on the
decay rate, p=0.47 on the intercept, p=0.65 jointly; robust to
adjustment for injected activity, p=0.88). See `03_sex_comparison.ipynb`.

The independent `statsmodels.MixedLM` cross-check reproduced a Teff close
to 87.8 min and a non-significant random-intercept/random-slope comparison.
In the joint weight-plus-activity model, the estimated p-value for the
weight-by-time term was approximately 0.72 and for the activity-by-time
term approximately 0.02. The latter is close to, but not numerically
identical to, the corresponding result from the hand-written
implementation (approximately 0.026).

Note: the Tbio confidence interval is very wide by construction — Tbio is
derived as 1/(1/Teff − 1/T_phys), a difference of two similar-sized rate
constants, which is numerically unstable when Teff is not much shorter
than the physical half-life (109.77 min for F-18). Teff should be treated
as the reliable headline quantity; Tbio should be reported with this
instability made explicit, or omitted.

## Data quality notes

* One duplicate row (patient P08: time=22, dose rate=29.2 entered twice)
  was identified and removed rather than imputed. The mixed-effects
  model naturally accommodates patients with fewer observations, so no
  substitute value was fabricated for the missing third measurement.
  The cleaning code identifies this patient dynamically from the
  duplicate itself rather than hardcoding an ID, since the raw file's
  patient-ID scheme has already changed once (plain integers -> "P01".."P36").

## Why a hand-written mixed model

`lmm.py` implements REML/ML estimation for a linear mixed model directly
via closed-form per-cluster marginal likelihoods (exact, not an
approximation, since clusters here are only 2–3 observations). It was
written because `statsmodels`/`lme4` were not available in the
development environment at the time. It is validated in
`00_validate_lmm.ipynb` against synthetic data with known parameters
before being trusted on the real data.

The later `04_statsmodels_reimplementation.ipynb` provides an independent
implementation of the same model using `statsmodels.MixedLM`. Its purpose
is to verify that the principal estimates are not artifacts of the
hand-written implementation.

Absolute likelihood-based information criteria should not be compared
blindly between the two implementations because the hand-written
likelihood uses a different constant convention. Model comparisons within
the same likelihood framework remain the relevant comparison.

## What belongs in the manuscript vs. this repository

This repository intentionally contains more than what should appear in
the paper — some of it exists only to show *how* the corrected numbers
were reached and verified, not because it should be reported as a result.

**Goes in the manuscript (Methods/Results):**

* Data cleaning approach (duplicate-row handling), described in prose
* The model comparison table (AIC/BIC/LRT) justifying the
  random-intercept-only model
* Final Teff (and, with the numerical-instability caveat, Tbio) with
  bootstrap 95% CIs
* The covariate finding (injected activity, not weight, associated with
  decay rate), reported as hypothesis-generating
* The sex-comparison result (no significant difference), reported briefly
* Descriptive statistics table
* Figure: population dose-rate curve with CI band
* Figure: raw per-patient trajectories vs. shrinkage-fitted trajectories
  (visual justification for the model choice)

**Stays in the repository only (not in the manuscript text):**

* `00_validate_lmm.ipynb` — cite it in one sentence in Methods
  ("implementation validated on simulated data; see repository"), don't
  reproduce it in the paper
* The naive pooled/binned/positional reconciliation methods in notebook
  01 — these were diagnostic, used only to confirm the original 43.5 min
  figure doesn't reproduce from the data. At most, one sentence in a
  robustness footnote ("a simple pooled fit ignoring clustering gave a
  consistent ~89 min") is defensible; the code itself doesn't belong in
  the manuscript
* Per-patient BLUP table — supplementary material at most
* Bootstrap-distribution histogram figure — supplementary material at
  most, not a main-text figure
* `04_statsmodels_reimplementation.ipynb` — independent implementation
  check and reproducibility evidence; not part of the primary analysis

**Doesn't belong anywhere (already dropped):**

* The earlier R script (`teff_reconciliation.R`) — fully superseded by
  this Python pipeline

## Limitations

* n = 36 patients; multiple comparisons were run in the covariate
  analysis (report as hypothesis-generating, not confirmatory)
* Age was not collected in this dataset
* Only 2–3 measurements per patient — cannot estimate individual
  per-patient decay rates, only a population-average rate with
  patient-specific starting levels
* The bootstrap confidence intervals are based on 36 patient-level
  sampling units; bootstrap estimates should therefore be interpreted
  in the context of the modest sample size.
* The independent `statsmodels` implementation is a cross-check rather
  than a second source of independent data; agreement between the two
  implementations supports the computational result but does not remove
  the underlying limitations of the dataset.

## License

TBD

## Citation

TBD once the associated manuscript is finalized.

