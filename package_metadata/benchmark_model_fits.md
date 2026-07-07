# EDI Exhaustive C++ Model Fit Benchmarks

_Generated: 2026-07-07 21:56:40 JST_

This report compares the performance of EDI's Rcpp-optimized model fitting paths against **low-level** canonical R implementations (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) where possible.

## Compilation Context

These rows are read from build metadata compiled into the loaded `EDI` shared object via `edi_build_info_cpp()`.

**Compilation warning:** EDI model-fit timings are sensitive to the compiler flags used to build the loaded `EDI.so`. If EDI is compiled without the proper optimized flags, or with flags that are known to degrade these kernels such as problematic LTO builds, the benchmark can show substantial performance regressions that reflect the binary build rather than the modeling algorithms.

*   **EDI shared object:** `/home/kapelner/R/x86_64-pc-linux-gnu-library/4.7/EDI/libs/EDI.so`
*   **EDI shared object mtime:** `2026-07-07 20:34:41`
*   **Capture method:** `configure-generated header compiled into EDI.so`
*   **Build timestamp:** `2026-07-07 08:06:50 JST`
*   **Build host:** `LAPTOP-J2T9TGGB`
*   **R version at build:** `R Under development (unstable) (2026-04-23 r89955) -- "Unsuffered Consequences"`
*   **R `CXX20` at build:** `g++`
*   **R `CXX20STD` at build:** `-std=gnu++20`
*   **R `CXX20FLAGS` at build:** `-O3 -march=native -funroll-loops -fno-math-errno`
*   **R `SHLIB_OPENMP_CXXFLAGS` at build:** `unavailable`
*   **Build env at build:** `EDI_PORTABLE=0`, `EDI_DISABLE_VECTORIZATION=0`, `EDI_NATIVE_SPEED=1`, `EDI_NATIVE_LTO=0`
*   **Package `PKG_CPPFLAGS` at build:** `-I../inst/include`
*   **Package `PKG_CXXFLAGS` at build:** `$(SHLIB_OPENMP_CXXFLAGS) -DNDEBUG -DEIGEN_NO_DEBUG -Wno-ignored-attributes -march=native -mtune=native -fno-lto; override CXXFLAGS+=-O3`
*   **Package `PKG_LIBS` at build:** `$(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -ltbb12 -fstack-protector`
*   **Compiler reported by binary:** `15.2.0`
*   **Compiler optimization macro enabled:** `TRUE`
*   **Compiler fast-math macro enabled:** `FALSE`
*   **Eigen vectorization disabled macro enabled:** `FALSE`

## Benchmark Dataset Specification

All benchmarks were performed on a synthetic clinical-trial-scale dataset generated for each response type. The data generation process ensures numerical stability and fair solver comparison by using the following parameters:

*   **Sample Size ($N$):** 1,000 subjects for most models; 500 subjects for survival models. Exact and trend tests may use smaller scaled samples (N=100-500) as noted in the results.
*   **Predictors ($p$):** 5 total predictors, including a global intercept, a balanced binary treatment assignment from fixed `iBCRD`, and 4 continuous covariates ($X \sim \text{Normal}(0, 1)$).
*   **Effect Sizes:** Covariate coefficients are sampled from $\text{Normal}(0, 0.5)$. The treatment coefficient is set to 0.5 in the linear predictor so the benchmarked treatment effect is meaningfully separated from zero.
*   **EDI Design Template:** EDI benchmark objects are instantiated on a fixed `iBCRD` design.
*   **Response Generation:**
    *   **Continuous:** Linear model with additive $\text{Normal}(0, 0.5)$ noise.
    *   **Incidence:** Binary outcomes via a Logistic link.
    *   **Count:** Integer outcomes via Poisson or Negative Binomial distributions with an exponential link.
    *   **Proportion:** Continuous outcomes in $(0, 1)$ via a Beta distribution with a logit link.
    *   **Survival:** Exponentially distributed event times with approximately 20% random censoring.
    *   **Ordinal:** 3-level categorical outcomes generated from the same ordinal construction used elsewhere in the benchmark suite.
*   **Stratified Cox Exception:** For `InferenceSurvivalStratCoxPHRegr`, the benchmark injects low-cardinality covariates before outcome generation so the row exercises a genuinely stratified Cox fit rather than the unstratified fallback.

## Methodology

*   **Bare Metal EDI Timing:** EDI rows call the exported C++ functions directly (e.g., `fast_logistic_regression_cpp`, `fast_ordinal_regression_cpp`) with all design matrices and fixed inputs pre-built outside the timed region. There is no R6 object instantiation, no cached state management, no warm start storage, and no standard error computation in the timed region — only the raw numerical solver.
*   **Apples-to-Apples Canonical Timing:** Canonical R timings likewise call the lowest-level publicly exposed interfaces (e.g., `glm.fit`, `lm.fit`, `coxph.fit`) with pre-built design matrices. If a canonical package exposes no low-level function, the formula-based API is used instead.
*   **Low-Level Comparison:** Both EDI and canonical timings are measured on pre-built numeric matrices, removing formula parsing, model-frame construction, and R6/S3/S4 dispatch overhead from the timed region wherever the API permits.
*   **Limitation:** Some canonical comparators only expose formula-based APIs. Those rows remain included but their canonical timings carry formula/model-frame overhead not present in the EDI bare-metal timing.
*   **Averaging:** All timings are medians over 30 cold estimate-only timing samples measured with adaptive batched `system.time`; paths below 0.01 ms use `microbenchmark(times = 5000)` instead.
*   **Timing P-Value:** `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.
*   **Row Highlighting:** Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons.
*   **Constraints**: Matched-pair/KK and highly custom paths are excluded as per user request.

## Results

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.16</td><td>stats</td><td>HL median pairwise diff</td><td>2.28</td><td>14.25x</td><td>5.42e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.04</td><td>stats</td><td>lm.fit</td><td>0.10</td><td>2.78x</td><td>3.99e-05</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>1.75</td><td>quantreg</td><td>rq.fit</td><td>1.61</td><td>0.92x</td><td>0.00157</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.10</td><td>MASS</td><td>rlm(MM)</td><td>45.50</td><td>457.25x</td><td>1.67e-44</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidBinomialIdentityRiskDiff</td><td>incidence</td><td>0.08</td><td>stats</td><td>glm.fit(ident)</td><td>12.08</td><td>150.01x</td><td>1.46e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.21</td><td>stats</td><td>glm.fit+gcomp(RD)</td><td>1.88</td><td>8.87x</td><td>9.22e-20</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.20</td><td>stats</td><td>glm.fit+gcomp(RR)</td><td>1.90</td><td>9.54x</td><td>2.15e-20</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>0.98</td><td>stats</td><td>glm.fit(log)</td><td>4.00</td><td>4.1x</td><td>8.72e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.17</td><td>stats</td><td>glm.fit</td><td>1.91</td><td>11.54x</td><td>1.98e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidModifiedPoisson</td><td>incidence</td><td>0.21</td><td>stats</td><td>glm.fit(modified)</td><td>2.68</td><td>12.97x</td><td>2.35e-21</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.27</td><td>stats</td><td>glm.fit(probit)</td><td>1.81</td><td>6.65x</td><td>2.42e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.02</td><td>stats</td><td>lm.fit(LPM)</td><td>0.13</td><td>5.89x</td><td>1.59e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>1.65</td><td>pscl</td><td>hurdle(nb)</td><td>47.62</td><td>28.81x</td><td>6.61e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>1.51</td><td>pscl</td><td>hurdle</td><td>17.82</td><td>11.82x</td><td>1.07e-39</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>0.47</td><td>MASS</td><td>glm.nb</td><td>57.63</td><td>123.03x</td><td>7.41e-33</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.15</td><td>stats</td><td>glm.fit</td><td>1.83</td><td>11.98x</td><td>4.69e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.14</td><td>stats</td><td>glm.fit(quasi)</td><td>1.63</td><td>11.5x</td><td>1.77e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.21</td><td>stats</td><td>glm.fit</td><td>1.87</td><td>8.93x</td><td>2.69e-19</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>1.03</td><td>pscl</td><td>zeroinfl(nb)</td><td>144.75</td><td>141.18x</td><td>1.71e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>3.54</td><td>pscl</td><td>zeroinfl</td><td>69.87</td><td>19.74x</td><td>7.19e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>1.02</td><td>betareg</td><td>betareg.fit</td><td>31.00</td><td>30.43x</td><td>8.52e-38</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropFractionalLogit</td><td>proportion</td><td>0.14</td><td>stats</td><td>glm.fit(quasi)</td><td>1.27</td><td>9.27x</td><td>1.06e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.18</td><td>stats</td><td>glm.fit(quasi)+gcomp</td><td>1.86</td><td>10.33x</td><td>6.25e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.34</td><td>survival</td><td>coxph.fit(breslow)</td><td>0.63</td><td>1.85x</td><td>4.47e-08</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>0.02</td><td>survival</td><td>survfit(median)</td><td>4.20</td><td>225.69x</td><td>1.53e-36</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.02</td><td>survival</td><td>survdiff</td><td>2.31</td><td>150.29x</td><td>7.22e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalRestrictedMeanDiff</td><td>survival</td><td>0.02</td><td>survival</td><td>survfit(rmean)</td><td>3.23</td><td>188.22x</td><td>5.03e-30</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.84</td><td>survival</td><td>coxph.fit(strat)</td><td>1.10</td><td>1.31x</td><td>6.59e-17</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.06</td><td>survival</td><td>survreg</td><td>3.64</td><td>66.11x</td><td>7e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.46</td><td>VGAM</td><td>vglm(acat)</td><td>14.00</td><td>30.25x</td><td>1.04e-34</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCauchitRegr</td><td>ordinal</td><td>0.84</td><td>ordinal</td><td>clm(cauchit)</td><td>17.65</td><td>21x</td><td>2.63e-12</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalCloglogRegr</td><td>ordinal</td><td>0.44</td><td>ordinal</td><td>clm(cloglog)</td><td>8.50</td><td>19.13x</td><td>1.11e-24</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.22</td><td>VGAM</td><td>vglm(cratio)</td><td>14.03</td><td>65.09x</td><td>2.7e-21</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.66</td><td>ordinal</td><td>clm+gcomp</td><td>16.07</td><td>24.51x</td><td>2.94e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalOrderedProbitRegr</td><td>ordinal</td><td>0.46</td><td>ordinal</td><td>clm(probit)</td><td>7.47</td><td>16.33x</td><td>4.28e-23</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.48</td><td>ordinal</td><td>clm</td><td>7.87</td><td>16.38x</td><td>2.83e-20</td><td>***</td></tr>
  </tbody>
</table>

## Wald Test Performance (Full Inference)

This table compares the performance of **Full Inference** (Model Fit + Standard Error calculation + P-value derivation).
Unlike the point-estimation table above, these results include the computational cost of the variance-covariance matrix (Hessian or Fisher Information) and the Wald test statistic calculation.
All paths (EDI and Canonical) use a reduced sample size ($N=200$) for this full-inference benchmark to ensure iterative stability.
**Stratified Cox Exception**: For `InferenceSurvivalStratCoxPHRegr`, the benchmark injects low-cardinality covariates before outcome generation so the row exercises a genuinely stratified Cox fit rather than the unstratified fallback.
EDI regression models (Logistic, Poisson) are benchmarked using the **IRLS** optimizer for these Wald tests.
**Solver-Only Prebuilds**: Benchmark setup prebuilds exposed observed-data design matrices, reduced design matrices, strata IDs, and other fixed working inputs outside the timed region when the implementation exposes those hooks. The timed region then measures the full-inference kernel on those fixed inputs.
**Limitation**: Some canonical comparators only expose formula-based APIs rather than comparable low-level fit kernels. Those rows remain included, but their canonical timings may still contain formula/model-frame overhead beyond the numerical solver, variance, and p-value work itself.
**Timing Note**: All timings are medians over 30 warmed runs measured with adaptive batched `system.time`; paths below 0.01 ms use `microbenchmark(times = 5000)` instead.
**Timing P-Value**: `Timing Pval` reports a Welch two-sample t-test comparing the EDI and canonical timing replicate distributions for each row. The unlabeled final column marks thresholds with `***` for p < 0.001, `**` for p < 0.01, and `*` for p < 0.05.
**Row Highlighting**: Light green rows indicate `Speedup > 1` and `Timing Pval < 0.05`; light grey rows indicate `NA` timing comparisons.

<table>
  <thead>
    <tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th><th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th><th>Speedup</th><th>Timing Pval</th><th></th></tr>
  </thead>
  <tbody>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleMeanDiffPooledVar</td><td>continuous</td><td>0.04</td><td>stats</td><td>t.test(pool)</td><td>0.15</td><td>4.05x</td><td>4.8e-35</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceAllSimpleWilcox</td><td>continuous</td><td>0.04</td><td>stats</td><td>wilcox.test</td><td>0.57</td><td>13.33x</td><td>5.07e-32</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinLin</td><td>continuous</td><td>0.22</td><td>stats</td><td>lm.fit(interact)+Wald</td><td>0.47</td><td>2.11x</td><td>6.24e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinOLS</td><td>continuous</td><td>0.01</td><td>stats</td><td>lm.fit+Wald</td><td>0.09</td><td>6.44x</td><td>6e-15</td><td>***</td></tr>
    <tr><td>InferenceContinQuantileRegr</td><td>continuous</td><td>3.26</td><td>quantreg</td><td>rq+summary</td><td>2.21</td><td>0.68x</td><td>2.4e-05</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceContinRobustRegr</td><td>continuous</td><td>0.05</td><td>MASS</td><td>rlm+summary</td><td>1.22</td><td>26.46x</td><td>3.74e-48</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidExactFisher</td><td>incidence</td><td>0.69</td><td>stats</td><td>fisher.test</td><td>0.81</td><td>1.18x</td><td>1.27e-22</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskDiff</td><td>incidence</td><td>0.07</td><td>stats</td><td>glm+gcomp(RD)+Wald</td><td>2.05</td><td>29x</td><td>2.54e-44</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidGCompRiskRatio</td><td>incidence</td><td>0.07</td><td>stats</td><td>glm+gcomp(RR)+Wald</td><td>1.90</td><td>27.52x</td><td>5.62e-44</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogBinomial</td><td>incidence</td><td>0.82</td><td>stats</td><td>glm.fit+Wald(log)</td><td>3.27</td><td>3.99x</td><td>8.93e-16</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidLogRegr</td><td>incidence</td><td>0.04</td><td>stats</td><td>glm.fit+Wald</td><td>0.65</td><td>14.9x</td><td>2.06e-26</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidMiettinenNurminenRiskDiff</td><td>incidence</td><td>0.01</td><td>DescTools</td><td>BinomDiffCI(mn)</td><td>0.59</td><td>53.83x</td><td>3.15e-42</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidNewcombeRiskDiff</td><td>incidence</td><td>0.08</td><td>DescTools</td><td>BinomDiffCI(score)</td><td>0.69</td><td>8.52x</td><td>3.56e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidProbitRegr</td><td>incidence</td><td>0.08</td><td>stats</td><td>glm.fit(probit)+Wald</td><td>0.85</td><td>10.97x</td><td>1.45e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceIncidRiskDiff</td><td>incidence</td><td>0.01</td><td>stats</td><td>prop.test</td><td>0.38</td><td>26.85x</td><td>1.94e-43</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdleNegBin</td><td>count</td><td>0.26</td><td>pscl</td><td>hurdle(nb)+summary</td><td>10.87</td><td>42.21x</td><td>1.4e-48</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountHurdlePoisson</td><td>count</td><td>0.21</td><td>pscl</td><td>hurdle+summary</td><td>8.42</td><td>39.86x</td><td>2.39e-48</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountNegBin</td><td>count</td><td>0.12</td><td>MASS</td><td>glm.nb+summary</td><td>12.31</td><td>101.55x</td><td>1.2e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountPoisson</td><td>count</td><td>0.05</td><td>stats</td><td>glm.fit+Wald</td><td>0.81</td><td>15.46x</td><td>1.17e-39</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountQuasiPoisson</td><td>count</td><td>0.05</td><td>stats</td><td>glm.fit+Wald(quasi)</td><td>0.89</td><td>16.96x</td><td>1.44e-27</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountRobustPoisson</td><td>count</td><td>0.10</td><td>sandwich</td><td>glm+vcovHC</td><td>3.57</td><td>37.41x</td><td>1.92e-48</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedNegBin</td><td>count</td><td>1.18</td><td>pscl</td><td>zeroinfl(nb)+summary</td><td>236.00</td><td>199.17x</td><td>1.83e-37</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceCountZeroInflatedPoisson</td><td>count</td><td>1.59</td><td>pscl</td><td>zeroinfl+summary</td><td>36.08</td><td>22.7x</td><td>5.8e-25</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropBetaRegr</td><td>proportion</td><td>0.33</td><td>betareg</td><td>betareg+summary</td><td>13.00</td><td>39.9x</td><td>5.25e-50</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferencePropGCompMeanDiff</td><td>proportion</td><td>0.06</td><td>stats</td><td>glm(quasi)+gcomp+Wald</td><td>1.80</td><td>30.56x</td><td>1.68e-48</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalCoxPHRegr</td><td>survival</td><td>0.11</td><td>survival</td><td>coxph.fit(breslow)+Wald</td><td>0.39</td><td>3.42x</td><td>2.78e-44</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalGehanWilcox</td><td>survival</td><td>1.40</td><td>survival</td><td>survdiff(rho=1)</td><td>1.47</td><td>1.05x</td><td>0.00132</td><td>**</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalKMDiff</td><td>survival</td><td>2.68</td><td>survival</td><td>survfit(median)+CI</td><td>2.93</td><td>1.1x</td><td>0.0454</td><td>*</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalLogRank</td><td>survival</td><td>0.01</td><td>survival</td><td>survdiff</td><td>1.44</td><td>135.06x</td><td>6.49e-44</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalStratCoxPHRegr</td><td>survival</td><td>0.42</td><td>survival</td><td>coxph.fit(strat)+Wald</td><td>0.55</td><td>1.32x</td><td>1.48e-40</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceSurvivalWeibullRegr</td><td>survival</td><td>0.06</td><td>survival</td><td>survreg+summary</td><td>2.91</td><td>46.57x</td><td>2.81e-45</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalAdjCatLogitRegr</td><td>ordinal</td><td>0.13</td><td>VGAM</td><td>vglm+summary</td><td>12.93</td><td>98.37x</td><td>2.47e-28</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalContRatioRegr</td><td>ordinal</td><td>0.06</td><td>VGAM</td><td>vglm+summary</td><td>12.11</td><td>209.36x</td><td>6.73e-43</td><td>***</td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalGCompMeanDiff</td><td>ordinal</td><td>0.31</td><td>ordinal</td><td>clm+gcomp+Wald</td><td>6.11</td><td>19.52x</td><td>1.76e-19</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceOrdinalJonckheereTerpstraTest</td><td>ordinal</td><td>9.10e-03</td><td>clinfun</td><td>jonckheere</td><td>0.36</td><td>39.56x</td><td>NA</td><td></td></tr>
    <tr style="background-color: #d9fdd3;"><td>InferenceOrdinalPropOddsRegr</td><td>ordinal</td><td>0.16</td><td>ordinal</td><td>clm+summary</td><td>4.97</td><td>31.12x</td><td>2.71e-32</td><td>***</td></tr>
    <tr style="background-color: #eceff1;"><td>InferenceOrdinalRidit</td><td>ordinal</td><td>8.90e-03</td><td>stats</td><td>mean(ridit)</td><td>0.12</td><td>13.99x</td><td>NA</td><td></td></tr>
  </tbody>
</table>

## Garbage Collection and Cache Management

To ensure that the benchmark results are highly precise, reproducible, and represent the actual computation speed of the numerical solvers, the benchmarking harness uses the following garbage collection and cache management strategies:

### 1. Garbage Collection (GC) Filtering
Garbage collection cycles run automatically by the R interpreter and can introduce significant, arbitrary pauses that skew timing measurements. To isolate the execution time of the code from R's GC overhead:
* **GC Disabling**: We disable R's memory stress-testing mode using `gctorture(FALSE)` before running timing loops.
* **Proactive Compaction**: In the `system.time()` path, we invoke `gc(verbose = FALSE)` immediately before timing each replicate. This starts the timer on a clean, compacted heap, minimizing the likelihood of triggering an automatic garbage collection cycle mid-replicate.
* **Automatic Filtering**: In the microbenchmarking path, we utilize the `bench::mark()` engine with the `filter_gc = TRUE` parameter, which automatically tracks and discards timing iterations during which a garbage collection event occurred.

### 2. Cold-Start Guarantee for EDI and Symmetric Warm-Up for Both Sides
Both EDI and canonical timing expressions receive a single **validation/warm-up call** executed once before the calibration loop begins. This puts the machine code and working data into the instruction and data caches in the same warmed state for both sides, so the official timed replicates start on equal footing.

EDI timings call exported C++ functions directly — no R6 objects are instantiated during benchmarking. As a result, **no R6 result caches exist to manage**. Each call to the C++ solver (e.g. `fast_logistic_regression_cpp`, `fast_ordinal_regression_cpp`) starts from a freshly zero-initialized parameter vector (or a model-specific data-driven initialization when `smart_cold_start = TRUE`). No prior-fit results are carried across timing repetitions, so every replication is a genuine cold start for the numerical optimizer.

<style>
    body, .markdown-body, .container {
        max-width: 1200px !important;
        width: 100% !important;
        margin: 0 auto !important;
    }
</style>
