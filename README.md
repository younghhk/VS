## 2. Stability Selection (to reduce sensitivity of LASSO / Elastic Net)

Penalized methods like LASSO and elastic net are powerful, but the **set of selected variables can be unstable**:

- Small changes in the data (adding or removing a few subjects),
- Small changes in the tuning parameter $\lambda$,
- Strong correlations among predictors,

can all lead to **different sets of selected variables**, even when prediction error is similar. This is a problem if we care about **interpretable biomarkers or exposures**, not just prediction accuracy.

Stability selection is a way to **stabilize** variable selection by combining subsampling with LASSO/elastic net. The idea goes back to Meinshausen & Bühlmann (2010).

### 2.1 Why LASSO / Elastic Net can be unstable

LASSO and elastic net choose coefficients $\beta$ by minimizing a penalized loss, for example

- loss: how well the model fits the data (e.g., least squares or partial likelihood),
- plus a penalty term that shrinks coefficients toward zero.

Because the penalty boundary is sharp:

- With **highly correlated variables**, LASSO may pick **one** variable in a correlated group somewhat arbitrarily.
- In high-dimensional settings (large $p$, modest $n$), small noise fluctuations can push a variable just above or just below the effective threshold.
- Different cross-validation splits or slightly different $\lambda$ values can lead to different selected models.

Result: two analysts, using similar choices, may get **different variable sets**, even if both models predict reasonably well.

### 2.2 Basic idea of stability selection

Stability selection reduces this sensitivity by repeatedly fitting the model on random subsamples and keeping only variables that are selected consistently.

Basic algorithm:

1. Choose a penalized method (e.g., LASSO or elastic net) and a grid of lambda values.
2. For `b = 1, ..., B` (e.g., `B = 50–100`):
   - Draw a random subsample of the data (commonly about 50–60% of subjects, without replacement).
   - Fit LASSO / elastic net on this subsample.
   - Record which variables are selected (non-zero coefficients).
3. For each variable `j`, compute its selection probability as  

   `pi_hat_j = s_j / B`,  

   where `s_j` is the number of subsamples in which variable `j` was selected and `B` is the total number of subsamples.
4. Choose a threshold `pi_thr` (for example, `pi_thr = 0.6`) and keep variables with  

   `pi_hat_j >= pi_thr`.

Interpretation:

- A variable is kept only if it is selected in a clear majority of subsamples.
- Variables that appear only occasionally (e.g., 5–10% of the time) are treated as unstable and dropped.
- Using a threshold around 0.6 is somewhat arbitrary but a reasonable starting point; higher thresholds give fewer, more conservative variables.

### 2.3 Error control

The original stability selection theory links:

- the threshold `pi_thr`,
- the typical number of variables selected in each subsample, and
- the expected number of false positives `E(V)` (sometimes described as a type of family-wise error control).

In plain language:

- If we choose a reasonably high selection-probability threshold and avoid very complex models in each subsample, we can bound the expected number of spurious variables.
- Stability selection is therefore not just a heuristic; it has a theoretical connection to controlling false discoveries.

