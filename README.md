# Variable Selection

## 1. LASSO vs Elastic Net: what they do and when to use them

Penalized regression adds a penalty to the usual loss to shrink coefficients and perform variable selection.

- **LASSO (L1 penalty)** solves something like  

  `min_beta { (1/(2n)) * ||Y - X beta||^2 + lambda * sum_j |beta_j| }`

  Main properties:

  - Sets many coefficients `beta_j` exactly to 0 → automatic variable selection.  
  - Works well when only a **small subset** of predictors are truly important.  
  - When predictors are **highly correlated**, it may pick one and drop the rest somewhat arbitrarily.

- **Elastic Net (L1 + L2 penalty)** solves  

  `min_beta { (1/(2n)) * ||Y - X beta||^2`
  `          + lambda * ( alpha * sum_j |beta_j| + (1 - alpha)/2 * sum_j beta_j^2 ) }`

  Main properties:

  - Combines LASSO (L1) and ridge (L2).  
  - Better behaved when predictors are **correlated** (e.g., groups of biomarkers, co-exposures):  
    tends to keep **groups** of related variables rather than forcing a single winner.  
  - Still performs variable selection as long as `alpha > 0`.

(Instability and sensitivity to tuning are common for both methods; see the *Stability selection* section for a way to address this.)



## 2. What do lambda and alpha do, and how are they chosen?

### Penalty strength: `lambda`

- `lambda >= 0` controls how strongly coefficients are shrunk.

  - Larger `lambda` → **more shrinkage**, more coefficients driven to 0 (simpler model).  
  - Smaller `lambda` → **less shrinkage**, more variables kept (higher risk of overfitting).


### Mixing parameter: `alpha`

- `alpha` controls the mix of L1 and L2 penalty:

  - `alpha = 1` → pure LASSO (L1 only).  
  - `alpha = 0` → pure ridge (L2 only; no variable selection).  
  - `0 < alpha < 1` → elastic net.

### Choosing `lambda` and `alpha` in practice

- Typically use **K-fold cross-validation**:

  - Pick a grid of `lambda` values (and sometimes `alpha` values).  
  - Fit models across folds and choose  

    - `lambda_min`: the value that gives minimum cross-validated error, or  
    - `lambda_1se`: the largest `lambda` whose error is within 1 standard error of the minimum (simpler, more regularized model).

- For `alpha`:

  - Often fix a small set (e.g., 0.25, 0.5, 0.75, 1) and choose by cross-validation, or  
  - Set based on scientific context (e.g., `alpha ~ 0.5` when strong collinearity is expected).



## 3. Why LASSO / Elastic Net over-shrink and sometimes select “too many” variables

Because these methods minimize a **penalized** objective,

`(1/(2n)) * ||Y - X beta||^2 + penalty_lambda(alpha, beta)`,

all coefficients are pulled toward zero:

- Even truly important predictors are **biased downward**.  
- For LASSO, the soft-thresholding formula shows this explicitly: every non-zero coefficient is reduced by roughly `lambda` (in the orthogonal case).  
- For elastic net, the additional L2 term spreads and smooths effects across correlated variables, often making individual coefficients smaller.

From theory, for **consistent variable selection**, LASSO needs:

- Reasonable conditions on the design matrix (so noise variables cannot mimic true ones), and  
- A **minimal signal strength**: the smallest true non-zero coefficient must be large enough relative to the penalty and noise level.

A common way to write this (informally) is:

`min_{j in S} |beta_j*| >= K * lambda`,

where:

- `S` is the set of truly non-zero predictors,  
- `beta_j*` is the true effect for variable `j`, and  
- `K > 0` is some constant.

**Interpretation:**

- If a true effect is weaker than a threshold on the order of `lambda`, LASSO will often shrink it to 0 and **fail to select it**.  
- In high-dimensional problems (large `p`), the “reasonable” `lambda` from theory is often of order `sigma * sqrt(log p / n)`, which can be fairly large.  
- With strong collinearity, true signal can be spread across groups of correlated variables, so LASSO / elastic net may:

  - keep **too many** variables (false positives), and  
  - assign each a **small, shrunken coefficient**.

**Practical takeaway:**  
LASSO and elastic net are very useful as **screening tools** for narrowing down a large set of predictors, but:

- effect sizes are biased toward zero, and  
- the exact selected set can be sensitive to tuning and sampling.

For interpretability and inference, it is often better to use LASSO / elastic net (possibly with stability selection) to pick a set of candidate variables
 
---

#  Stability Selection (to reduce sensitivity of LASSO / Elastic Net)

Penalized methods like LASSO and elastic net are powerful, but the **set of selected variables can be unstable**:

- Small changes in the data (adding or removing a few subjects),
- Small changes in the tuning parameter $\lambda$,
- Strong correlations among predictors,

can all lead to **different sets of selected variables**, even when prediction error is similar. This is a problem if we care about **interpretable biomarkers or exposures**, not just prediction accuracy.

Stability selection is a way to **stabilize** variable selection by combining subsampling with LASSO/elastic net. The idea goes back to Meinshausen & Bühlmann (2010).

### 1 Why LASSO / Elastic Net can be unstable

LASSO and elastic net choose coefficients $\beta$ by minimizing a penalized loss, for example

- loss: how well the model fits the data (e.g., least squares or partial likelihood),
- plus a penalty term that shrinks coefficients toward zero.

Because the penalty boundary is sharp:

- With **highly correlated variables**, LASSO may pick **one** variable in a correlated group somewhat arbitrarily.
- In high-dimensional settings (large $p$, modest $n$), small noise fluctuations can push a variable just above or just below the effective threshold.
- Different cross-validation splits or slightly different $\lambda$ values can lead to different selected models.

Result: two analysts, using similar choices, may get **different variable sets**, even if both models predict reasonably well.

### 2 Basic idea of stability selection

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

### 3 Error control

The original stability selection theory links:

- the threshold `pi_thr`,
- the typical number of variables selected in each subsample, and
- the expected number of false positives `E(V)` (sometimes described as a type of family-wise error control).

In plain language:

- If we choose a reasonably high selection-probability threshold and avoid very complex models in each subsample, we can bound the expected number of spurious variables.
- Stability selection is therefore not just a heuristic; it has a theoretical connection to controlling false discoveries.

