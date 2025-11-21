
#  Variable Selection and Prediction Pipeline



1. Step 1: Selection of alpha and lambda for penalized regression
2. Step 2: Stability selection to identify reproducible metabolites
3. Step 3: Linear modeling using the stable metabolites
4. Step 4: Prediction of metabolite based scores for new data

Each step is summarized below.

---

## ⭐ Step 1: Selection of Alpha and Lambda

The first step identifies a reasonable combination of alpha and lambda to be used in stability selection. Alpha controls the elastic net mixing parameter. Lambda controls the amount of shrinkage.

The procedure:

1. Fit elastic net models across a grid of alpha values.
2. For each alpha, use cross validation to obtain the cross validated error curve.
3. Choose the alpha that gives the smallest cross validated error.
4. For this alpha, select lambda using the one standard error rule to obtain a stable estimate.
5. Check how many metabolites are selected at this alpha lambda combination.
6. If too few metabolites are selected, relax lambda to avoid selecting an empty model.

Confounders are included in the model and are kept unpenalized. This ensures proper epidemiologic adjustment.

The result of Step 1 is a pair of tuning values:
**best_alpha** and **best_lambda**.

These values define a reasonable amount of shrinkage before stability selection.

---

Below is a polished **GitHub-ready `README.md`** section for your Stability Selection step.
Math is formatted with `$$` for GitHub Markdown rendering.

Clear, concise, with no dashes.

You can paste directly into your repo.

---

# ⭐ Step 2: Stability Selection

Stability selection identifies metabolites that appear consistently across random subsamples of the data. It uses repeated subsampling and elastic net fitting with the tuning parameters from Step 1.

The procedure:

1. Draw many half subsamples of the data.
2. Fit an elastic net model to each subsample using the chosen alpha and lambda.
3. Record how often each metabolite is selected.
4. Estimate the selection probability for each metabolite.
5. Retain metabolites whose selection probability exceeds a user defined cutoff.

Stability selection introduces the **Per Family Error Rate** (PFER). PFER is the expected number of false positive metabolites among the final selected set. The control of PFER follows from a theoretical inequality.

Let

* (p) be the total number of candidate metabolites
* (q) be the average number of metabolites selected in each subsample
* (\pi_{\text{thr}}) be the selection probability cutoff
* (V) be the number of false positives

The stability selection inequality states

$$
\mathbb{E}[V] \le \frac{q^{2}}{2p - q,} \cdot \frac{1}{\pi_{\text{thr}}^{2}}.
$$

This bound connects the cutoff and the expected number of false positives.
The user specifies a target PFER and chooses a cutoff large enough so that

$$
\mathbb{E}[V] \le \text{PFER}.
$$

For example, PFER equal to 1 means the expected number of false positive metabolites in the final selection is no greater than one.

The output of Step 2 is the list of **stable metabolites** that satisfy the cutoff and the estimated selection probabilities for all metabolites. These metabolites are passed to Step 3 for final model fitting.



---

##  Step 3: Final Linear Model

The stable metabolites are used to fit a simple linear regression model on the training dataset. This step produces the coefficients needed for prediction.

Only the stable metabolites identified in Step 2 are included. Confounders are included as covariates in the model. This model produces a clean and interpretable set of coefficients for score computation.

The output of Step 3 is:

* the fitted linear model
* the coefficient vector
* the list of stable metabolites

These coefficients are used to generate prediction scores.

---

##  Step 4: Prediction

The final step computes prediction scores for both training and test datasets. The score is computed as:

$\text{intercept} + X\beta$

where $X$ contains the confounders and the stable metabolites in the correct order. The function returns both training and test predictions.

---



