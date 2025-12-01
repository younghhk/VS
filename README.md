[![Cancer Research Software Hub](https://img.shields.io/badge/Back_to-Hub-blue)](https://github.com/younghhk/NCI)

# Variable Selection and Stability in High-Dimensional Data

*Using LASSO, Elastic Net, and Stability Selection*

---

## 1. What LASSO and Elastic Net Do

Penalized regression shrinks coefficients and selects variables. This is useful when the number of predictors is large or predictors are correlated.

### LASSO (L1 penalty)

* Performs automatic variable selection.
* Effective when only a small number of predictors matter.
* Can behave unpredictably when predictors are correlated (may select one and drop others arbitrarily).

### Elastic Net (L1 + L2 penalty)

* More stable when predictors are correlated.
* Tends to keep groups of related metabolites or exposures.
* Still performs variable selection as long as some L1 penalty is used.

### Choosing tuning parameters

* Tuning is done via K-fold cross-validation.
* Two common choices:

  * `lambda_min`: gives the lowest prediction error.
  * `lambda_1se`: simpler model; better when interpretability is a priority.
* `alpha` controls the mix of L1 vs L2; values around 0.5 work well in highly correlated data.

---

## 2. Why Penalized Methods Can Be Unstable

LASSO and elastic net are powerful, but the selected variables can change when:

* a few participants are added or removed,
* tuning parameters change slightly, or
* predictors are highly correlated.

Because all coefficients are shrunk toward zero, weak but true signals may be missed, and correlated predictors may be selected inconsistently.



---

## 3. Stability Selection

*A method to identify variables that are selected consistently*

Stability selection improves reproducibility by combining subsampling with penalized regression.

Process:

1. Take many random subsamples of the data.
2. Fit LASSO or elastic net to each subsample.
3. Count how often each variable is selected.

Variables with high selection frequency are more likely to be robust, biologically meaningful predictors.

### Why it helps

* Reduces sensitivity to sampling noise.
* Prioritizes reproducible biomarkers.
* More reliable for interpretation than a single LASSO or elastic net fit.


### Reference

Meinshausen N, Bühlmann P. *Stability selection.* **Journal of the Royal Statistical Society: Series B.** 2010;72(4):417-473.
[https://doi.org/10.1111/j.1467-9868.2010.00740.x](https://doi.org/10.1111/j.1467-9868.2010.00740.x)


---

## Group LASSO: Brief Overview

Group LASSO extends LASSO to situations where predictors naturally belong to **groups**—for example:

* metabolites within the same biochemical pathway
* SNPs within a gene
* repeated measures of the same exposure
* correlated biomarkers from the same assay platform

### What it does

Instead of penalizing each coefficient separately, group LASSO applies a penalty to the **entire group** of coefficients. As a result:

* **Groups are selected or removed as a whole**
* Works well when predictors within a group are **correlated**
* Encourages **biologically coherent** variable selection

Mathematically, the penalty is:

$$
\lambda \sum_{g=1}^{G} |\beta_g|_2,
$$

which shrinks an entire group of coefficients to zero if the evidence for the group is weak.

### When it is useful in epidemiologic studies

Group LASSO is helpful when:

* scientific knowledge defines meaningful groups (pathways, genes, chemical classes)
* strong correlation occurs within groups
* interest lies in **process-level** or **pathway-level** signals, not individual variables
* you want more stable selection than standard LASSO, which may pick one variable and drop the rest

### Limitations

* All-or-none selection: it does **not** pick individual predictors within a group
* Requires pre-defined groups
* Coefficients remain biased toward zero, so refitting models after selection may help with interpretation


---

