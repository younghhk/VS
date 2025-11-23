# You need glmnet
library(glmnet)

choose_lambda_for_stability <- function(dat,
                                        y_var,
                                        method = c("lasso", "elastic_net"),
                                        alpha_enet = 0.5,
                                        family = "gaussian",
                                        q_target = 20) {
  # 1. Match method argument
  method <- match.arg(method)
  
  # 2. Extract outcome Y
  Y <- dat[[y_var]]
  
  # 3. Build predictor matrix X (all other columns)
  x_vars <- setdiff(names(dat), y_var)
  
  # model.matrix() creates a numeric design matrix and handles factors
  X <- model.matrix(~ . - 1, data = dat[, x_vars, drop = FALSE])
  #   ~ . - 1  means: use all variables, no intercept (glmnet adds its own)
  
  # 4. Set alpha:
  #    - LASSO: alpha = 1
  #    - Elastic net: alpha between 0 and 1 (e.g., 0.5)
  if (method == "lasso") {
    alpha_value <- 1
  } else {
    alpha_value <- alpha_enet  # e.g., 0.5
  }
  
  # 5. Fit full glmnet path (no cross-validation)
  #    family:
  #      "gaussian"  -> continuous outcome
  #      "binomial"  -> binary outcome
  #      "cox"       -> survival outcome (with Surv object, not handled here)
  fit_full <- glmnet(X, Y,
                     alpha  = alpha_value,
                     family = family,
                     standardize = TRUE)
  
  # 6. Count number of non-zero coefficients at each lambda
  # coef(fit_full) returns a matrix: rows = coefficients, cols = lambda values
  # We drop the intercept row (first row)
  beta_mat <- coef(fit_full)[-1, , drop = FALSE]
  
  n_nonzero <- apply(beta_mat, 2, function(b) sum(b != 0))
  
  # 7. Choose lambda where the number of selected vars is closest to q_target
  idx <- which.min(abs(n_nonzero - q_target))
  lambda_stab <- fit_full$lambda[idx]
  
  # 8. Return useful info
  list(
    lambda_stab  = lambda_stab,
    n_nonzero    = n_nonzero[idx],
    alpha        = alpha_value,
    method       = method,
    family       = family,
    fit_full     = fit_full,
    n_nonzero_all = n_nonzero
  )
}

           
