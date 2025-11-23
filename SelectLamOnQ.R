## ------------------------------------------------------------
## Requirements
## ------------------------------------------------------------
# install.packages("glmnet")
library(glmnet)

## ------------------------------------------------------------
## Helper: compute PFER bound
## ------------------------------------------------------------
compute_pfer_bound <- function(p, q, pi_thr) {
  if (pi_thr <= 0.5) {
    return(Inf)  # theory requires pi_thr > 0.5
  }
  (q^2) / ((2 * pi_thr - 1) * p)
}

## ------------------------------------------------------------
## Helper: choose lambda for stability selection
## ------------------------------------------------------------
choose_lambda_for_stability <- function(dat,
                                        y_var,
                                        method      = c("lasso", "elastic_net"),
                                        alpha_enet  = 0.5,
                                        family      = "gaussian",
                                        q_target,
                                        standardize = TRUE) {
  method <- match.arg(method)
  
  # Outcome
  Y <- dat[[y_var]]
  
  # Predictors: all other columns
  x_vars <- setdiff(names(dat), y_var)
  X <- model.matrix(~ . - 1, data = dat[, x_vars, drop = FALSE])
  
  # Alpha: 1 for LASSO, alpha_enet for elastic net
  if (method == "lasso") {
    alpha_value <- 1
  } else {
    alpha_value <- alpha_enet
  }
  
  # Fit full glmnet path (no CV)
  fit_full <- glmnet(
    x          = X,
    y          = Y,
    alpha      = alpha_value,
    family     = family,
    standardize = standardize
  )
  
  # Count non-zero coefficients for each lambda (ignore intercept)
  beta_mat  <- coef(fit_full)[-1, , drop = FALSE]
  n_nonzero <- apply(beta_mat, 2, function(b) sum(b != 0))
  
  # Choose lambda where number of selected vars is closest to q_target
  idx <- which.min(abs(n_nonzero - q_target))
  lambda_stab <- fit_full$lambda[idx]
  
  list(
    lambda_stab   = lambda_stab,
    n_nonzero     = n_nonzero[idx],
    alpha         = alpha_value,
    method        = method,
    family        = family,
    fit_full      = fit_full,
    n_nonzero_all = n_nonzero,
    x_vars        = colnames(X)
  )
}

## ------------------------------------------------------------
## Main function: stability selection with optional PFER control
## ------------------------------------------------------------
stability_select <- function(dat,
                             y_var,
                             method         = c("lasso", "elastic_net"),
                             alpha_enet     = 0.5,
                             family         = "gaussian",
                             # EITHER set q_target OR pfer_target (not both)
                             q_target       = NULL,
                             pfer_target    = NULL,
                             pi_thr         = 0.6,
                             B              = 100,
                             subsample_frac = 0.5,
                             standardize    = TRUE,
                             seed           = NULL) {
  method <- match.arg(method)
  
  # Outcome and predictors
  Y <- dat[[y_var]]
  x_vars <- setdiff(names(dat), y_var)
  X <- model.matrix(~ . - 1, data = dat[, x_vars, drop = FALSE])
  
  n <- nrow(X)
  p <- ncol(X)
  
  # ----------------------------------------------------------
  # 0) Handle q_target vs pfer_target
  # ----------------------------------------------------------
  if (!is.null(q_target) && !is.null(pfer_target)) {
    stop("Specify either q_target OR pfer_target, not both.")
  }
  if (is.null(q_target) && is.null(pfer_target)) {
    stop("You must specify either q_target or pfer_target.")
  }
  
  # If user provided pfer_target, derive q_target from theory
  if (!is.null(pfer_target)) {
    if (pi_thr <= 0.5) {
      stop("For PFER control, pi_thr must be > 0.5.")
    }
    q_max <- floor(sqrt(pfer_target * (2 * pi_thr - 1) * p))
    if (q_max < 1) {
      warning("Computed q_target < 1 from PFER; setting q_target = 1.")
      q_target <- 1
    } else {
      q_target <- q_max
    }
  }
  
  # Compute implied PFER bound (even if user gave q_target directly)
  pfer_bound <- compute_pfer_bound(p = p, q = q_target, pi_thr = pi_thr)
  
  # ----------------------------------------------------------
  # 1) Choose lambda for stability selection (given q_target)
  # ----------------------------------------------------------
  lambda_info <- choose_lambda_for_stability(
    dat         = dat,
    y_var       = y_var,
    method      = method,
    alpha_enet  = alpha_enet,
    family      = family,
    q_target    = q_target,
    standardize = standardize
  )
  
  lambda_stab <- lambda_info$lambda_stab
  alpha_value <- lambda_info$alpha
  x_vars      <- lambda_info$x_vars
  
  # ----------------------------------------------------------
  # 2) Stability selection loop
  # ----------------------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)  # reproducible subsampling
  }
  
  n_sub <- floor(subsample_frac * n)
  p <- ncol(X)
  
  s_counts <- numeric(p)
  names(s_counts) <- x_vars
  
  for (b in seq_len(B)) {
    # Subsample indices
    idx_sub <- sample(seq_len(n), size = n_sub, replace = FALSE)
    
    X_sub <- X[idx_sub, , drop = FALSE]
    Y_sub <- Y[idx_sub]
    
    fit_sub <- glmnet(
      x          = X_sub,
      y          = Y_sub,
      alpha      = alpha_value,
      lambda     = lambda_stab,   # single lambda
      family     = family,
      standardize = standardize
    )
    
    beta_sub <- as.vector(coef(fit_sub)[-1, 1])  # drop intercept
    selected <- (beta_sub != 0)
    
    s_counts[selected] <- s_counts[selected] + 1
  }
  
  # ----------------------------------------------------------
  # 3) Selection probabilities and stable set
  # ----------------------------------------------------------
  pi_hat <- s_counts / B
  pi_hat_sorted <- sort(pi_hat, decreasing = TRUE)
  
  stable_vars <- names(pi_hat_sorted)[pi_hat_sorted >= pi_thr]
  
  list(
    stable_vars    = stable_vars,
    pi_hat         = pi_hat_sorted,
    lambda_stab    = lambda_stab,
    alpha          = alpha_value,
    method         = method,
    family         = family,
    q_target       = q_target,
    pi_thr         = pi_thr,
    B              = B,
    subsample_frac = subsample_frac,
    p              = p,
    pfer_bound     = pfer_bound    # theoretical upper bound on E[#false positives]
  )
}

## ------------------------------------------------------------
## Example usage
## ------------------------------------------------------------
# Assume:
#   dat: data.frame
#   Y: outcome column (continuous here)
#
# Example 1: Fix q_target = 20 (no explicit PFER)
# res1 <- stability_select(
#   dat           = dat,
#   y_var         = "Y",
#   method        = "lasso",
#   family        = "gaussian",
#   q_target      = 20,
#   pi_thr        = 0.6,
#   B             = 100,
#   subsample_frac = 0.5,
#   seed          = 123
# )
# res1$pfer_bound   # implied PFER bound
# res1$stable_vars
#
# Example 2: Specify PFER target instead of q_target
#   Target: E[#false positives] <= 1, pi_thr = 0.9
#
# res2 <- stability_select(
#   dat           = dat,
#   y_var         = "Y",
#   method        = "lasso",
#   family        = "gaussian",
#   pfer_target   = 1,
#   pi_thr        = 0.9,
#   B             = 100,
#   subsample_frac = 0.5,
#   seed          = 123
# )
# res2$q_target    # q derived from pfer_target
# res2$pfer_bound  # should be <= 1 (up to rounding)
# res2$stable_vars

