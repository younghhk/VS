fit_full <- glmnet(X, Y, alpha = alpha_value)  # no CV, just the path

# For each lambda, how many vars are non-zero?
n_nonzero <- apply(coef(fit_full)[-1, ], 2, function(b) sum(b != 0))

# Choose lambda where number of selected vars is about q
q_target <- 20
idx <- which.min(abs(n_nonzero - q_target))
lambda_stab <- fit_full$lambda[idx]
