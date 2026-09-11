LDPE_func <- function(ID, X, OT, status) {
  coi <- ID
  x <- as.matrix(X)
  d <- ncol(x)
  n <- nrow(x)

  PF <- matrix(1, 1, d)
  PF[ID] <- 1

  fit <- glmnet::glmnet(
    x,
    survival::Surv(OT, status),
    family = "cox",
    alpha = 1,
    standardize = FALSE,
    penalty.factor = PF
  )
  cv.fit <- glmnet::cv.glmnet(
    x,
    survival::Surv(OT, status),
    family = "cox",
    alpha = 1,
    standardize = FALSE,
    penalty.factor = PF
  )
  betas <- as.numeric(stats::coef(fit, s = cv.fit$lambda.min))[seq_len(d)]

  stime <- sort(OT)
  otime <- order(OT)

  Hs <- matrix(0, nrow = d, ncol = d)
  la <- numeric(n)
  lb <- matrix(0, nrow = n, ncol = d - 1)

  i <- 1
  while (i <= n) {
    if (status[otime[i]] == 1) {
      ind <- which(OT >= stime[i])
      S0 <- 0
      S1 <- rep(0, d)
      S2 <- matrix(0, nrow = d, ncol = d)

      if (length(ind) > 0) {
        for (j in ind) {
          tmp <- exp(x[j, ] %*% betas)
          S0 <- S0 + tmp
          S1 <- S1 + tmp %*% t(x[j, ])
          tmp_num <- as.numeric(tmp)
          S2 <- S2 + tmp_num * x[j, ] %*% t(x[j, ])
        }
      }

      S0 <- as.numeric(S0)
      la[i] <- -(x[otime[i], coi] - S1[coi] / S0)

      nuisance_idx <- setdiff(seq_len(d), coi)
      lb[i, ] <- -(x[otime[i], nuisance_idx] - S1[nuisance_idx] / S0)
      V <- S0 * S2 - t(S1) %*% S1
      Hs <- Hs + V / (n * S0^2)
    }
    i <- i + 1
  }

  fit <- glmnet::glmnet(
    lb,
    la,
    alpha = 1,
    standardize = FALSE,
    intercept = FALSE,
    lambda = sqrt(log(d) / n)
  )
  what <- as.numeric(stats::coef(fit)[-1])
  nuisance_idx <- setdiff(seq_len(d), coi)

  correction <- mean(la) - t(what) %*% colMeans(lb)
  variance_term <- Hs[coi, coi] - t(what) %*% Hs[nuisance_idx, coi]
  S <- betas[coi] - correction / variance_term
  var_est <- variance_term

  c(est = as.numeric(S), se = as.numeric(sqrt(1 / (n * var_est))))
}
