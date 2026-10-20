## New simulations LSPIM

library(multcomp)
library(nparLD)
library(pim)
library(tidyverse)
library(aod)
library(geessbin)
library(mvtnorm)
library(dplyr)


## Functions
geessbin_adjusted = function (formula, data = parent.frame(), id = NULL, corstr = "independence", 
                              repeated = NULL, beta.method = "PGEE", SE.method = "MB", 
                              b = NULL, maxitr = 50, tol = 1e-05, scale.fix = FALSE, conf.level = 0.95) 
{
  Call <- match.call()
  mc <- match.call(expand.dots = FALSE)
  mc$corstr <- mc$beta.method <- mc$SE.method <- mc$b <- mc$maxitr <- mc$tol <- mc$scale.fix <- mc$conf.level <- NULL
  mc[[1]] <- as.name("model.frame")
  dat <- eval(mc, parent.frame())
  id <- model.extract(dat, "id")
  repeated <- model.extract(dat, "repeated")
  names(id) <- names(repeated) <- NULL
  if (is.null(id) & !is.null(repeated)) {
    stop("'id' must be specified when 'repeated' is not NULL")
  }
  if (is.null(id) & is.null(repeated)) {
    message(paste("'id' and 'repeated' are not specified", 
                  "\n", "all observations are assumed to be independent"))
    idseq <- 1:nrow(dat)
    repval <- NULL
    repseq <- rep(1, nrow(dat))
  }
  if (!is.null(id) & is.null(repeated)) {
    id <- deparse(substitute(id))
    idval <- dat[, "(id)"]
    chg <- (1:length(idval))[c(TRUE, idval[-length(idval)] != 
                                 idval[-1])]
    nidat <- c(chg[-1], length(idval) + 1) - chg
    idseq <- rep(1:length(nidat), time = nidat)
    repseq <- unlist(tapply(nidat, unique(idseq), function(x) 1:x))
    names(repseq) <- repval <- NULL
  }
  if (!is.null(id) & !is.null(repeated)) {
    dat <- dat[order(id, repeated), ]
    idseq <- as.numeric(factor(dat[, "(id)"]))
    repval <- dat[, "(repeated)"]
    repseq <- as.numeric(factor(dat[, "(repeated)"]))
  }
  n <- length(unique(repseq))
  K <- length(unique(idseq))
  ndat <- as.numeric(table(idseq))
  replst <- split(repseq, idseq)
  Terms <- attr(dat, "terms")
  y <- as.matrix(model.extract(dat, "response"))
  X <- model.matrix(Terms, dat)
  p <- ncol(X)
  if (!is.numeric(y) || anyNA(y) ||
      !all(y %in% c(0, 0.5, 1))) {
    stop("outcome vector must be numeric and take values in {0, 0.5, 1}")
  }
  for (v in c("corstr", "beta.method", "SE.method")) {
    if (eval(parse(text = paste0("length(", v, ")"))) > 1) {
      stop(paste0("'", v, "'", " has length > 1"))
    }
  }
  corstrs <- c("independence", "exchangeable", "ar1", "unstructured")
  if (is.na(match(corstr, corstrs))) {
    stop(paste(c("invalid correlation structure", "\n", "'corstr' must be specified from the following list:", 
                 "\n", paste(paste0("\"", corstrs, "\""), collapse = ", "))))
  }
  beta.methods <- c("GEE", "BCGEE", "PGEE")
  if (is.na(match(beta.method, beta.methods))) {
    stop(paste(c("invalid estimation method", "\n", "'beta.method' must be specified from the following list:", 
                 "\n", paste(paste0("\"", beta.methods, "\""), collapse = ", "))))
  }
  SE.methods <- c("SA", "MK", "KC", "MD", "FG", "PA", "GS", 
                  "MB", "WL", "WB", "FW", "FZ")
  if (is.na(match(SE.method, SE.methods))) {
    stop(paste(c("invalid SE estimator", "\n", "'SE.method' must be specified from the following list:", 
                 "\n", paste(paste0("\"", SE.methods, "\""), collapse = ", "))))
  }
  if (!is.null(b) & length(b) != p) {
    stop(paste("length of 'b' must be ncol(X) =", p))
  }
  comp <- unlist(lapply(replst, function(x) length(x) == n))
  if (sum(!comp) > 0) {
    if (!is.na(match(SE.method, c("PA", "GS", "WL", "WB")))) {
      stop(paste0("\"", SE.method, "\"", " method cannot be used for incomplete data"))
    }
  }
  if (conf.level <= 0 | conf.level >= 1) {
    stop("'conf.level' must be in interval (0,1)")
  }
  if (is.null(b)) {
    if (beta.method == "PGEE") {
      b <- numeric(p)
      del <- 100
      nitr <- 0
      while (del > 1e-05) {
        mu <- 1/(1 + exp(-X %*% b))
        I <- t(c(mu * (1 - mu)) * X) %*% X
        U <- t(X) %*% (y - mu + diag(X %*% ginv(I) %*% 
                                       t(X)) * mu * (1 - mu) * (0.5 - mu))
        del <- max(abs(U))
        if (del > 1e-05) 
          b <- b + ginv(I) %*% U
        nitr <- nitr + 1
        if (nitr == 50) 
          break
      }
    }
    else {
      b <- glm.fit(X, y, family = binomial(link = "logit"))$coefficients
    }
  }
  else {
    if (anyNA(b) | sum(is.infinite(b)) > 0) 
      stop("'b' contains Na/NaN/Inf")
    names(b) <- colnames(X)
  }
  conv <- "converged"
  nitr <- 0
  del <- 100
  while (del > tol) {
    mu <- 1/(1 + exp(-X %*% b))
    r <- (y - mu)/sqrt(mu * (1 - mu))
    if (min(mu) < 1e-04 | max(mu) > 0.9999) {
      conv <- "fitted probabilities numerically 0 or 1 occurred."
      warning(conv)
      break
    }
    if (scale.fix == TRUE) 
      phi <- 1
    if (scale.fix == FALSE) 
      phi <- sum(r^2)/(sum(ndat) - p)
    if (is.infinite(phi)) {
      conv <- "infinite scale parameter"
      warning(conv)
      break
    }
    if (corstr == "independence") 
      R <- diag(n)
    if (corstr == "exchangeable") {
      a0 <- 0
      for (i in 1:K) {
        ri <- r[idseq == i]
        pmat <- tcrossprod(ri)
        a0 <- a0 + sum(pmat[upper.tri(pmat)])
      }
      alpha <- a0/((0.5 * sum(ndat * (ndat - 1)) - p) * 
                     phi)
      R <- matrix(alpha, n, n) + diag(1 - alpha, n, n)
    }
    if (corstr == "ar1") {
      a0 <- d0 <- 0
      for (i in 1:K) {
        ti <- replst[[i]]
        ri <- numeric(n)
        ri[replst[[i]]] <- r[idseq == i]
        a0 <- a0 + sum(ri[-1] * ri[-n])
        d0 <- d0 + sum((ti[-1] - ti[-length(ti)]) == 
                         1)
      }
      alpha <- a0/((d0 - p) * phi)
      R <- alpha^abs(matrix(0:(n - 1), nrow = n, ncol = n, 
                            byrow = TRUE) - 0:(n - 1))
    }
    if (corstr == "unstructured") {
      m <- count <- matrix(0, n, n)
      for (i in 1:K) {
        ri <- ci <- numeric(n)
        ri[replst[[i]]] <- r[idseq == i]
        ci[replst[[i]]] <- 1
        m <- m + tcrossprod(ri)
        count <- count + tcrossprod(ci)
      }
      R <- m/(phi * (count - p))
      diag(R) <- 1
    }
    U <- numeric(p)
    I <- matrix(0, p, p)
    dI <- array(0, c(p, p, p))
    for (i in 1:K) {
      mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == 
                                                         i], b, R[replst[[i]], replst[[i]]], phi)
      U <- U + t(mat$VD) %*% mat$e
      I <- I + t(mat$D) %*% mat$VD
      if (beta.method == "PGEE") {
        dI <- dI + array(apply(X[idseq == i, , drop = FALSE], 
                               2, function(x) {
                                 t((c(1 - 2 * mat$mu) * x) * mat$D) %*% mat$VD
                               }), c(p, p, p))
      }
    }
    Iinv <- ginv(I)
    if (beta.method == "PGEE") {
      U <- U + 0.5 * apply(dI, 3, function(x) sum(diag(Iinv %*% 
                                                         x)))
    }
    del <- max(abs(U))
    if (del > tol) 
      b <- b + Iinv %*% U
    nitr <- nitr + 1
    if (nitr == maxitr) {
      if (del > tol) 
        conv <- "maximum number of iterations consumed"
      warning(conv)
      break
    }
  }
  if (conv == "converged" & del > tol) {
    conv <- "convergence failure"
    warning(conv)
  }
  if (conv == "converged") {
    if (beta.method == "BCGEE") {
      k11 <- array(0, c(p, p))
      k21 <- k3 <- array(0, c(p, p, p))
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        k11 <- k11 + t(mat$VD) %*% mat$emat %*% mat$VD
        px <- (1 - 2 * mat$mu) * X[idseq == i, ]
        dD <- array(matrix(rep(t(px), p), ncol = p, byrow = TRUE) * 
                      rep(mat$D, p), c(ndat[i], p, p))
        dDV <- array(0.5 * t(mat$D) %*% matrix((matrix(rep(t(px), 
                                                           ndat[i]), ncol = p, byrow = TRUE) - rep(px, 
                                                                                                   each = ndat[i])) * rep(mat$Vinv, p), nrow = ndat[i]), 
                     c(p, ndat[i], p))
        for (u in 1:p) {
          k21[, u, ] <- k21[, u, ] + dDV[, , u] %*% mat$emat %*% 
            mat$VD
          dDV_D <- dDV[, , u] %*% mat$D
          k3[, u, ] <- k3[, u, ] - dDV_D
          k3[, , u] <- k3[, , u] - dDV_D - t(mat$VD) %*% 
            dD[, , u]
        }
      }
      bhat0 <- numeric(p)
      for (u in 1:p) {
        bhat0 <- bhat0 + (k21[, , u] + 0.5 * k3[, , u] %*% 
                            Iinv %*% k11) %*% Iinv[, u]
      }
      b <- b - Iinv %*% bhat0
      I <- matrix(0, p, p)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        I <- I + t(mat$D) %*% mat$VD
      }
      Iinv <- ginv(I)
    }
    J <- matrix(0, p, p)
    if (SE.method == "SA" | SE.method == "MK") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        J <- J + t(mat$VD) %*% mat$emat %*% mat$VD
      }
      if (SE.method == "MK") 
        J <- J * K/(K - p)
    }
    if (SE.method == "KC") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        HKC <- sqrtmat(ginv(diag(ndat[i]) - mat$D %*% 
                              Iinv %*% t(mat$VD)))
        J <- J + t(mat$VD) %*% HKC %*% mat$emat %*% t(HKC) %*% 
          mat$VD
      }
    }
    if (SE.method == "MD") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        HMD <- ginv(diag(ndat[i]) - mat$D %*% Iinv %*% 
                      t(mat$VD))
        J <- J + t(mat$VD) %*% HMD %*% mat$emat %*% t(HMD) %*% 
          mat$VD
      }
    }
    if (SE.method == "FG") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        Fi <- diag((1 - pmin(0.75, diag(t(mat$VD) %*% 
                                          mat$D %*% Iinv)))^(-0.5), p, p)
        J <- J + Fi %*% t(mat$VD) %*% mat$emat %*% mat$VD %*% 
          Fi
      }
    }
    if (SE.method == "PA" | SE.method == "GS") {
      M <- matrix(0, n, n)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        M <- M + sqrt(1/mat$nu) * mat$emat * rep(sqrt(1/mat$nu), 
                                                 each = ndat[i])
      }
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        J <- J + t(mat$VD) %*% (sqrt(mat$nu) * M) %*% 
          (sqrt(mat$nu) * mat$VD)
      }
      if (SE.method == "PA") 
        J <- J/K
      if (SE.method == "GS") 
        J <- J/(K - p)
    }
    if (SE.method == "MB") {
      d <- matrix(0, K, p)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        d[i, ] <- t(mat$VD) %*% mat$e
      }
      I1 <- (sum(ndat) - 1) * K * cov(d)/(sum(ndat) - p)
      q <- min(0.5, p/(K - p)) * max(1, sum(diag(Iinv %*% 
                                                   I1))/p)
      J <- I1 + q * I
    }
    if (SE.method == "WL") {
      M <- matrix(0, n, n)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        HMD <- ginv(diag(ndat[i]) - mat$D %*% Iinv %*% 
                      t(mat$VD))
        M <- M + (sqrt(1/mat$nu) * HMD) %*% mat$emat %*% 
          t((sqrt(1/mat$nu) * HMD))
      }
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        J <- J + t(mat$VD) %*% (sqrt(mat$nu) * M) %*% 
          (sqrt(mat$nu) * mat$VD)/K
      }
    }
    if (SE.method == "WB") {
      M <- matrix(0, n, n)
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        HKC <- sqrtmat(ginv(diag(ndat[i]) - mat$D %*% 
                              Iinv %*% t(mat$VD)))
        M <- M + (sqrt(1/mat$nu) * HKC) %*% mat$emat %*% 
          t((sqrt(1/mat$nu) * HKC))
      }
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        J <- J + t(mat$VD) %*% (sqrt(mat$nu) * M) %*% 
          (sqrt(mat$nu) * mat$VD)/K
      }
    }
    if (SE.method == "FW") {
      for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], 
                        y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                        phi)
        Hi <- mat$D %*% Iinv %*% t(mat$VD)
        HKC <- sqrtmat(ginv(diag(ndat[i]) - Hi))
        HMD <- ginv(diag(ndat[i]) - Hi)
        J <- J + 0.5 * t(mat$VD) %*% (HKC %*% mat$emat %*% 
                                        t(HKC) + HMD %*% mat$emat %*% t(HMD)) %*% mat$VD
      }
    }
    if (SE.method == "FZ") {
      for (i in 1:K) {
        mati <- calc_mat(X[idseq == i, , drop = FALSE], 
                         y[idseq == i], b, R[replst[[i]], replst[[i]]], 
                         phi)
        HMD <- ginv(diag(ndat[i]) - mati$D %*% Iinv %*% 
                      t(mati$VD))
        M <- matrix(0, p, p)
        for (j in setdiff(1:K, i)) {
          matj <- calc_mat(X[idseq == j, , drop = FALSE], 
                           y[idseq == j], b, R[replst[[j]], replst[[j]]], 
                           phi)
          M <- M + t(matj$VD) %*% matj$emat %*% matj$VD
        }
        J <- J + t(mati$VD) %*% HMD %*% (mati$emat - 
                                           mati$D %*% Iinv %*% M %*% Iinv %*% t(mati$D)) %*% 
          t(HMD) %*% mati$VD
      }
    }
    covb <- Iinv %*% J %*% Iinv
  }
  else {
    covb <- matrix(NA, p, p)
  }
  if (!exists("phi")) 
    phi <- NA
  if (!exists("b")) 
    b <- rep(NA, p)
  if (!exists("R")) 
    R <- matrix(NA, n, n)
  lin <- c(X %*% b)
  mu <- c(1/(1 + exp(-lin)))
  resid <- c(y - mu)
  b <- as.vector(b)
  names(b) <- colnames(X)
  if (is.null(repval)) {
    colnames(R) <- rownames(R) <- NULL
  }
  else {
    rep_unique <- unique(data.frame(repseq = repseq, repval = repval))
    rep_unique <- rep_unique[order(rep_unique$repseq), ]
    colnames(R) <- rownames(R) <- rep_unique$repval
  }
  structure(class = "geessbin", list(call = Call, coefficients = b, 
                                     linear.predictors = lin, fitted.values = mu, residuals = resid, 
                                     scale = phi, covb = covb, wcorr = R, iterations = nitr, 
                                     beta.method = beta.method, SE.method = SE.method, K = K, 
                                     max.ni = n, corstr = corstr, convergence = conv, conf.level = conf.level, 
                                     model.matrix = X, data = data))
}

require(geessbin)
environment(geessbin_adjusted) <- asNamespace('geessbin')
assignInNamespace("geessbin", geessbin_adjusted, ns = "geessbin")


lachin_constancy_test <- function(
    data,
    group_col,
    time_cols,
    alpha = 0.05
) {
  
  K <- length(time_cols)
  
  if (K < 2L)
    stop("At least two visits are required.")
  
  if (!all(c(group_col, time_cols) %in% names(data)))
    stop("Requested columns are missing.")
  
  group <- data[[group_col]]
  
  if (anyNA(group) || !all(group %in% c(0, 1)))
    stop("Group must contain only 0 and 1.")
  
  if (!all(vapply(data[time_cols], is.numeric, logical(1))))
    stop("Visit outcomes must be numeric.")
  
  Y <- as.matrix(data[time_cols])
  
  if (any(!is.finite(Y)))
    stop("Complete, finite outcomes are required.")
  
  # Population 1 = control; population 2 = treatment
  Y0 <- Y[group == 0, , drop = FALSE]
  Y1 <- Y[group == 1, , drop = FALSE]
  
  n0 <- nrow(Y0)
  n1 <- nrow(Y1)
  N <- n0 + n1
  
  if (min(n0, n1) < 2L)
    stop("At least two subjects per group are required.")
  
  # Subject contributions from equations (33)–(34)
  contribution <- function(own, other) {
    
    # Y_i(x_j): number of own-group outcomes >= x_j
    own_ge <- rowSums(
      outer(own, own, FUN = function(x, y) y >= x)
    )
    
    # Y_other(x_j)
    other_ge <- rowSums(
      outer(own, other, FUN = function(x, y) y >= x)
    )
    
    # mu_i(x_j)
    mu <- other_ge / N
    
    # eta_i(x_j): sum over own-group observations <= x_j
    own_le <- outer(
      own, own,
      FUN = function(x, y) y <= x
    )
    
    eta <- drop(own_le %*% (mu / own_ge))
    
    mu - eta
  }
  
  H0 <- matrix(NA_real_, n0, K)
  H1 <- matrix(NA_real_, n1, K)
  theta <- numeric(K)
  
  for (k in seq_len(K)) {
    
    # Equations (30)–(31), including ties
    theta[k] <- mean(outer(
      Y1[, k],
      Y0[, k],
      FUN = function(treated, control) sign(treated - control)
    ))
    
    H0[, k] <- contribution(Y0[, k], Y1[, k])
    H1[, k] <- contribution(Y1[, k], Y0[, k])
  }
  
  # Equation (33): uncentered contribution cross-products
  sigma0 <- crossprod(H0) / n0
  sigma1 <- crossprod(H1) / n1
  
  # Equation (32), complete-data case
  V_theta <- N^2 * (n0 * sigma0 + n1 * sigma1) /
    (n0^2 * n1^2)
  
  V_theta <- (V_theta + t(V_theta)) / 2
  
  names(theta) <- time_cols
  dimnames(V_theta) <- list(time_cols, time_cols)
  
  # SAS homogeneity contrasts:
  # theta_1 - theta_2, theta_2 - theta_3, ...
  L <- matrix(0, nrow = K - 1L, ncol = K)
  
  for (k in seq_len(K - 1L)) {
    L[k, k] <- 1
    L[k, k + 1L] <- -1
  }
  
  delta <- drop(L %*% theta)
  V_delta <- L %*% V_theta %*% t(L)
  V_delta <- (V_delta + t(V_delta)) / 2
  
  test <- tryCatch({
    
    eigenvalues <- eigen(
      V_delta,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    
    tolerance <- 1e-10 * max(abs(eigenvalues))
    
    if (any(!is.finite(eigenvalues)) ||
        min(eigenvalues) <= tolerance) {
      stop("Interaction covariance is singular or not positive definite.")
    }
    
    R <- chol(V_delta)
    standardized <- forwardsolve(t(R), delta)
    statistic <- sum(standardized^2)
    
    data.frame(
      procedure = "Lachin_Homogeneity",
      statistic = statistic,
      df = K - 1L,
      p_value = pchisq(
        statistic,
        df = K - 1L,
        lower.tail = FALSE
      ),
      status = "ok",
      error = NA_character_
    )
    
  }, error = function(e) {
    
    data.frame(
      procedure = "Lachin_Homogeneity",
      statistic = NA_real_,
      df = K - 1L,
      p_value = NA_real_,
      status = "failed",
      error = conditionMessage(e)
    )
  })
  
  test$reject <- test$p_value < alpha
  
  list(
    theta = theta,
    PI = (theta + 1) / 2,
    vcov_theta = V_theta,
    contrasts = L,
    contrast_estimates = delta,
    vcov_contrasts = V_delta,
    global = test
  )
}



simulate_PO_PIM <- function(
    n_per_group = 30,
    n_visits = 5,
    cutpoints = c(-1.5, -0.5, 0.5, 1.5),
    treatment_PI = 0.65,
    temporal_PI = 0.60,
    time_shift_grid = seq(0.05, 2, length.out = 100),
    seed = NULL,
    tolerance = 1e-7
) {
  if (!requireNamespace("lpSolve", quietly = TRUE))
    stop("Install lpSolve first.")
  
  stopifnot(
    n_per_group >= 1,
    n_visits >= 2,
    n_visits == as.integer(n_visits),
    length(cutpoints) >= 1,
    all(is.finite(cutpoints)),
    all(diff(cutpoints) > 0),
    treatment_PI > 0, treatment_PI < 1,
    temporal_PI > 0, temporal_PI < 1,
    length(time_shift_grid) > 0,
    all(is.finite(time_shift_grid))
  )
  
  if (!is.null(seed)) set.seed(seed)
  
  visits <- seq_len(n_visits)
  time <- visits - 1
  K <- length(cutpoints) + 1L
  
  po_probs <- function(eta) {
    diff(c(0, plogis(cutpoints - eta), 1))
  }
  
  PI <- function(p0, p1) {
    below <- c(0, head(cumsum(p0), -1))
    sum(p1 * (below + 0.5 * p0))
  }
  
  calibrate_delta <- function(eta0) {
    p0 <- po_probs(eta0)
    
    if (treatment_PI <= 0.5 * p0[1] ||
        treatment_PI >= 1 - 0.5 * p0[K])
      return(NA_real_)
    
    f <- function(delta) {
      PI(p0, po_probs(eta0 + delta)) - treatment_PI
    }
    
    bound <- 1
    while (f(-bound) > 0 || f(bound) < 0) {
      bound <- 2 * bound
      if (bound > 1e6) return(NA_real_)
    }
    
    uniroot(f, c(-bound, bound), tol = 1e-12)$root
  }
  
  # All K^n_visits possible trajectories
  trajectory_grid <- rep(list(seq_len(K)), n_visits)
  names(trajectory_grid) <- paste0("visit", visits)
  
  trajectories <- as.matrix(do.call(
    expand.grid,
    c(trajectory_grid, list(
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ))
  ))
  M <- nrow(trajectories)
  
  pairs <- t(combn(visits, 2))
  betaT <- qlogis(temporal_PI)
  
  temporal_targets <- plogis(
    betaT * (pairs[, 2] - pairs[, 1])
  )
  
  # Only K-1 categories per visit are needed:
  # the final category follows from normalization.
  marginal_rows <- do.call(rbind, lapply(visits, function(t) {
    do.call(rbind, lapply(seq_len(K - 1L), function(k) {
      as.numeric(trajectories[, t] == k)
    }))
  }))
  
  temporal_rows <- t(vapply(seq_len(nrow(pairs)), function(j) {
    earlier <- trajectories[, pairs[j, 1]]
    later <- trajectories[, pairs[j, 2]]
    
    as.numeric(later > earlier) +
      0.5 * as.numeric(later == earlier)
  }, numeric(M)))
  
  Aeq <- rbind(rep(1, M), marginal_rows, temporal_rows)
  
  # Write q_j = z + r_j, with z >= 0 and r_j >= 0.
  # Maximizing z maximizes the minimum trajectory probability.
  # This avoids constructing an M-by-M identity matrix.
  constraint_matrix <- cbind(Aeq, rowSums(Aeq))
  
  solve_joint <- function(margins) {
    rhs <- c(
      1,
      as.vector(margins[seq_len(K - 1L), , drop = FALSE]),
      temporal_targets
    )
    
    fit <- lpSolve::lp(
      direction = "max",
      objective.in = c(rep(0, M), 1),
      const.mat = constraint_matrix,
      const.dir = rep("=", nrow(Aeq)),
      const.rhs = rhs
    )
    
    if (fit$status != 0) return(NULL)
    
    z <- fit$solution[M + 1L]
    probabilities <- pmax(fit$solution[seq_len(M)] + z, 0)
    probabilities <- probabilities / sum(probabilities)
    
    if (max(abs(drop(Aeq %*% probabilities) - rhs)) > tolerance)
      return(NULL)
    
    probabilities
  }
  
  chosen <- NULL
  
  for (time_shift in time_shift_grid) {
    eta0 <- time_shift * time
    delta <- vapply(eta0, calibrate_delta, numeric(1))
    
    if (anyNA(delta)) next
    
    margins0 <- vapply(eta0, po_probs, numeric(K))
    margins1 <- vapply(eta0 + delta, po_probs, numeric(K))
    
    joint0 <- solve_joint(margins0)
    if (is.null(joint0)) next
    
    joint1 <- solve_joint(margins1)
    if (is.null(joint1)) next
    
    chosen <- list(
      time_shift = time_shift,
      delta = delta,
      margins0 = margins0,
      margins1 = margins1,
      joint0 = joint0,
      joint1 = joint1
    )
    break
  }
  
  if (is.null(chosen)) {
    stop(paste(
      "No feasible joint distribution found on this grid.",
      "Change the grid, cutpoints, or target PIs.",
      "This is not proof of infeasibility for every PO construction."
    ))
  }
  
  sample_trajectories <- function(probabilities) {
    index <- sample.int(
      M, n_per_group,
      replace = TRUE,
      prob = probabilities
    )
    trajectories[index, , drop = FALSE]
  }
  
  Y <- rbind(
    sample_trajectories(chosen$joint0),
    sample_trajectories(chosen$joint1)
  )
  
  treatment <- rep(0:1, each = n_per_group)
  
  dat <- data.frame(
    id = rep(seq_len(nrow(Y)), each = n_visits),
    treatment = rep(treatment, each = n_visits),
    visit = rep(visits, times = nrow(Y)),
    outcome = as.vector(t(Y))
  )
  
  treatment_check <- data.frame(
    visit = visits,
    target = treatment_PI,
    achieved = vapply(visits, function(t) {
      PI(chosen$margins0[, t], chosen$margins1[, t])
    }, numeric(1))
  )
  
  temporal_check <- data.frame(
    earlier_visit = pairs[, 1],
    later_visit = pairs[, 2],
    target = temporal_targets,
    control = drop(temporal_rows %*% chosen$joint0),
    treatment = drop(temporal_rows %*% chosen$joint1)
  )
  
  # Check all population constraints before returning the dataset
  errors <- c(
    treatment_check$achieved - treatment_check$target,
    temporal_check$control - temporal_check$target,
    temporal_check$treatment - temporal_check$target
  )
  
  if (max(abs(errors)) > tolerance)
    stop("Population PI verification failed.")
  
  list(
    data = dat,
    truth = list(
      beta0 = 0,
      betaAt = rep(qlogis(treatment_PI), n_visits),
      betaT0 = betaT,
      betaT1 = betaT
    ),
    treatment_check = treatment_check,
    temporal_check = temporal_check,
    PO_parameters = data.frame(
      visit = visits,
      control_location = chosen$time_shift * time,
      treatment_shift = chosen$delta
    ),
    trajectories = trajectories,
    joint_probabilities = data.frame(
      control = chosen$joint0,
      treatment = chosen$joint1
    )
  )
}



test_contrasts <- function(b, V, L, nu, alpha = 0.05) {
  
  b <- as.numeric(b)
  V <- as.matrix(V)
  L <- as.matrix(L)
  
  stopifnot(
    ncol(L) == length(b),
    all(dim(V) == length(b)),
    is.finite(nu), nu > 0,
    all(is.finite(b)),
    all(is.finite(V)),
    all(is.finite(L))
  )
  
  V <- (V + t(V)) / 2
  
  # Elementary contrasts: retain every row for Holm
  estimates <- drop(L %*% b)
  contrast_V <- L %*% V %*% t(L)
  variances <- diag(contrast_V)
  
  if (any(variances <= 0))
    stop("Non-positive contrast variance.")
  
  SE <- sqrt(variances)
  statistics <- estimates / SE
  
  p_z <- 2 * pnorm(-abs(statistics))
  p_t <- 2 * pt(-abs(statistics), df = nu)
  
  p_Holm_z <- p.adjust(p_z, method = "holm")
  p_Holm_t <- p.adjust(p_t, method = "holm")
  
  # Independent rows for the omnibus test
  q <- qr(t(L))
  r <- q$rank
  
  if (r == 0)
    stop("Contrast matrix has rank zero.")
  
  L_global <- L[q$pivot[seq_len(r)], , drop = FALSE]
  
  delta <- drop(L_global %*% b)
  V_global <- L_global %*% V %*% t(L_global)
  V_global <- (V_global + t(V_global)) / 2
  
  # Reject invalid covariance rather than silently repairing it
  global_result <- tryCatch({
    
    R <- chol(V_global)
    standardized <- forwardsolve(t(R), delta)
    W <- sum(standardized^2)
    
    list(
      W = W,
      p_chisq = pchisq(W, df = r, lower.tail = FALSE),
      p_F = pf(W / r, df1 = r, df2 = nu, lower.tail = FALSE),
      status = "ok",
      error = NA_character_
    )
    
  }, error = function(e) {
    
    list(
      W = NA_real_,
      p_chisq = NA_real_,
      p_F = NA_real_,
      status = "failed",
      error = conditionMessage(e)
    )
  })
  
  W <- global_result$W
  p_chisq <- global_result$p_chisq
  p_F <- global_result$p_F
  
  list(
    global = data.frame(
      procedure = c(
        "Asymptotic_Global",
        "Asymptotic_Holm",
        "Finite_Global",
        "Finite_Holm"
      ),
      statistic = c(W, NA_real_, W / r, NA_real_),
      df1 = c(r, NA_real_, r, NA_real_),
      df2 = c(NA_real_, NA_real_, nu, nu),
      p_global = c(
        p_chisq,
        min(p_Holm_z),
        p_F,
        min(p_Holm_t)
      ),
      reject = c(
        p_chisq,
        min(p_Holm_z),
        p_F,
        min(p_Holm_t)
      ) < alpha,
      
      status = c(
        global_result$status,
        "ok",
        global_result$status,
        "ok"
      ),
      
      error = c(
        global_result$error,
        NA_character_,
        global_result$error,
        NA_character_
      )
    ),
    
    elementary = data.frame(
      contrast = seq_len(nrow(L)),
      estimate = estimates,
      SE = SE,
      z_or_t = statistics,
      df_t = nu,
      p_z = p_z,
      p_Holm_z = p_Holm_z,
      p_t = p_t,
      p_Holm_t = p_Holm_t
    )
  )
}

sandwich.estimator_clustered_mod1_mod2 <- function(U, U.diff, 
                                                   g1, g2, 
                                                   shared.factor = 1, switched.factor = 1,
                                                   self.factor = 1) {
  # Note: this estimator is greatly optimized, based on T Lumley's code!
  
  # We transform the original g1 and g2 from the comparison of two groups to indicators corresponding to
  # the overarching subject ids. This leads to additional terms that are counted twice in utu (as compared to
  # only the diagonal elements in the original code). This is accounted for using the subtract_df below. All shared/switched/selffactors are fixed to 1
  
  g1 = individuals_var1
  g2 = individuals_var2
  
  
  use_n1 = length(unique(g1))
  
  tim = proc.time()[3]
  
  usum1.tmp <- rowsum(U, g1, reorder = FALSE)
  usum2.tmp <- rowsum(U, g2, reorder = FALSE)
  
  ng <- length(union(g1, g2))
  
  usum1 <- matrix(nrow = ng, ncol = ncol(usum1.tmp), 0)
  usum2 <- matrix(nrow = ng, ncol = ncol(usum2.tmp), 0)
  usum1[unique(g1), ] <- usum1.tmp
  usum2[unique(g2), ] <- usum2.tmp
  
  
  utushared <- crossprod(usum1) + crossprod(usum2) 
  utuswitched <- crossprod(usum1, usum2) + crossprod(usum2, usum1)
  
  
  
  subtract_df = data.frame(g1,g2)
  subtract_df = subtract_df %>%group_by(g1, g2)%>% mutate(ID = cur_group_id()) 
  if(length(unique(subtract_df$ID)) == length(g1))
  {uDiag <- crossprod(U)  
  utu <- shared.factor * utushared + 
    (switched.factor) * utuswitched +
    (self.factor - 2 * shared.factor) * (uDiag)
  }
  
  
  if(length(unique(subtract_df$ID)) < length(g1))
  { 
    usum3.tmp <- rowsum(U, subtract_df$ID, reorder = FALSE)
    subtract = crossprod(usum3.tmp[1:(use_n1/2*use_n1/2),]) 
    utu <- utushared - subtract 
    utu[ncol(U),ncol(U)] = utu[ncol(U),ncol(U)]+utuswitched[ncol(U),ncol(U)]
    subtract = crossprod(usum3.tmp[(use_n1/2*use_n1/2+1):(use_n1+use_n1/2*use_n1/2),])  
    utu <- utu - subtract
    
  }
  
  
  
  mx <- ng + 1
  
  
  # Skip this last part if U.diff is identity... (note: NULL here will mean identity)
  if (is.null(U.diff)) return(utu)
  
  U.diff.inv<-solve(U.diff)
  return(U.diff.inv%*%utu%*%U.diff.inv)
}



## Original sandwich estimator of the PIM (only sums U over shared indices)
sandwich.estimator<-function(U, U.diff, 
                             g1, g2, 
                             shared.factor=1, switched.factor=1,
                             self.factor=1)
{
  
  #note: this estimator is greatly optimized, based on T Lumley's code!
  
  if(any(fullyequal <- g1==g2 ))
  {
    #the contribution of these terms to variance and covariance is zero anyway, because
    #I_{ii} = 0.5 always (ie: constant), so we drop them out!
    U<-U[!fullyequal,,drop=FALSE]
    g1<-g1[!fullyequal]
    g2<-g2[!fullyequal]
  }
  
  usum1.tmp <- rowsum(U,g1,reorder=FALSE)
  usum2.tmp <- rowsum(U,g2,reorder=FALSE)
  
  ng <- length(union(g1,g2))
  
  usum1  <- matrix(nrow = ng, ncol = ncol(usum1.tmp),0)
  usum2  <- matrix(nrow = ng, ncol = ncol(usum2.tmp),0)
  usum1[unique(g1),] <- usum1.tmp
  usum2[unique(g2),] <- usum2.tmp
  
  utushared <- crossprod(usum1) + crossprod(usum2) 
  utuswitched <- crossprod(usum1,usum2) + crossprod(usum2,usum1)
  uDiag<-crossprod(U) #Is counted twice as shared.factor, but needs to be counted as self.factor
  
  
  utu<-shared.factor*utushared  + 
    (switched.factor)*utuswitched +
    (self.factor-2*shared.factor)*(uDiag)
  
  #if the inverses occur (ij, ji), they are counted doubly as switched!!
  #However: they should be counted as - self
  
  mx<-ng+1
  
  uids<-g1*mx+g2
  invuids<-g2*mx+g1
  invrowperrow<-match(uids, invuids)
  rowsWithInv<-(!is.na(invrowperrow)) 
  if(any(rowsWithInv))
  {
    rowsWithInv<-which(rowsWithInv)
    tmp<-sapply(rowsWithInv, function(ri){
      res<-U[ri,] %*% t(U[invrowperrow[ri],]) 
      return(res)
    })
    if(is.null(dim(tmp)))
    {
      if(length(rowsWithInv)==1)
      {
        tmp<-matrix(tmp, ncol=1)
      }
      else
      {
        tmp<-matrix(tmp, nrow=1)
      }
    }
    
    uijji<-matrix(rowSums(tmp), ncol=ncol(usum1.tmp))
    
    utu<-utu-(2 *  switched.factor + self.factor)*uijji
  }
  
  
  
  individuals_var2 = g2
  individuals_var1 = g1
  #
  # utu_test = matrix(0,nrow=ncol(U),ncol=ncol(U))
  # phi = ((outer(individuals_var1, individuals_var1, FUN = "==")+outer(individuals_var2, individuals_var2, FUN = "==")+outer(individuals_var1, individuals_var2, FUN = "==")+outer(individuals_var2, individuals_var1, FUN = "=="))>0)*1
  # for(j in 1:nrow(U)){
  #   for(k in 1:nrow(U)){
  #     utu_test = utu_test+ (phi[j,k]*U[j,]%*%t(U[k,]))
  #   }
  # }
  # rownames(utu)=rownames(utu_test)=colnames(utu_test) = colnames(utu) = c()
  # 
  
  #skip this last part if U.diff is identity... (note: NULL here will mean identity)
  if(is.null(U.diff)) return (utu)
  
  U.diff.inv<-solve(U.diff)
  
  return(U.diff.inv%*%utu%*%U.diff.inv)
}

# Proportional Odds model
nr_of_times <- 4L

treatment_PI <- 0.65
temporal_PI <- 0.60

po_template <- simulate_PO_PIM(
  n_per_group = 1,
  n_visits = nr_of_times,
  cutpoints = c(-3, -1, 1, 3),
  treatment_PI = treatment_PI,
  temporal_PI = temporal_PI,
  time_shift_grid = seq(0.05, 1, by = 0.01),
  seed = 123
)

# Check the calibrated population PIs
po_template$treatment_check
po_template$temporal_check

trajectories <- po_template$trajectories

prob_control <- po_template$joint_probabilities$control
prob_treatment <- po_template$joint_probabilities$treatment

# Order matches your model:
# temporal treatment, temporal control, treatment at visits 1–4
truth <- c(
  po_template$truth$betaT1,
  po_template$truth$betaT0,
  po_template$truth$betaAt
)



L_trt = matrix(c(c(0,0,c(3/4,rep(-1/4,3))),
                 c(0,0,c(-1/4,3/4,rep(-1/4,2))),
                 c(0,0,c(rep(-1/4,2),3/4,-1/4)),
                 c(0,0,c(rep(-1/4,3),3/4))),nrow=4,byrow = TRUE)
L_time = rbind(c(1, -1,rep(0,4)))


sample_sizes <- c(16, 18, 20, 24, 30, 40, 50, 60)
nr_of_iters <- 10000

n_datasets <- length(sample_sizes ) * nr_of_iters

test_store <- vector("list", 8L * n_datasets)
contrast_store <- vector("list", 8L * n_datasets)
model_store <- vector("list", 4L * n_datasets)
benchmark_store <- vector("list", n_datasets)
lachin_store <- vector("list", n_datasets)


test_index <- 0L
store_index <- 0L
model_index <- 0L
benchmark_index <- 0L
lachin_index <- 0L

for (N in sample_sizes) {
  
  nr_of_subjects <- N
  id = rep(1:nr_of_subjects,each=nr_of_times)
  time = rep(1:nr_of_times,nr_of_subjects)
  group = rep(rep(c(0,1),each=nr_of_subjects/2),each=nr_of_times)
  
  design  = data.frame(id,time,group)
  
  for(l in 1:nr_of_iters){
    if(l%%100==0){cat("Iteration: ",l,"\n")}
    set.seed(l)
    
    index0 <- sample.int(
      nrow(po_template$trajectories),
      N / 2,
      replace = TRUE,
      prob = po_template$joint_probabilities$control
    )
    
    index1 <- sample.int(
      nrow(po_template$trajectories),
      N / 2,
      replace = TRUE,
      prob = po_template$joint_probabilities$treatment
    )
    
    design$Y <- as.vector(t(rbind(
      po_template$trajectories[index0, , drop = FALSE],
      po_template$trajectories[index1, , drop = FALSE]
    )))
    
    # Fit Lachin test
    
    lachin_rows <- tryCatch({
      
      wide <- tidyr::pivot_wider(
        design[, c("id", "group", "time", "Y")],
        id_cols = c(id, group),
        names_from = time,
        values_from = Y,
        names_prefix = "visit"
      )
      
      fit_lachin <- lachin_constancy_test(
        data = as.data.frame(wide),
        group_col = "group",
        time_cols = paste0("visit", seq_len(nr_of_times))
      )
      
      fit_lachin$global
      
    }, error = function(e) {
      
      data.frame(
        procedure = "Lachin_Homogeneity",
        statistic = NA_real_,
        df = nr_of_times - 1L,
        p_value = NA_real_,
        status = "failed",
        error = conditionMessage(e),
        reject = NA
      )
    })
    
    lachin_index <- lachin_index + 1L
    
    lachin_store[[lachin_index]] <- data.frame(
      N = N,
      run = l,
      model = "Lachin",
      hypothesis = "treatment",
      lachin_rows,
      row.names = NULL
    )
    
    # Fit benchmark nparLD
    
    benchmark_result <- tryCatch({
      
      dat_nparLD <- design
      dat_nparLD$time <- factor(dat_nparLD$time)
      dat_nparLD$group <- factor(dat_nparLD$group)
      
      mod_nparLD <- nparLD(Y ~ time * group, data = design,
                           subject = "id", description = FALSE)
      
      # Extract the interaction row explicitly
      interaction_p <- function(tab) {
        tab <- as.matrix(tab)
        idx <- grep(":", rownames(tab), fixed = TRUE)
        
        if (length(idx) != 1L) {
          stop("Could not identify a unique interaction row.")
        }
        
        # Third column is the p-value in these test tables
        as.numeric(tab[idx, 3L])
      }
      
      p <- c(
        interaction_p(mod_nparLD$Wald.test),
        interaction_p(mod_nparLD$ANOVA.test)
      )
      
      valid <- is.finite(p) & p >= 0 & p <= 1
      
      data.frame(
        procedure = c("nparLD_WTS", "nparLD_ATS"),
        p_value = p,
        reject = ifelse(valid, p < 0.05, NA),
        status = ifelse(valid, "ok", "invalid p-value"),
        error = NA_character_
      )
      
    }, error = function(e) {
      
      data.frame(
        procedure = c("nparLD_WTS", "nparLD_ATS"),
        p_value = NA_real_,
        reject = NA,
        status = "failed",
        error = conditionMessage(e)
      )
    })
    
    benchmark_index <- benchmark_index + 1L
    
    benchmark_store[[benchmark_index]] <- data.frame(
      N = N,
      run = l,
      hypothesis = "RTE interaction",
      benchmark_result,
      row.names = NULL
    )
    
    ## Model poset, defined via compare.
    
    id.fac1 <- which((design[,"group"] == 0)&(design[,"time"] == 1))
    id.nonfac1 <- which((design[,"group"] == 1)&(design[,"time"] == 1))
    compare1 <- expand.grid(id.fac1,id.nonfac1)
    
    id.fac2 <- which((design[,"group"] == 0)&(design[,"time"] == 2))
    id.nonfac2 <- which((design[,"group"] == 1)&(design[,"time"] == 2))
    compare2 <- expand.grid(id.fac2,id.nonfac2)
    
    id.fac3 <- which((design[,"group"] == 0)&(design[,"time"] == 3))
    id.nonfac3 <- which((design[,"group"] == 1)&(design[,"time"] == 3))
    compare3 <- expand.grid(id.fac3,id.nonfac3)
    
    id.fac4 <- which((design[,"group"] == 0)&(design[,"time"] == 4))
    id.nonfac4 <- which((design[,"group"] == 1)&(design[,"time"] == 4))
    compare4 <- expand.grid(id.fac4,id.nonfac4)
    
    compare_between = rbind(compare1,compare2,compare3,compare4)
    compare_between = compare_between[order(compare_between[,"Var2"]),]
    
    start_timepoints <- c()
    later_timepoints <- c()
    k=1
    for (i in 1:(nrow(design) - 1)) {
      if(i%%4 != 0){
        start <- i           # Current time point
        later <- (i+1) : (4*k)  # Later time points
        # Store the results in vectors
        start_timepoints <- c(start_timepoints, rep(start, length(later)))
        later_timepoints <- c(later_timepoints, later)}
      if(i%%4 == 0){
        k=k+1}
    }
    within1 = start_timepoints
    within2 = later_timepoints
    
    compare_within = cbind(within1,within2)
    colnames(compare_within) = c("Var1","Var2")
    compare_within = compare_within[order(compare_within[,2]),]
    compare = rbind(compare_between,compare_within)
    
    
    individuals_var1=design[compare$Var1,"id"]
    individuals_var2=design[compare$Var2,"id"]
    
    # Fit original PIM
    assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")
    mod2_orig  <- pim(Y~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) +I((R(group)-L(group))*(R(time)==4)*(L(time)==4)),data=design,compare=compare,link="logit",start = truth)
    
    
    V_PIM = vcov(mod2_orig)
    b_PIM = coef(mod2_orig)
    
    # Fit sandwich adjusted PIM
    assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")
    mod2_adj  <- pim(Y~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) +I((R(group)-L(group))*(R(time)==4)*(L(time)==4)),data=design,compare=compare,link="logit",start=truth)
    
    V_PIM_adj = vcov(mod2_adj)
    b_PIM_adj = coef(mod2_adj)
    
    # Fit 3 GEE models with PGEE for small sample inference
    
    # Prepare data
    
    C1 = individuals_var1
    C2 = individuals_var2
    
    x = model.matrix(mod2_adj)
    y = response(mod2_adj)
    
    dat_GEE = data.frame(y,x)
    
    names(dat_GEE) = c("y",c(paste0("x",1:6)))
    dat_GEE$C1=C1
    dat_GEE$C2=C2
    dat_GEE$C3 <- paste(dat_GEE$C1, dat_GEE$C2,sep="_")
    
    
    # Fit models
    
    dat_GEE=dat_GEE[order(dat_GEE$C1),]
    mod1 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C1, corstr = "independence",beta.method="GEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C2),]
    mod2 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C2, corstr = "independence",beta.method="GEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C3),]
    mod3 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C3, corstr = "independence",beta.method="GEE",SE.method = "FW")
    
    
    V_GEE=mod1$covb+mod2$covb-mod3$covb
    b_GEE = coef(mod1)
    
    # Fit 3 GEE models with PGEE for small sample inference
    
    # Prepare data
    
    C1 = individuals_var1
    C2 = individuals_var2
    
    x = model.matrix(mod2_adj)
    y = response(mod2_adj)
    
    dat_GEE = data.frame(y,x)
    
    names(dat_GEE) = c("y",c(paste0("x",1:6)))
    dat_GEE$C1=C1
    dat_GEE$C2=C2
    dat_GEE$C3 <- paste(dat_GEE$C1, dat_GEE$C2,sep="_")
    
    
    # Fit models
    
    dat_GEE=dat_GEE[order(dat_GEE$C1),]
    mod1 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C1, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C2),]
    mod2 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C2, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C3),]
    mod3 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C3, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    
    b2 <- rbind(
      C1 = coef(mod1),
      C2 = coef(mod2),
      C3 = coef(mod3)
    )
    
    # Largest absolute difference across any two fits
    max_difference <- max(apply(b2, 2, function(x) diff(range(x))))
    
    # Equal within numerical tolerance?
    if(max_difference > 1e-5){print("parameter estimates not equal")}
    
    
    
    V_PGEE=  mod1$covb+mod2$covb-mod3$covb
    b_PGEE = coef(mod1)
    
    # 
    # ### Interaction tests based on the different models
    
    
    models <- list(
      PIM = list(
        b = b_PIM,
        V = V_PIM
      ),
      PIM_adj = list(
        b = b_PIM_adj,
        V = V_PIM_adj
      ),
      GEE = list(
        b = b_GEE,
        V = V_GEE
      ),
      PGEE = list(
        b = b_PGEE,
        V = V_PGEE
      )
    )
    
    
    for (model_name in names(models)) {
      
      model_index <- model_index + 1L
      
      model_store[[model_index]] <- list(
        N = N,
        run = l,
        model = model_name,
        coefficients = models[[model_name]]$b,
        vcov = as.matrix(models[[model_name]]$V)
      )
    }
    
    
    contrast_matrices <- list(
      treatment = L_trt,
      time = L_time
    )
    
    # Candidate approximations; specify these deliberately
    nu_trt <- nr_of_subjects - 1
    nu_time <- nr_of_subjects - 1
    
    denominator_df <- c(
      treatment = nu_trt,
      time = nu_time
    )
    
    
    global_rows <- list()
    elementary_rows <- list()
    k <- 0L
    
    for (model_name in names(models)) {
      for (hypothesis in names(contrast_matrices)) {
        
        b <- models[[model_name]]$b
        V <- models[[model_name]]$V
        L <- contrast_matrices[[hypothesis]]
        
        # Use one consistent definition of degrees of freedom
        nu <- if (hypothesis == "treatment") nu_trt else nu_time
        
        # Store contrast estimates independently of global-test success
        contrast_result <- tryCatch({
          
          estimate <- drop(L %*% b)
          contrast_V <- L %*% V %*% t(L)
          estimated_variance <- diag(contrast_V)
          
          data.frame(
            contrast = seq_len(nrow(L)),
            estimate = estimate,
            variance = estimated_variance,
            status = ifelse(
              is.finite(estimate) &
                is.finite(estimated_variance) &
                estimated_variance > 0,
              "ok",
              "invalid estimate or variance"
            ),
            error = NA_character_
          )
          
        }, error = function(e) {
          
          data.frame(
            contrast = seq_len(nrow(L)),
            estimate = NA_real_,
            variance = NA_real_,
            status = "failed",
            error = conditionMessage(e)
          )
        })
        
        # Run the tests once; catch errors and continue
        result <- tryCatch({
          
          res <- test_contrasts(
            b = b,
            V = V,
            L = L,
            nu = nu,
            alpha = 0.05
          )
          
          if (!"status" %in% names(res$global)) {
            res$global$status <- ifelse(
              is.na(res$global$reject), "failed", "ok"
            )
          }
          
          if (!"error" %in% names(res$global)) {
            res$global$error <- NA_character_
          }
          res
          
        }, error = function(e) {
          
          list(
            global = data.frame(
              procedure = c(
                "Asymptotic_Global",
                "Asymptotic_Holm",
                "Finite_Global",
                "Finite_Holm"
              ),
              statistic = NA_real_,
              df1 = NA_real_,
              df2 = NA_real_,
              p_global = NA_real_,
              reject = NA,
              status = "failed",
              error = conditionMessage(e)
            ),
            
            elementary = data.frame(
              contrast = seq_len(nrow(L)),
              estimate = contrast_result$estimate,
              SE = sqrt(ifelse(
                contrast_result$variance > 0,
                contrast_result$variance,
                NA_real_
              )),
              z_or_t = NA_real_,
              df_t = nu,
              p_z = NA_real_,
              p_Holm_z = NA_real_,
              p_t = NA_real_,
              p_Holm_t = NA_real_
            )
          )
        })
        
        # Store test results across all simulation runs
        test_index <- test_index + 1L
        
        test_store[[test_index]] <- data.frame(
          N = N,
          run = l,
          model = model_name,
          hypothesis = hypothesis,
          result$global,
          row.names = NULL
        )
        
        # Store contrast estimates and variances across all runs
        store_index <- store_index + 1L
        
        contrast_store[[store_index]] <- data.frame(
          N = N,
          run = l,
          model = model_name,
          hypothesis = hypothesis,
          contrast_result,
          row.names = NULL
        )
        
        # Results for the current simulation run
        k <- k + 1L
        
        global_rows[[k]] <- test_store[[test_index]]
        
        elementary_rows[[k]] <- data.frame(
          N = N,
          run = l,
          model = model_name,
          hypothesis = hypothesis,
          result$elementary,
          row.names = NULL
        )
      }
    }
    
    global_results <- dplyr::bind_rows(global_rows)
    elementary_results <- dplyr::bind_rows(elementary_rows)
    
  }
}

test_runs <- dplyr::bind_rows(test_store)
contrast_runs <- dplyr::bind_rows(contrast_store)

table(test_runs$N)
table(contrast_runs$N)

model_levels <- c("GEE", "PGEE", "PIM", "PIM_adj")

procedure_levels <- c(
  "Asymptotic_Global", "Asymptotic_Holm",
  "Finite_Global", "Finite_Holm"
)

# ------------------------------------------------------------
# 1. Variance comparison
# ------------------------------------------------------------

contrast_runs <- bind_rows(contrast_store) |>
  mutate(
    valid = is.finite(estimate) &
      is.finite(variance) &
      variance > 0
  )

# Detect accidental duplicate storage
stopifnot(
  !anyDuplicated(
    contrast_runs[c("N", "run", "model", "hypothesis", "contrast")]
  )
)

variance_comparison <- contrast_runs |>
  group_by(N, model, hypothesis, contrast) |>
  summarise(
    n_total = n(),
    n_valid = sum(valid),
    n_failed = sum(!valid),
    
    mean_estimate = if (any(valid)) {
      mean(estimate[valid])
    } else {
      NA_real_
    },
    
    empirical_variance = if (sum(valid) >= 2L) {
      var(estimate[valid])
    } else {
      NA_real_
    },
    
    average_estimated_variance = if (any(valid)) {
      mean(variance[valid])
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) |>
  mutate(
    failure_rate = n_failed / n_total,
    empirical_SD = sqrt(empirical_variance),
    estimated_SD = sqrt(average_estimated_variance),
    variance_ratio =
      empirical_variance / average_estimated_variance
  )

print(variance_comparison, n = Inf)

# ------------------------------------------------------------
# 2. Rejection rates
# ------------------------------------------------------------

test_runs <- bind_rows(test_store)

stopifnot(
  !anyDuplicated(
    test_runs[c("N", "run", "model", "hypothesis", "procedure")]
  )
)

rejection_summary <- test_runs |>
  group_by(N, model, hypothesis, procedure) |>
  summarise(
    n_total = n(),
    n_valid = sum(!is.na(reject)),
    n_failed = sum(is.na(reject)),
    n_rejected = sum(reject, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    failure_rate = n_failed / n_total,
    rejection_rate = n_rejected / na_if(n_valid, 0L),
    
    MCSE = sqrt(
      rejection_rate * (1 - rejection_rate) /
        na_if(n_valid, 0L)
    )
  )

# Wilson intervals: also work when there are zero rejections
z <- qnorm(0.975)

rejection_summary <- rejection_summary |>
  mutate(
    denominator = na_if(n_valid, 0L),
    interval_center =
      (rejection_rate + z^2 / (2 * denominator)) /
      (1 + z^2 / denominator),
    
    interval_halfwidth =
      z * sqrt(
        rejection_rate * (1 - rejection_rate) / denominator +
          z^2 / (4 * denominator^2)
      ) / (1 + z^2 / denominator),
    
    lower = pmax(0, interval_center - interval_halfwidth),
    upper = pmin(1, interval_center + interval_halfwidth)
  ) |>
  dplyr::select(-denominator, -interval_center, -interval_halfwidth)

print(rejection_summary, n = Inf)

# ------------------------------------------------------------
# 3. Variance-ratio plot
# ------------------------------------------------------------

variance_plot_data <- variance_comparison |>
  mutate(
    model = factor(model, levels = model_levels),
    hypothesis = factor(
      hypothesis,
      levels = c("time", "treatment"),
      labels = c("Equal temporal effects", "Treatment-effect constancy")
    ),
    contrast_label = paste("Contrast", contrast)
  )

p_variance <- ggplot(
  variance_plot_data,
  aes(
    x = N,
    y = variance_ratio,
    shape = model,
    linetype = model,
    group = model
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line() +
  geom_point(size = 2.5) +
  facet_wrap(vars(hypothesis, contrast_label), ncol = 2) +
  scale_shape_manual(
    values = c(GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2)
  ) +
  scale_linetype_manual(
    values = c(
      GEE = "solid",
      PGEE = "dashed",
      PIM = "dotted",
      PIM_adj = "dotdash"
    )
  ) +
  scale_x_continuous(breaks = sample_sizes) +
  labs(
    x = "Number of subjects",
    y = "Empirical variance / mean estimated variance",
    shape = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_variance)


p_bias <- ggplot(
  variance_plot_data,
  aes(
    x = N,
    y = mean_estimate,
    shape = model,
    linetype = model,
    group = model
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line() +
  geom_point(size = 2.5) +
  facet_wrap(vars(hypothesis, contrast_label), ncol = 2) +
  scale_shape_manual(
    values = c(GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2)
  ) +
  scale_linetype_manual(
    values = c(
      GEE = "solid",
      PGEE = "dashed",
      PIM = "dotted",
      PIM_adj = "dotdash"
    )
  ) +
  scale_x_continuous(breaks = sample_sizes) +
  labs(
    x = "Number of subjects",
    y = "Empirical mean of the contrasts",
    shape = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_bias)


# ------------------------------------------------------------
# 4. Rejection-rate plot
# ------------------------------------------------------------

rejection_plot_data <- rejection_summary |>
  mutate(
    model = factor(model, levels = model_levels),
    procedure = factor(procedure, levels = procedure_levels)
  )

procedure_labels <- c(
  Asymptotic_Global = "Global: chi-square",
  Asymptotic_Holm   = "Holm: z",
  Finite_Global    = "Global: F",
  Finite_Holm      = "Holm: t"
)

pd <- position_dodge(width = 0.8)

p_rejection <- ggplot(
  rejection_plot_data,
  aes(
    x = N,
    y = rejection_rate,
    shape = procedure,
    linetype = procedure,
    group = procedure
  )
) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = pd,
    width = 0.5
  ) +
  facet_grid(hypothesis ~ model) +
  scale_shape_manual(
    values = setNames(c(16, 17, 1, 2), procedure_levels),
    labels = procedure_labels
  ) +
  scale_linetype_manual(
    values = setNames(
      c("solid", "dashed", "dotted", "dotdash"),
      procedure_levels
    ),
    labels = procedure_labels
  ) +
  scale_x_continuous(breaks = sample_sizes) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    x = "Number of subjects",
    y = "Rejection rate",
    shape = "Procedure",
    linetype = "Procedure",
    caption = "Error bars: 95% Wilson Monte Carlo intervals."
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_rejection)

# ------------------------------------------------------------
# 5. Test-failure plot
# ------------------------------------------------------------

p_failure <- ggplot(
  rejection_plot_data,
  aes(
    x = N,
    y = failure_rate,
    shape = procedure,
    linetype = procedure,
    group = procedure
  )
) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(hypothesis ~ model) +
  scale_shape_manual(
    values = setNames(c(16, 17, 1, 2), procedure_levels),
    labels = procedure_labels
  ) +
  scale_linetype_manual(
    values = setNames(
      c("solid", "dashed", "dotted", "dotdash"),
      procedure_levels
    ),
    labels = procedure_labels
  ) +
  scale_x_continuous(breaks = sample_sizes) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    x = "Number of subjects",
    y = "Test failure rate",
    shape = "Procedure",
    linetype = "Procedure"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_failure)




library(dplyr)
library(ggplot2)

M <- nr_of_iters
alpha <- 0.05

mc_lower <- qbinom(0.025, M, alpha) / M
mc_upper <- qbinom(0.975, M, alpha) / M

plot_data <- rejection_summary |>
  mutate(
    model = factor(
      model, levels = c("GEE", "PGEE", "PIM", "PIM_adj")
    ),
    procedure = factor(
      procedure,
      levels = c(
        "Asymptotic_Global", "Asymptotic_Holm",
        "Finite_Global", "Finite_Holm"
      ),
      labels = c(
        "Global: chi-square", "Holm: z",
        "Global: F", "Holm: t"
      )
    )
  )

p_rejection <- ggplot(
  plot_data,
  aes(
    x = N, y = rejection_rate,
    shape = model, linetype = model, group = model
  )
) +
  geom_hline(
    yintercept = alpha,
    linetype = "dashed", colour = "grey40"
  ) +
  geom_hline(
    yintercept = c(mc_lower, mc_upper),
    linetype = "dotted", colour = "grey50"
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  facet_wrap(~ procedure, ncol = 2) +
  scale_shape_manual(
    values = c(GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2)
  ) +
  scale_linetype_manual(
    values = c(
      GEE = "solid", PGEE = "dashed",
      PIM = "dotted", PIM_adj = "dotdash"
    )
  ) +
  scale_x_continuous(
    breaks = sort(unique(plot_data$N)),
    guide = guide_axis(angle = 45)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, NA)
  ) +
  labs(
    x = "Number of subjects",
    y = "Rejection rate",
    shape = "Model",
    linetype = "Model",
    caption = paste0(
      "Dashed: 5%. Dotted: pointwise 95% Monte Carlo range; M = ", M, "."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_treatment <- p_rejection %+%
  filter(plot_data, hypothesis == "treatment") +
  labs(title = "Treatment-effect constancy")

p_time <- p_rejection %+%
  filter(plot_data, hypothesis == "time") +
  labs(title = "Equal temporal effects")

print(p_treatment)
print(p_time)






library(dplyr)

parameter_runs <- bind_rows(
  lapply(model_store, function(x) {
    
    b <- x$coefficients
    V <- x$vcov
    
    parameter_names <- names(b)
    if (is.null(parameter_names)) {
      parameter_names <- paste0("beta", seq_along(b))
    }
    
    stopifnot(
      nrow(V) == length(b),
      ncol(V) == length(b)
    )
    
    data.frame(
      N = x$N,
      run = x$run,
      model = x$model,
      parameter = seq_along(b),
      parameter_name = parameter_names,
      estimate = as.numeric(b),
      model_variance = diag(V),
      row.names = NULL
    )
  })
) |>
  mutate(
    valid = is.finite(estimate) &
      is.finite(model_variance) &
      model_variance > 0
  )

parameter_variance_comparison <- parameter_runs |>
  group_by(N, model, parameter, parameter_name) |>
  summarise(
    n_total = n(),
    n_valid = sum(valid),
    n_failed = sum(!valid),
    
    mean_estimate = if (any(valid)) {
      mean(estimate[valid])
    } else {
      NA_real_
    },
    
    empirical_variance = if (sum(valid) >= 2L) {
      var(estimate[valid])
    } else {
      NA_real_
    },
    
    average_model_variance = if (any(valid)) {
      mean(model_variance[valid])
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) |>
  mutate(
    failure_rate = n_failed / n_total,
    empirical_SD = sqrt(empirical_variance),
    model_SD = sqrt(average_model_variance),
    variance_ratio = empirical_variance / average_model_variance
  ) |>
  arrange(N, model, parameter)

print(parameter_variance_comparison, n = Inf)




library(ggplot2)

ggplot(
  parameter_variance_comparison,
  aes(
    x = N,
    y = variance_ratio,
    shape = model,
    linetype = model,
    group = model
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_point(size = 2.5) +
  facet_wrap(~ parameter, ncol = 2, labeller = label_both) +
  scale_shape_manual(
    values = c(GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2)
  ) +
  scale_linetype_manual(
    values = c(
      GEE = "solid", PGEE = "dashed",
      PIM = "dotted", PIM_adj = "dotdash"
    )
  ) +
  scale_x_continuous(
    breaks = sort(unique(parameter_variance_comparison$N)),
    guide = guide_axis(angle = 45)
  ) +
  labs(
    x = "Number of subjects",
    y = "Empirical variance / mean model variance",
    shape = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")


benchmark_runs <- dplyr::bind_rows(benchmark_store)

benchmark_summary <- benchmark_runs |>
  dplyr::group_by(N, hypothesis, procedure) |>
  dplyr::summarise(
    n_total = dplyr::n(),
    n_valid = sum(!is.na(reject)),
    n_failed = sum(is.na(reject)),
    n_rejected = sum(reject, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    failure_rate = n_failed / n_total,
    rejection_rate = n_rejected / dplyr::na_if(n_valid, 0L),
    MCSE = sqrt(
      rejection_rate * (1 - rejection_rate) /
        dplyr::na_if(n_valid, 0L)
    )
  )

print(benchmark_summary, n = Inf)



## Appendix benchmark plots

library(dplyr)
library(ggplot2)

M <- nr_of_iters
alpha <- 0.05

mc_lower <- qbinom(0.025, M, alpha) / M
mc_upper <- qbinom(0.975, M, alpha) / M

# PIM treatment-constancy results
pim_plot_data <- rejection_summary |>
  filter(hypothesis == "treatment") |>
  transmute(
    N,
    rejection_rate,
    method = as.character(model),
    panel = recode(
      as.character(procedure),
      Asymptotic_Global = "Global: chi-square",
      Asymptotic_Holm   = "Holm: z",
      Finite_Global    = "Global: F",
      Finite_Holm      = "Holm: t"
    )
  )

# nparLD interaction results
benchmark_plot_data <- benchmark_summary |>
  transmute(
    N,
    rejection_rate,
    method = recode(
      as.character(procedure),
      nparLD_WTS = "nparLD WTS",
      nparLD_ATS = "nparLD ATS"
    ),
    panel = "nparLD: RTE interaction"
  )

appendix_plot_data <- bind_rows(
  pim_plot_data,
  benchmark_plot_data
) |>
  mutate(
    panel = factor(
      panel,
      levels = c(
        "Global: chi-square", "Holm: z",
        "Global: F", "Holm: t",
        "nparLD: RTE interaction"
      )
    )
  )

p_appendix <- ggplot(
  appendix_plot_data,
  aes(
    x = N,
    y = rejection_rate,
    shape = method,
    linetype = method,
    group = method
  )
) +
  geom_hline(
    yintercept = alpha,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_hline(
    yintercept = c(mc_lower, mc_upper),
    linetype = "dotted",
    colour = "grey50"
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  facet_wrap(~ panel, ncol = 2) +
  scale_shape_manual(
    values = c(
      GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2,
      "nparLD WTS" = 15, "nparLD ATS" = 0
    )
  ) +
  scale_linetype_manual(
    values = c(
      GEE = "solid",
      PGEE = "dashed",
      PIM = "dotted",
      PIM_adj = "dotdash",
      "nparLD WTS" = "solid",
      "nparLD ATS" = "dashed"
    )
  ) +
  scale_x_continuous(
    breaks = sort(unique(appendix_plot_data$N)),
    guide = guide_axis(angle = 45)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, NA)
  ) +
  labs(
    x = "Number of subjects",
    y = "Rejection rate",
    shape = "Method",
    linetype = "Method",
    caption = paste0(
      "PIM panels: treatment-effect constancy; nparLD panel: RTE interaction.\n",
      "Dashed: 5%. Dotted: pointwise 95% binomial reference range; M = ",
      M, "."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(
    shape = guide_legend(nrow = 2),
    linetype = guide_legend(nrow = 2)
  )

print(p_appendix)




lachin_runs <- dplyr::bind_rows(lachin_store)

lachin_summary <- lachin_runs |>
  dplyr::group_by(N, model, hypothesis, procedure) |>
  dplyr::summarise(
    n_total = dplyr::n(),
    n_valid = sum(!is.na(reject)),
    n_failed = sum(is.na(reject)),
    n_rejected = sum(reject, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    failure_rate = n_failed / n_total,
    rejection_rate = n_rejected / dplyr::na_if(n_valid, 0L),
    MCSE = sqrt(
      rejection_rate * (1 - rejection_rate) /
        dplyr::na_if(n_valid, 0L)
    )
  )

print(lachin_summary, n = Inf)





library(dplyr)
library(ggplot2)

# These must match the degrees of freedom used in your tests
df_trt  <- function(N) N - 1
df_time <- function(N) N - 1

density_data <- bind_rows(contrast_store) |>
  filter(
    is.finite(estimate),
    is.finite(variance),
    variance > 0
  ) |>
  mutate(
    statistic = estimate / sqrt(variance),
    nu = if_else(
      hypothesis == "treatment",
      df_trt(N),
      df_time(N)
    ),
    contrast_label = paste("Contrast", contrast)
  )

# Reference densities for each panel
reference_data <- density_data |>
  distinct(N, hypothesis, contrast_label, nu) |>
  group_by(N, hypothesis, contrast_label, nu) |>
  reframe(
    x = seq(-5, 5, length.out = 501),
    t_density = dt(x, df = first(nu)),
    normal_density = dnorm(x)
  )

plot_contrast_density <- function(hypothesis_name) {
  ggplot(
    filter(density_data, hypothesis == hypothesis_name),
    aes(x = statistic, colour = model)
  ) +
    geom_density(linewidth = 0.7) +
    geom_line(
      data = filter(reference_data, hypothesis == hypothesis_name),
      aes(x = x, y = t_density),
      inherit.aes = FALSE,
      colour = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_line(
      data = filter(reference_data, hypothesis == hypothesis_name),
      aes(x = x, y = normal_density),
      inherit.aes = FALSE,
      colour = "black",
      linetype = "dotted",
      linewidth = 0.8
    ) +
    facet_grid(contrast_label ~ N, labeller = label_both) +
    coord_cartesian(xlim = c(-5, 5)) +
    labs(
      x = "Estimated contrast / estimated SE",
      y = "Density",
      colour = "Model",
      caption = "Dashed: t reference; dotted: standard normal reference."
    ) +
    theme_bw()
}

plot_contrast_density("treatment")
plot_contrast_density("time")


save.image("sim_PO_test.Rdata")



library(dplyr)
library(ggplot2)

# Labels follow your coefficient order
parameter_labels <- c(
  "Time: treated", "Time: control",
  "Treatment: visit 1", "Treatment: visit 2",
  "Treatment: visit 3", "Treatment: visit 4"
)

records <- Filter(function(x) {
  !is.null(x) && x$model %in% c("PIM", "GEE", "PGEE")
}, model_store)

keys <- bind_rows(lapply(records, function(x) {
  data.frame(N = x$N, run = x$run, model = x$model)
}))

stopifnot(nrow(keys) > 0L, !anyDuplicated(keys))

groups <- split(
  records,
  interaction(keys$N, keys$model, drop = TRUE)
)

empirical_covariances <- bind_rows(lapply(groups, function(g) {
  
  p <- length(parameter_labels)
  
  valid <- vapply(g, function(x) {
    length(x$coefficients) == p &&
      all(is.finite(x$coefficients))
  }, logical(1))
  
  good <- g[valid]
  n_valid <- length(good)
  
  V <- matrix(NA_real_, p, p)
  
  if (n_valid >= 2L) {
    B <- do.call(rbind, lapply(good, function(x) {
      as.numeric(x$coefficients)
    }))
    V <- cov(B)
  }
  
  ij <- expand.grid(i = seq_len(p), j = seq_len(p))
  
  data.frame(
    N = g[[1]]$N,
    model = g[[1]]$model,
    n_valid = n_valid,
    n_excluded = length(g) - n_valid,
    parameter_i = parameter_labels[ij$i],
    parameter_j = parameter_labels[ij$j],
    covariance = V[cbind(ij$i, ij$j)]
  )
})) |>
  mutate(
    parameter_i = factor(parameter_i, levels = parameter_labels),
    parameter_j = factor(parameter_j, levels = rev(parameter_labels)),
    model = factor(model, levels = c("PIM", "GEE", "PGEE")),
    N = factor(N, levels = sort(unique(N)))
  )



ggplot(
  empirical_covariances,
  aes(x = parameter_i, y = parameter_j, fill = covariance)
) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Empirical\ncovariance"
  ) +
  facet_grid(model ~ N) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )



ggplot(
  empirical_covariances,
  aes(
    x = parameter_i,
    y = parameter_j,
    fill = as.numeric(as.character(N)) * covariance
  )
) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    name = "N × empirical\ncovariance"
  ) +
  facet_grid(model ~ N) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )


tail_check <- density_data |>
  dplyr::group_by(N, model, hypothesis, contrast_label) |>
  dplyr::summarise(
    n_valid = dplyr::n(),
    mean_T = mean(statistic),
    SD_T = sd(statistic),
    
    # Individual test: expected 0.05 under calibration
    rejection_05 = mean(
      abs(statistic) > qt(0.975, df = nu)
    ),
    
    # First Holm threshold for four treatment contrasts:
    # expected 0.0125 under calibration
    rejection_0125 = mean(
      abs(statistic) > qt(0.99375, df = nu)
    ),
    .groups = "drop"
  )

tail_check |>
  dplyr::filter(hypothesis == "treatment") |>
  print(n = Inf)



tail_check |>
  dplyr::filter(model == "PGEE", hypothesis == "treatment") |>
  print(n = Inf)

rejection_summary |>
  dplyr::filter(model == "PGEE", hypothesis == "treatment") |>
  print(n = Inf)



library(dplyr)

raw_contrast_check <- bind_rows(contrast_store) |>
  filter(
    model == "PGEE",
    hypothesis == "treatment",
    is.finite(estimate)
  ) |>
  group_by(N, contrast) |>
  summarise(
    n_valid = n(),
    mean_contrast = mean(estimate),
    empirical_SD = sd(estimate),
    MCSE_mean = empirical_SD / sqrt(n_valid),
    standardized_bias = mean_contrast / empirical_SD,
    .groups = "drop"
  ) |>
  mutate(
    lower_mean = mean_contrast - 1.96 * MCSE_mean,
    upper_mean = mean_contrast + 1.96 * MCSE_mean
  )

print(raw_contrast_check, n = Inf)


# Select independent rows without changing the hypothesis
qr_L <- qr(t(L_trt))
r <- qr_L$rank

stopifnot(r == 3L)

selected_rows <- qr_L$pivot[seq_len(r)]
L3 <- L_trt[selected_rows, , drop = FALSE]

records <- Filter(function(x) {
  !is.null(x) && x$model == "PGEE"
}, model_store)

stopifnot(length(records) > 0L)

N_values <- sort(unique(vapply(records, function(x) x$N, numeric(1))))

covariance_details <- setNames(
  lapply(N_values, function(N_current) {
    
    g <- Filter(function(x) x$N == N_current, records)
    
    stopifnot(!anyDuplicated(vapply(
      g, function(x) x$run, numeric(1)
    )))
    
    # Use the same runs for empirical and estimated covariance
    valid <- vapply(g, function(x) {
      length(x$coefficients) == ncol(L3) &&
        all(is.finite(x$coefficients)) &&
        all(dim(as.matrix(x$vcov)) == c(ncol(L3), ncol(L3))) &&
        all(is.finite(x$vcov))
    }, logical(1))
    
    good <- g[valid]
    
    if (length(good) < 2L)
      stop("Fewer than two valid runs for N = ", N_current)
    
    estimates <- do.call(rbind, lapply(good, function(x) {
      drop(L3 %*% x$coefficients)
    }))
    
    contrast_V <- lapply(good, function(x) {
      V <- as.matrix(x$vcov)
      V <- (V + t(V)) / 2
      L3 %*% V %*% t(L3)
    })
    
    empirical_V <- cov(estimates)
    average_V <- Reduce(`+`, contrast_V) / length(contrast_V)
    
    labels <- paste0("Contrast ", selected_rows)
    dimnames(empirical_V) <- dimnames(average_V) <-
      list(labels, labels)
    
    list(
      N = N_current,
      n_valid = length(good),
      n_excluded = sum(!valid),
      empirical = empirical_V,
      average_estimated = average_V,
      difference = average_V - empirical_V,
      
      # Run-specific positive-definiteness diagnostic
      minimum_eigenvalues = vapply(contrast_V, function(V) {
        min(eigen(V, symmetric = TRUE, only.values = TRUE)$values)
      }, numeric(1))
    )
  }),
  as.character(N_values)
)

# Inspect matrices at the smallest sample size
covariance_details[[as.character(min(N_values))]][
  c("N", "n_valid", "n_excluded",
    "empirical", "average_estimated", "difference")
]

  
  covariance_details[["16"]][
    c("empirical", "average_estimated", "difference")
  ]

  x <- covariance_details[["16"]]
  
  # Standardize by the average estimated covariance
  R <- chol(x$average_estimated)
  R_inv <- backsolve(R, diag(nrow(R)))
  
  relative_V <- t(R_inv) %*% x$empirical %*% R_inv
  
  eigen(
    relative_V,
    symmetric = TRUE,
    only.values = TRUE
  )$values  
  

  
  
  
  library(dplyr)
  
  parameter_runs <- bind_rows(
    lapply(model_store, function(x) {
      
      b <- x$coefficients
      V <- x$vcov
      
      parameter_names <- names(b)
      if (is.null(parameter_names)) {
        parameter_names <- paste0("beta", seq_along(b))
      }
      
      stopifnot(
        nrow(V) == length(b),
        ncol(V) == length(b)
      )
      
      data.frame(
        N = x$N,
        run = x$run,
        model = x$model,
        parameter = seq_along(b),
        parameter_name = parameter_names,
        estimate = as.numeric(b),
        model_variance = diag(V),
        row.names = NULL
      )
    })
  ) |>
    mutate(
      valid = is.finite(estimate) &
        is.finite(model_variance) &
        model_variance > 0
    )
  
  parameter_variance_comparison <- parameter_runs |>
    group_by(N, model, parameter, parameter_name) |>
    summarise(
      n_total = n(),
      n_valid = sum(valid),
      n_failed = sum(!valid),
      
      mean_estimate = if (any(valid)) {
        mean(estimate[valid])
      } else {
        NA_real_
      },
      
      empirical_variance = if (sum(valid) >= 2L) {
        var(estimate[valid])
      } else {
        NA_real_
      },
      
      average_model_variance = if (any(valid)) {
        mean(model_variance[valid])
      } else {
        NA_real_
      },
      
      .groups = "drop"
    ) |>
    mutate(
      failure_rate = n_failed / n_total,
      empirical_SD = sqrt(empirical_variance),
      model_SD = sqrt(average_model_variance),
      variance_ratio = empirical_variance / average_model_variance
    ) |>
    arrange(N, model, parameter)
  
  print(parameter_variance_comparison, n = Inf)
  
  
  
  
  library(ggplot2)
  
  ggplot(
    parameter_variance_comparison,
    aes(
      x = N,
      y = variance_ratio,
      shape = model,
      linetype = model,
      group = model
    )
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line() +
    geom_point(size = 2.5) +
    facet_wrap(~ parameter, ncol = 2, labeller = label_both) +
    scale_shape_manual(
      values = c(GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2)
    ) +
    scale_linetype_manual(
      values = c(
        GEE = "solid", PGEE = "dashed",
        PIM = "dotted", PIM_adj = "dotdash"
      )
    ) +
    scale_x_continuous(
      breaks = sort(unique(parameter_variance_comparison$N)),
      guide = guide_axis(angle = 45)
    ) +
    labs(
      x = "Number of subjects",
      y = "Empirical variance / mean model variance",
      shape = "Model",
      linetype = "Model"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  
  
  
  
  
## Benchmarks
  
  library(dplyr)
  library(ggplot2)
  
  # Summarise benchmark rejection rates
  summarise_benchmark <- function(store, benchmark_name) {
    bind_rows(store) |>
      group_by(N, procedure) |>
      summarise(
        n_total = n(),
        n_valid = sum(!is.na(reject)),
        n_failed = sum(is.na(reject)),
        n_rejected = sum(reject, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        rejection_rate = n_rejected / na_if(n_valid, 0L),
        model = benchmark_name
      )
  }
  
  lachin_summary <- summarise_benchmark(lachin_store, "Lachin")
  nparld_summary <- summarise_benchmark(benchmark_store, "nparLD")
  
  # Check the stored procedure names
  print(lachin_summary, n = Inf)
  print(nparld_summary, n = Inf)
  
  # Monte Carlo reference band for a true rejection probability of 5%
  M <- nr_of_iters
  
  mc_lower <- qbinom(0.025, M, 0.05) / M
  mc_upper <- qbinom(0.975, M, 0.05) / M
  
  
  
  p_pim <- ggplot(
    filter(rejection_summary, hypothesis == "treatment"),
    aes(
      x = N, y = rejection_rate,
      group = model, shape = model, linetype = model
    )
  ) +
    geom_line(colour = "black") +
    geom_point(colour = "black", size = 2) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_hline(
      yintercept = c(mc_lower, mc_upper),
      linetype = "dotted"
    ) +
    scale_shape_manual(values = c(
      GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2
    )) +
    scale_linetype_manual(values = c(
      GEE = "solid", PGEE = "dashed",
      PIM = "dotted", PIM_adj = "dotdash"
    )) +
    facet_wrap(~ procedure, ncol = 2) +
    labs(
      title = "Treatment-PI constancy",
      x = "Number of subjects",
      y = "Rejection rate",
      shape = "Model",
      linetype = "Model"
    ) +
    theme_bw()
  
  benchmark_plot_data <- bind_rows(
    lachin_summary |>
      mutate(panel = "Lachin: treatment-PI constancy"),
    nparld_summary |>
      mutate(panel = "nparLD: relative-effect interaction")
  ) |>
    mutate(
      curve = paste(model, procedure, sep = ": ")
    )
  
  p_benchmarks <- ggplot(
    benchmark_plot_data,
    aes(
      x = N, y = rejection_rate,
      group = curve, shape = curve, linetype = curve
    )
  ) +
    geom_line(colour = "black") +
    geom_point(colour = "black", size = 2) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_hline(
      yintercept = c(mc_lower, mc_upper),
      linetype = "dotted"
    ) +
    facet_wrap(~ panel, ncol = 2) +
    labs(
      title = "Benchmark tests",
      x = "Number of subjects",
      y = "Rejection rate",
      shape = "Test",
      linetype = "Test",
      caption = paste(
        "Dashed horizontal line: 5%. Dotted lines: Monte Carlo reference bounds.",
        "The nparLD interaction null differs from treatment-PI constancy."
      )
    ) +
    theme_bw()
  
  p_pim
  p_benchmarks

  
  
  
  
  # Add Lachin benchmark
  
  library(dplyr)
  library(ggplot2)
  
  # Lachin: one result per dataset
  lachin_summary <- bind_rows(lachin_store) |>
    group_by(N) |>
    summarise(
      n_valid = sum(!is.na(reject)),
      n_failed = sum(is.na(reject)),
      n_rejected = sum(reject, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      rejection_rate = n_rejected / na_if(n_valid, 0L)
    )
  
  # Repeat the same Lachin benchmark in each treatment procedure panel
  procedures <- rejection_summary |>
    distinct(procedure)
  
  lachin_overlay <- cross_join(lachin_summary, procedures) |>
    mutate(
      hypothesis = "treatment",
      model = "Lachin"
    )
  
  plot_data <- bind_rows(
    rejection_summary |>
      filter(hypothesis %in% c("treatment", "time")) |>
      dplyr::select(N, hypothesis, procedure, model, rejection_rate),
    
    lachin_overlay |>
      dplyr::select(N, hypothesis, procedure, model, rejection_rate)
  ) |>
    mutate(
      hypothesis = factor(
        hypothesis,
        levels = c("treatment", "time"),
        labels = c(
          "Treatment-effect constancy",
          "Equal temporal effects"
        )
      ),
      model = factor(
        model,
        levels = c("GEE", "PGEE", "PIM", "PIM_adj", "Lachin")
      )
    )
  
  M <- nr_of_iters
  mc_lower <- qbinom(0.025, M, 0.05) / M
  mc_upper <- qbinom(0.975, M, 0.05) / M
  
  ggplot(
    plot_data,
    aes(
      x = N, y = rejection_rate,
      group = model, shape = model, linetype = model
    )
  ) +
    geom_line(colour = "black", linewidth = 0.6) +
    geom_point(colour = "black", size = 2) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_hline(
      yintercept = c(mc_lower, mc_upper),
      linetype = "dotted"
    ) +
    scale_shape_manual(values = c(
      GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2, Lachin = 8
    )) +
    scale_linetype_manual(values = c(
      GEE = "solid",
      PGEE = "dashed",
      PIM = "dotted",
      PIM_adj = "dotdash",
      Lachin = "longdash"
    )) +
    scale_x_continuous(breaks = sort(unique(plot_data$N))) +
    facet_grid(hypothesis ~ procedure) +
    labs(
      x = "Number of subjects",
      y = "Rejection rate",
      shape = "Method",
      linetype = "Method",
      caption = paste(
        "Lachin's treatment-constancy benchmark is repeated unchanged",
        "across the treatment panels.",
        "Dotted horizontal lines: 95% Monte Carlo reference bounds at 5%."
      )
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  
  
  
  
  
  

  library(multcomp)
  library(dplyr)
  
  records <- Filter(function(x) {
    !is.null(x) && x$N >= 16 && x$N <= 60
  }, model_store)
  
  # Run once per stored model fit
  test_stored_fit <- function(j) {
    
    fit <- records[[j]]
    b <- as.numeric(fit$coefficients)
    V <- as.matrix(fit$vcov)
    L <- as.matrix(L_trt)
    
    setup <- tryCatch({
      stopifnot(
        length(b) == ncol(L),
        identical(dim(V), c(length(b), length(b))),
        all(is.finite(b)),
        all(is.finite(V))
      )
      
      V <- (V + t(V)) / 2
      S <- L %*% V %*% t(L)
      S <- (S + t(S)) / 2
      
      if (any(diag(S) <= 0))
        stop("Non-positive contrast variance.")
      
      ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
      
      # Allow the expected singularity from redundant contrasts
      if (min(ev) < -1e-10 * max(abs(ev)))
        stop("Contrast covariance is not positive semidefinite.")
      
      TRUE
    }, error = function(e) e)
    
    family_rows <- list()
    elementary_rows <- list()
    k <- 0L
    
    for (reference in c("Asymptotic", "Finite")) {
      
      nu <- if (reference == "Asymptotic") 0 else fit$N - 1
      
      g <- tryCatch({
        if (inherits(setup, "error"))
          stop(conditionMessage(setup))
        
        glht(
          parm(coef = b, vcov = V, df = nu),
          linfct = L,
          alternative = "two.sided"
        )
      }, error = function(e) e)
      
      for (method in c("single-step", "Westfall")) {
        
        k <- k + 1L
        warnings <- character()
        
        result <- tryCatch({
          if (inherits(g, "error"))
            stop(conditionMessage(g))
          
          # Reproducible numerical integration
          set.seed(100000L + 4L * j + k)
          
          out <- withCallingHandlers(
            summary(
              g,
              test = adjusted(
                method,
                maxpts = 10000L,
                abseps = 1e-4
              )
            ),
            warning = function(w) {
              warnings <<- c(warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          )
          
          p <- as.numeric(out$test$pvalues)
          
          if (length(p) != nrow(L) ||
              any(!is.finite(p)) ||
              any(p < 0 | p > 1)) {
            stop("Invalid adjusted p-values.")
          }
          
          integration_error <- out$test$error
          if (!is.numeric(integration_error) ||
              length(integration_error) != 1L) {
            integration_error <- NA_real_
          }
          
          list(
            p = p,
            status = if (length(warnings)) "warning" else "ok",
            error = NA_character_,
            integration_error = integration_error
          )
          
        }, error = function(e) {
          list(
            p = rep(NA_real_, nrow(L)),
            status = "failed",
            error = conditionMessage(e),
            integration_error = NA_real_
          )
        })
        
        method_label <- if (method == "single-step") {
          "SingleStep"
        } else {
          "Westfall"
        }
        
        procedure <- paste(reference, method_label, sep = "_")
        valid <- all(is.finite(result$p))
        
        family_rows[[k]] <- data.frame(
          N = fit$N,
          run = fit$run,
          model = fit$model,
          hypothesis = "treatment",
          procedure = procedure,
          df = nu,
          p_family = if (valid) min(result$p) else NA_real_,
          reject = if (valid) any(result$p <= 0.05) else NA,
          status = result$status,
          error = result$error,
          warning = paste(unique(warnings), collapse = "; "),
          integration_error = result$integration_error
        )
        
        elementary_rows[[k]] <- data.frame(
          N = fit$N,
          run = fit$run,
          model = fit$model,
          hypothesis = "treatment",
          procedure = procedure,
          contrast = seq_len(nrow(L)),
          p_adjusted = result$p,
          reject = result$p <= 0.05
        )
      }
    }
    
    if (j %% 100L == 0L)
      cat("Processed", j, "of", length(records), "stored fits\n")
    
    list(
      family = bind_rows(family_rows),
      elementary = bind_rows(elementary_rows)
    )
  }
  
  # Preserve RNG state, including if processing stops with an error
  run_stored_tests <- function() {
    had_seed <- exists(
      ".Random.seed", envir = .GlobalEnv, inherits = FALSE
    )
    
    if (had_seed)
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
    
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(
        ".Random.seed", envir = .GlobalEnv, inherits = FALSE
      )) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    
    lapply(seq_along(records), test_stored_fit)
  }
  
  results <- run_stored_tests()
  
  joint_test_runs <- bind_rows(lapply(results, `[[`, "family"))
  joint_elementary_runs <- bind_rows(
    lapply(results, `[[`, "elementary")
  )
  
  
  
  # Install once if needed:
  # install.packages(c("future", "future.apply", "parallelly"))
  
  library(future)
  library(future.apply)
  library(parallelly)
  
  # Leave one available core free
  n_workers <- max(1L, parallelly::availableCores() - 1L)
  
  # Start conservatively: each worker needs memory
  n_workers <- min(n_workers, 4L)
  
  future::plan(
    future::multisession,
    workers = n_workers
  )
  
  results <- tryCatch(
    future.apply::future_lapply(
      seq_along(records),
      FUN = test_stored_fit,
      future.seed = TRUE,
      future.packages = c("multcomp", "dplyr"),
      future.scheduling = 2
    ),
    finally = {
      future::plan(future::sequential)
    }
  )
  
  joint_test_runs <- dplyr::bind_rows(
    lapply(results, `[[`, "family")
  )
  
  joint_elementary_runs <- dplyr::bind_rows(
    lapply(results, `[[`, "elementary")
  )
  
  

  library(dplyr)
  
  holm_runs <- bind_rows(test_store) |>
    filter(
      N >= 16, N <= 60,
      hypothesis == "treatment",
      procedure %in% c("Asymptotic_Holm", "Finite_Holm")
    ) |>
    rename(p_family = p_global)
  
  multiplicity_runs <- bind_rows(
    holm_runs,
    joint_test_runs
  )
  
  stopifnot(
    !anyDuplicated(
      multiplicity_runs[
        c("N", "run", "model", "hypothesis", "procedure")
      ]
    )
  )
  
  multiplicity_summary <- multiplicity_runs |>
    group_by(N, model, hypothesis, procedure) |>
    summarise(
      n_total = n(),
      n_valid = sum(!is.na(reject)),
      n_failed = sum(is.na(reject)),
      n_rejected = sum(reject, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      failure_rate = n_failed / n_total,
      rejection_rate = n_rejected / na_if(n_valid, 0L),
      MCSE = sqrt(
        rejection_rate * (1 - rejection_rate) /
          na_if(n_valid, 0L)
      )
    )
  
  print(multiplicity_summary, n = Inf)
  
  
  

  
  library(dplyr)
  library(ggplot2)
  
  plot_data <- multiplicity_summary |>
    mutate(
      reference = if_else(
        grepl("^Asymptotic_", procedure),
        "Normal reference",
        "t reference"
      ),
      adjustment = sub("^(Asymptotic|Finite)_", "", procedure),
      adjustment = factor(
        adjustment,
        levels = c("Holm", "SingleStep", "Westfall")
      ),
      reference = factor(
        reference,
        levels = c("Normal reference", "t reference")
      ),
      model = factor(
        model,
        levels = c("PIM", "PIM_adj", "GEE", "PGEE")
      ),
      
      # Pointwise MC reference bounds under a true error rate of 5%.
      # Use the number of valid runs for each result.
      mc_lower = qbinom(0.025, size = n_valid, prob = 0.05) /
        na_if(n_valid, 0L),
      mc_upper = qbinom(0.975, size = n_valid, prob = 0.05) /
        na_if(n_valid, 0L)
    )
  
  ggplot(
    plot_data,
    aes(
      x = N, y = rejection_rate,
      group = adjustment,
      shape = adjustment,
      linetype = adjustment
    )
  ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      colour = "grey40"
    ) +
    geom_line(
      aes(y = mc_lower),
      linetype = "dotted",
      colour = "grey65",
      linewidth = 0.4,
      show.legend = FALSE
    ) +
    geom_line(
      aes(y = mc_upper),
      linetype = "dotted",
      colour = "grey65",
      linewidth = 0.4,
      show.legend = FALSE
    ) +
    geom_line(colour = "black", linewidth = 0.6) +
    geom_point(colour = "black", size = 2.5) +
    facet_grid(reference ~ model) +
    scale_shape_manual(
      values = c(Holm = 16, SingleStep = 1, Westfall = 17)
    ) +
    scale_linetype_manual(
      values = c(
        Holm = "solid",
        SingleStep = "dashed",
        Westfall = "dotdash"
      )
    ) +
    scale_x_continuous(
      breaks = sort(unique(plot_data$N)),
      guide = guide_axis(angle = 45)
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, NA)
    ) +
    labs(
      title = "Treatment-effect constancy",
      x = "Number of subjects",
      y = "Empirical Type I error",
      shape = "Adjustment",
      linetype = "Adjustment",
      caption = paste(
        "Dashed horizontal line: 5%.",
        "Dotted curves: pointwise 95% Monte Carlo reference bounds.",
        "Rates use valid runs."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  
  
  library(dplyr)
  library(ggplot2)
  
  # Existing global tests
  global_summary <- rejection_summary |>
    filter(
      N >= 16, N <= 60,
      hypothesis == "treatment",
      procedure %in% c("Asymptotic_Global", "Finite_Global")
    ) |>
    transmute(
      N, model, n_valid, rejection_rate,
      reference = if_else(
        procedure == "Asymptotic_Global",
        "Asymptotic", "Finite approximation"
      ),
      test = "Global"
    )
  
  # Holm, single-step and Westfall
  multiple_summary <- multiplicity_summary |>
    filter(
      N >= 16, N <= 60,
      hypothesis == "treatment"
    ) |>
    transmute(
      N, model, n_valid, rejection_rate,
      reference = if_else(
        grepl("^Asymptotic_", procedure),
        "Asymptotic", "Finite approximation"
      ),
      test = sub("^(Asymptotic|Finite)_", "", procedure)
    )
  
  main_plot_data <- bind_rows(global_summary, multiple_summary)
  
  # Lachin benchmark
  lachin_plot_data <- bind_rows(lachin_store) |>
    filter(N >= 16, N <= 60) |>
    group_by(N) |>
    summarise(
      n_valid = sum(!is.na(reject)),
      n_rejected = sum(reject, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      rejection_rate = n_rejected / na_if(n_valid, 0L),
      model = "Lachin"
    )
  
  # Repeat Lachin in each panel without changing its test
  panels <- main_plot_data |>
    distinct(reference, test)
  
  plot_data <- bind_rows(
    main_plot_data,
    cross_join(lachin_plot_data, panels) |>
      dplyr::select(
        N, model, n_valid, rejection_rate, reference, test
      )
  ) |>
    mutate(
      reference = factor(
        reference,
        levels = c("Asymptotic", "Finite approximation")
      ),
      test = factor(
        test,
        levels = c("Global", "Holm", "SingleStep", "Westfall"),
        labels = c("Global", "Holm", "Single-step", "Westfall")
      ),
      model = factor(
        model,
        levels = c("GEE", "PGEE", "PIM", "PIM_adj", "Lachin")
      )
    )
  
  # Horizontal Monte Carlo reference lines
  M <- nr_of_iters
  mc_lower <- qbinom(0.025, M, 0.05) / M
  mc_upper <- qbinom(0.975, M, 0.05) / M
  
  p_rejection <- ggplot(
    plot_data%>%filter(test!="Westfall"),
    aes(
      x = N, y = rejection_rate,
      group = model, shape = model, linetype = model
    )
  ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      colour = "grey40"
    ) +
    geom_hline(
      yintercept = c(mc_lower, mc_upper),
      linetype = "dotted",
      colour = "grey60"
    ) +
    geom_line(colour = "black", linewidth = 0.6) +
    geom_point(colour = "black", size = 2.3) +
    facet_grid(reference ~ test) +
    scale_shape_manual(values = c(
      GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2, Lachin = 8
    )) +
    scale_linetype_manual(values = c(
      GEE = "solid",
      PGEE = "dashed",
      PIM = "dotted",
      PIM_adj = "dotdash",
      Lachin = "longdash"
    )) +
    scale_x_continuous(
      breaks = sort(unique(plot_data$N)),
      guide = guide_axis(angle = 45)
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, NA)
    ) +
    labs(
      title = "Treatment-effect constancy",
      x = "Number of subjects",
      y = "Empirical Type I error",
      shape = "Method",
      linetype = "Method",
      caption = paste0(
        "Global references: chi-square (top), F (bottom). ",
        "Multiplicity references: z (top), t (bottom).\n",
        "Lachin remains the same asymptotic benchmark in every panel. ",
        "Dotted lines: 95% Monte Carlo reference range for M = ", M, "."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  print(p_rejection)
  
  save.image("PO_Single_Step.Rdata")

  
  ## Add all Tukey comparisons with single step approach
  
  future::plan(
       future::multisession,
       workers = n_workers
    )
  
  
  library(multcomp)
  library(dplyr)
  
  # Current coefficient order:
  # time treated, time control, treatment at each visit
  p <- ncol(L_trt)
  treatment_columns <- 3:p
  
  visit_pairs <- t(combn(seq_along(treatment_columns), 2))
  
  L_tukey <- matrix(
    0,
    nrow = nrow(visit_pairs),
    ncol = p
  )
  
  for (i in seq_len(nrow(visit_pairs))) {
    L_tukey[i, treatment_columns[visit_pairs[i, 1]]] <-  1
    L_tukey[i, treatment_columns[visit_pairs[i, 2]]] <- -1
  }
  
  rownames(L_tukey) <- paste0(
    "Visit ", visit_pairs[, 1],
    " - Visit ", visit_pairs[, 2]
  )
  
  L_tukey  
  
  
  
  test_stored_fit <- function(j) {
    
    fit <- records[[j]]
    b <- as.numeric(fit$coefficients)
    V <- as.matrix(fit$vcov)
    
    contrast_families <- list(
      SingleStep = L_trt,     # Visit versus average
      Tukey = L_tukey         # All pairwise visit effects
    )
    
    family_rows <- list()
    elementary_rows <- list()
    k <- 0L
    
    for (family_name in names(contrast_families)) {
      
      L <- contrast_families[[family_name]]
      
      for (reference in c("Asymptotic", "Finite")) {
        
        k <- k + 1L
        warnings <- character()
        nu <- if (reference == "Asymptotic") 0 else fit$N - 1
        
        result <- tryCatch({
          
          stopifnot(
            length(b) == ncol(L),
            identical(dim(V), c(length(b), length(b))),
            all(is.finite(b)),
            all(is.finite(V))
          )
          
          V_symmetric <- (V + t(V)) / 2
          
          S <- L %*% V_symmetric %*% t(L)
          S <- (S + t(S)) / 2
          
          if (any(diag(S) <= 0))
            stop("Non-positive contrast variance.")
          
          ev <- eigen(
            S, symmetric = TRUE, only.values = TRUE
          )$values
          
          # Redundant contrasts produce expected zero eigenvalues
          if (min(ev) < -1e-10 * max(abs(ev)))
            stop("Contrast covariance is not positive semidefinite.")
          
          g <- glht(
            parm(coef = b, vcov = V_symmetric, df = nu),
            linfct = L,
            alternative = "two.sided"
          )
          
          set.seed(100000L + 4L * j + k)
          
          out <- withCallingHandlers(
            summary(
              g,
              test = adjusted(
                "single-step",
                maxpts = 10000L,
                abseps = 1e-4
              )
            ),
            warning = function(w) {
              warnings <<- c(warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          )
          
          adjusted_p <- as.numeric(out$test$pvalues)
          
          if (length(adjusted_p) != nrow(L) ||
              any(!is.finite(adjusted_p)) ||
              any(adjusted_p < 0 | adjusted_p > 1)) {
            stop("Invalid adjusted p-values.")
          }
          
          integration_error <- out$test$error
          
          if (!is.numeric(integration_error) ||
              length(integration_error) != 1L) {
            integration_error <- NA_real_
          }
          
          list(
            p = adjusted_p,
            status = if (length(warnings)) "warning" else "ok",
            error = NA_character_,
            integration_error = integration_error
          )
          
        }, error = function(e) {
          list(
            p = rep(NA_real_, nrow(L)),
            status = "failed",
            error = conditionMessage(e),
            integration_error = NA_real_
          )
        })
        
        valid <- all(is.finite(result$p))
        procedure <- paste(reference, family_name, sep = "_")
        
        family_rows[[k]] <- data.frame(
          N = fit$N,
          run = fit$run,
          model = fit$model,
          hypothesis = "treatment",
          procedure = procedure,
          df = nu,
          p_family = if (valid) min(result$p) else NA_real_,
          reject = if (valid) any(result$p <= 0.05) else NA,
          status = result$status,
          error = result$error,
          warning = paste(unique(warnings), collapse = "; "),
          integration_error = result$integration_error
        )
        
        labels <- rownames(L)
        if (is.null(labels)) {
          labels <- paste("Visit", seq_len(nrow(L)), "- mean")
        }
        
        elementary_rows[[k]] <- data.frame(
          N = fit$N,
          run = fit$run,
          model = fit$model,
          hypothesis = "treatment",
          procedure = procedure,
          contrast = seq_len(nrow(L)),
          contrast_label = labels,
          p_adjusted = result$p,
          reject = result$p <= 0.05
        )
      }
    }
    
    list(
      family = bind_rows(family_rows),
      elementary = bind_rows(elementary_rows)
    )
  }
  
  
  
  
  results <- tryCatch(
    future.apply::future_lapply(
      seq_along(records),
      FUN = test_stored_fit,
      future.seed = TRUE,
      future.packages = c("multcomp", "dplyr"),
      future.scheduling = 2
    ),
    finally = {
      future::plan(future::sequential)
    }
  )
  
  joint_test_runs <- bind_rows(lapply(results, `[[`, "family"))
  joint_elementary_runs <- bind_rows(
    lapply(results, `[[`, "elementary")
  )
  
  
  
  
  
  
  
  
  
  
  library(dplyr)
  library(ggplot2)
  
  alpha <- 0.05
  M <- nr_of_iters
  
  # ------------------------------------------------------------
  # 1. Combine the four procedures
  # ------------------------------------------------------------
  
  procedure_labels <- c(
    Asymptotic_Global     = "Global",
    Finite_Global         = "Global",
    Asymptotic_Holm       = "Holm",
    Finite_Holm           = "Holm",
    Asymptotic_SingleStep = "Single-step",
    Finite_SingleStep     = "Single-step",
    Asymptotic_Tukey      = "Tukey",
    Finite_Tukey          = "Tukey"
  )
  
  core_runs <- bind_rows(test_store) |>
    filter(
      hypothesis == "treatment",
      procedure %in% c(
        "Asymptotic_Global", "Finite_Global",
        "Asymptotic_Holm", "Finite_Holm"
      )
    ) |>
    dplyr::select(N, run, model, hypothesis, procedure, reject)
  
  joint_runs <- joint_test_runs |>
    filter(
      hypothesis == "treatment",
      procedure %in% c(
        "Asymptotic_SingleStep", "Finite_SingleStep",
        "Asymptotic_Tukey", "Finite_Tukey"
      )
    ) |>
    dplyr::select(N, run, model, hypothesis, procedure, reject)
  
  plot_runs <- bind_rows(core_runs, joint_runs)
  
  # Prevent accidental double counting
  stopifnot(
    !anyDuplicated(
      plot_runs[c("N", "run", "model", "hypothesis", "procedure")]
    )
  )
  
  # Check that all eight procedure/reference combinations are present
  missing_procedures <- setdiff(
    names(procedure_labels),
    unique(plot_runs$procedure)
  )
  
  if (length(missing_procedures)) {
    stop(
      "Missing procedures: ",
      paste(missing_procedures, collapse = ", "),
      ". Check the labels in joint_test_runs$procedure."
    )
  }
  
  # ------------------------------------------------------------
  # 2. Summarise rejection rates and failures
  # ------------------------------------------------------------
  
  plot_summary <- plot_runs |>
    group_by(N, model, hypothesis, procedure) |>
    summarise(
      n_stored = n(),
      n_valid = sum(!is.na(reject)),
      n_failed = sum(is.na(reject)),
      n_rejected = sum(reject, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      n_missing = M - n_stored,
      rejection_rate = n_rejected / na_if(n_valid, 0L),
      MCSE = sqrt(
        rejection_rate * (1 - rejection_rate) /
          na_if(n_valid, 0L)
      ),
      adjustment = unname(procedure_labels[procedure]),
      reference = if_else(
        grepl("^Asymptotic_", procedure),
        "Asymptotic",
        "Finite approximation"
      )
    )
  
  stopifnot(all(plot_summary$n_missing >= 0))
  
  print(
    plot_summary |>
      dplyr::select(
        N, model, procedure, n_valid, n_failed,
        n_missing, rejection_rate
      ),
    n = Inf
  )
  
  # ------------------------------------------------------------
  # 3. Optional Lachin overlay
  #    The same benchmark is repeated across panels for comparison.
  # ------------------------------------------------------------
  
  benchmark_summary <- NULL
  
  if (exists("lachin_store")) {
    lachin_runs <- bind_rows(lachin_store) |>
      filter(hypothesis == "treatment")
    
    stopifnot(
      !anyDuplicated(lachin_runs[c("N", "run")])
    )
    
    lachin_summary <- lachin_runs |>
      group_by(N) |>
      summarise(
        n_valid = sum(!is.na(reject)),
        n_failed = sum(is.na(reject)),
        rejection_rate =
          sum(reject, na.rm = TRUE) /
          na_if(sum(!is.na(reject)), 0L),
        .groups = "drop"
      )
    
    print(lachin_summary, n = Inf)
    
    panels <- expand.grid(
      adjustment = c("Global", "Holm", "Single-step", "Tukey"),
      reference = c("Asymptotic", "Finite approximation"),
      stringsAsFactors = FALSE
    )
    
    benchmark_summary <- bind_rows(
      lapply(seq_len(nrow(panels)), function(i) {
        lachin_summary |>
          mutate(
            model = "Lachin",
            adjustment = panels$adjustment[i],
            reference = panels$reference[i]
          )
      })
    )
  }
  
  plot_data <- bind_rows(plot_summary, benchmark_summary) |>
    mutate(
      adjustment = factor(
        adjustment,
        levels = c("Global", "Holm", "Single-step", "Tukey")
      ),
      reference = factor(
        reference,
        levels = c("Asymptotic", "Finite approximation")
      ),
      model = factor(
        model,
        levels = c("GEE", "PGEE", "PIM", "PIM_adj", "Lachin")
      )
    )
  
  # ------------------------------------------------------------
  # 4. Nominal Monte Carlo bounds
  # ------------------------------------------------------------
  
  mc_lower <- max(
    0, alpha - qnorm(0.975) * sqrt(alpha * (1 - alpha) / M)
  )
  
  mc_upper <- min(
    1, alpha + qnorm(0.975) * sqrt(alpha * (1 - alpha) / M)
  )
  
  # ------------------------------------------------------------
  # 5. Plot
  # ------------------------------------------------------------
  
  # Use fewer axis labels to avoid overlap around N = 16–26.
  # Every simulated N still appears as a point.
  N_breaks <- sort(unique(plot_data$N))
  N_breaks <- N_breaks[
    N_breaks %in% c(16, 20, 26, 30, 40, 50, 60)
  ]
  
  p_type1 <- ggplot(
    plot_data,
    aes(
      x = N,
      y = rejection_rate,
      shape = model,
      linetype = model,
      group = model
    )
  ) +
    geom_hline(
      yintercept = alpha,
      colour = "grey45",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_hline(
      yintercept = c(mc_lower, mc_upper),
      colour = "grey65",
      linetype = "dotted",
      linewidth = 0.5
    ) +
    geom_line(linewidth = 0.65, colour = "black", na.rm = TRUE) +
    geom_point(size = 2.6, colour = "black", na.rm = TRUE) +
    facet_grid(reference ~ adjustment) +
    scale_shape_manual(
      values = c(
        GEE = 16, PGEE = 17, PIM = 1,
        PIM_adj = 2, Lachin = 8
      )
    ) +
    scale_linetype_manual(
      values = c(
        GEE = "solid",
        PGEE = "longdash",
        PIM = "dotted",
        PIM_adj = "dotdash",
        Lachin = "twodash"
      )
    ) +
    scale_x_continuous(breaks = N_breaks) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      breaks = seq(0, 1, by = 0.05)
    ) +
    labs(
      title = "Treatment-effect constancy",
      x = "Total number of subjects",
      y = "Empirical Type I error",
      shape = "Method",
      linetype = "Method",
      caption = paste0(
        "Dashed: nominal 5%; dotted: pointwise 95% Monte Carlo bounds",
        " for M = ", format(M, big.mark = ","), ".",
        "\nRates exclude failed runs. Lachin is the same benchmark",
        " repeated across panels."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.caption = element_text(hjust = 0, size = 9),
      legend.key.width = grid::unit(1.1, "cm")
    )
  
  print(p_type1)
  
  ggsave(
    "TypeI_treatment_global_holm_single_step_tukey.pdf",
    plot = p_type1,
    width = 14,
    height = 7
  )
  
  
  
  
  
  library(dplyr)
  library(ggplot2)
  
  time_type1 <- test_runs |>
    filter(hypothesis == "time") |>
    mutate(
      valid = status == "ok" & !is.na(reject),
      reference = case_when(
        grepl("^Asymptotic_", procedure) ~ "Asymptotic",
        grepl("^Finite_", procedure) ~ "Finite approximation"
      ),
      test = case_when(
        grepl("_Global$", procedure) ~ "Global",
        grepl("_Holm$", procedure) ~ "Holm"
      )
    ) |>
    filter(!is.na(reference), !is.na(test)) |>
    group_by(N, model, reference, test) |>
    summarise(
      n_total = n(),
      n_valid = sum(valid),
      n_failed = sum(!valid),
      n_rejected = sum(reject[valid]),
      .groups = "drop"
    ) |>
    mutate(
      rejection_rate = ifelse(
        n_valid > 0, n_rejected / n_valid, NA_real_
      ),
      MCSE = sqrt(
        rejection_rate * (1 - rejection_rate) /
          ifelse(n_valid > 0, n_valid, NA_real_)
      ),
      reference = factor(
        reference,
        levels = c("Asymptotic", "Finite approximation")
      ),
      test = factor(test, levels = c("Global", "Holm")),
      model = factor(
        model, levels = c("GEE", "PGEE", "PIM", "PIM_adj")
      )
    )
  
  
  M <- 10000
  alpha <- 0.05
  
  mc_band <- alpha + c(-1, 1) *
    qnorm(0.975) * sqrt(alpha * (1 - alpha) / M)
  
  ggplot(
    time_type1,
    aes(
      x = N, y = rejection_rate,
      shape = model, linetype = model, group = model
    )
  ) +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = mc_band[1], ymax = mc_band[2],
      fill = "grey85", alpha = 0.6
    ) +
    geom_hline(yintercept = alpha, linetype = "dashed") +
    geom_line() +
    geom_point(size = 2.5) +
    facet_grid(reference ~ test) +
    scale_shape_manual(
      values = c(GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2),
      labels = c("GEE+FW", "PGEE+FW", "PIM", "PIM_adj")
    ) +
    scale_linetype_manual(
      values = c(
        GEE = "solid",
        PGEE = "longdash",
        PIM = "dotted",
        PIM_adj = "dotdash"
      ),
      labels = c("GEE+FW", "PGEE+FW", "PIM", "PIM_adj")
    ) +
    scale_x_continuous(breaks = sort(unique(time_type1$N))) +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(
      title = "Equal temporal effects between groups",
      x = "Total number of subjects",
      y = "Empirical Type I error",
      shape = "Method",
      linetype = "Method",
      caption = paste(
        "Dashed line: nominal 5%.",
        "Shaded band: pointwise 95% Monte Carlo reference range",
        "under 5% rejection with 10,000 valid replications.",
        "Rates exclude failed tests."
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  
  
  
  
  
  library(dplyr)
  
  stored_models <- Filter(Negate(is.null), model_store)
  
  # Ordering in your fitted formula: treatment-group time effect first
  parameter_names <- c(
    "beta_T1", "beta_T0",
    "beta_A1", "beta_A2", "beta_A3", "beta_A4"
  )
  
  metadata <- bind_rows(lapply(stored_models, function(r) {
    data.frame(N = r$N, run = r$run, model = r$model)
  }))
  
  stopifnot(
    !anyDuplicated(metadata[c("N", "run", "model")])
  )
  
  groups <- split(
    seq_along(stored_models),
    interaction(metadata$N, metadata$model, drop = TRUE)
  )
  
  parameter_covariance_summary <- bind_rows(
    lapply(groups, function(idx) {
      
      records <- stored_models[idx]
      p <- length(parameter_names)
      
      extracted <- lapply(records, function(r) {
        b <- as.numeric(r$coefficients)
        V <- as.matrix(r$vcov)
        
        stopifnot(
          length(b) == p,
          identical(dim(V), c(p, p))
        )
        
        list(
          b = b,
          V = (V + t(V)) / 2
        )
      })
      
      # Same replications for empirical and estimated covariance
      valid <- vapply(extracted, function(z) {
        all(is.finite(z$b)) &&
          all(is.finite(z$V)) &&
          all(diag(z$V) > 0)
      }, logical(1))
      
      fits <- extracted[valid]
      n_valid <- length(fits)
      
      empirical <- estimated <- matrix(NA_real_, p, p)
      mean_b <- rep(NA_real_, p)
      
      if (n_valid > 0L) {
        B <- do.call(rbind, lapply(fits, function(z) z$b))
        
        mean_b <- colMeans(B)
        
        estimated <- Reduce(
          `+`, lapply(fits, function(z) z$V)
        ) / n_valid
        
        if (n_valid >= 2L) {
          empirical <- cov(B)
        }
      }
      
      empirical_cor <- estimated_cor <- matrix(NA_real_, p, p)
      
      if (all(is.finite(empirical)) &&
          all(diag(empirical) > 0)) {
        empirical_cor <- cov2cor(empirical)
      }
      
      if (all(is.finite(estimated)) &&
          all(diag(estimated) > 0)) {
        estimated_cor <- cov2cor(estimated)
      }
      
      ij <- which(
        upper.tri(estimated, diag = TRUE),
        arr.ind = TRUE
      )
      
      data.frame(
        N = records[[1]]$N,
        model = records[[1]]$model,
        parameter_1 = parameter_names[ij[, 1]],
        parameter_2 = parameter_names[ij[, 2]],
        n_total = length(records),
        n_valid = n_valid,
        n_excluded = sum(!valid),
        mean_estimate_1 = mean_b[ij[, 1]],
        mean_estimate_2 = mean_b[ij[, 2]],
        empirical_covariance = empirical[ij],
        average_estimated_covariance = estimated[ij],
        covariance_difference = estimated[ij] - empirical[ij],
        empirical_correlation = empirical_cor[ij],
        estimated_correlation = estimated_cor[ij],
        correlation_difference =
          estimated_cor[ij] - empirical_cor[ij]
      )
    })
  ) |>
    arrange(N, model, parameter_1, parameter_2)
  
  # Parameter variance comparison
  parameter_variance_summary <- parameter_covariance_summary |>
    filter(parameter_1 == parameter_2) |>
    mutate(
      empirical_SD = sqrt(empirical_covariance),
      estimated_SD = sqrt(average_estimated_covariance),
      variance_ratio =
        empirical_covariance / average_estimated_covariance
    )
  
  # Off-diagonal parameter covariance comparison
  parameter_covariance_offdiag <- parameter_covariance_summary |>
    filter(parameter_1 != parameter_2)
  
  print(parameter_variance_summary, n = Inf)
  print(parameter_covariance_offdiag, n = Inf)
  
  # Treatment-parameter covariances for the original PIM
  parameter_covariance_offdiag |>
    filter(
      model == "PIM",
      grepl("^beta_A", parameter_1),
      grepl("^beta_A", parameter_2)
    ) |>
    print(n = Inf)
  
  library(ggplot2)
  library(tidyr)
  
  parameter_covariance_offdiag |>
    filter(
      grepl("^beta_T", parameter_1),
      grepl("^beta_T", parameter_2)
    ) |>
    mutate(
      parameter_pair = paste(parameter_1, parameter_2, sep = " / ")
    ) |>
    pivot_longer(
      cols = c(empirical_correlation, estimated_correlation),
      names_to = "source",
      values_to = "correlation"
    ) |>
    mutate(
      source = recode(
        source,
        empirical_correlation = "Empirical",
        estimated_correlation = "Estimated"
      )
    ) |>
    ggplot(aes(
      x = N, y = correlation,
      colour = model, linetype = source
    )) +
    geom_line() +
    geom_point() +
    facet_wrap(~ parameter_pair) +
    labs(
      x = "Total number of subjects",
      y = "Parameter correlation",
      colour = "Model",
      linetype = "Covariance source"
    ) +
    theme_bw()
  
  
  
  
## Additional runs with GEE only
  
  # ============================================================
  # 1. Required existing PO objects and functions
  # ============================================================
  
  required_po_sa <- c(
    "po_template", "sample_sizes", "nr_of_iters", "nr_of_times",
    "geessbin_adjusted", "test_contrasts", "test_stored_fit",
    "L_trt", "L_time", "L_tukey",
    "model_store", "test_store", "contrast_store",
    "joint_test_runs", "joint_elementary_runs"
  )
  
  missing_po_sa <- required_po_sa[
    !vapply(required_po_sa, exists, logical(1), inherits = TRUE)
  ]
  
  if (length(missing_po_sa)) {
    stop(
      "Missing objects: ",
      paste(missing_po_sa, collapse = ", ")
    )
  }
  
  # Reuse the ORIGINAL calibrated distribution.
  # Do not regenerate po_template.
  stopifnot(
    ncol(po_template$trajectories) == nr_of_times,
    length(po_template$joint_probabilities$control) ==
      nrow(po_template$trajectories),
    length(po_template$joint_probabilities$treatment) ==
      nrow(po_template$trajectories)
  )
  
  # Assumes the original PO simulation used R's default RNG.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  
  sample_sizes <- c(16, 18, 20, 24,26, 30, 40, 50, 60)
  
  n_po_sa <- length(sample_sizes) * nr_of_iters
  
  # New stores: original stores are preserved
  model_store_PO_GEE_SA <- vector("list", n_po_sa)
  test_store_PO_GEE_SA <- vector("list", 2L * n_po_sa)
  contrast_store_PO_GEE_SA <- vector("list", 2L * n_po_sa)
  joint_store_PO_GEE_SA <- vector("list", n_po_sa)
  
  model_index_po_sa <- 0L
  test_index_po_sa <- 0L
  
  contrast_families_po_sa <- list(
    treatment = L_trt,
    time = L_time
  )
  
  # Local copy: does not overwrite your global 'records'
  joint_env_po_sa <- new.env(
    parent = environment(test_stored_fit)
  )
  
  test_joint_po_sa <- test_stored_fit
  environment(test_joint_po_sa) <- joint_env_po_sa
  
  # ============================================================
  # 2. Regenerate PO datasets and fit ONLY ordinary GEE + SA
  # ============================================================

  for (N in sample_sizes) {
    
    stopifnot(N %% 2L == 0L, nr_of_times >= 2L)
    
    design_po_sa <- data.frame(
      id = rep(seq_len(N), each = nr_of_times),
      time = rep(seq_len(nr_of_times), N),
      group = rep(
        rep(c(0, 1), each = N / 2),
        each = nr_of_times
      )
    )
    
    # ----------------------------------------------------------
    # Fixed comparisons and predictors: once per sample size
    # ----------------------------------------------------------
    
    between_po_sa <- do.call(
      rbind,
      lapply(seq_len(nr_of_times), function(t) {
        expand.grid(
          Var1 = which(
            design_po_sa$group == 0 &
              design_po_sa$time == t
          ),
          Var2 = which(
            design_po_sa$group == 1 &
              design_po_sa$time == t
          )
        )
      })
    )
    
    between_po_sa <- between_po_sa[
      order(between_po_sa$Var2), , drop = FALSE
    ]
    
    within_po_sa <- do.call(
      rbind,
      lapply(seq_len(N), function(s) {
        idx <- which(design_po_sa$id == s)
        pairs <- t(combn(idx, 2L))
        
        data.frame(
          Var1 = pairs[, 1],
          Var2 = pairs[, 2]
        )
      })
    )
    
    within_po_sa <- within_po_sa[
      order(within_po_sa$Var2), , drop = FALSE
    ]
    
    comparisons_po_sa <- rbind(
      between_po_sa, within_po_sa
    )
    
    left_po_sa <- comparisons_po_sa$Var1
    right_po_sa <- comparisons_po_sa$Var2
    
    dt_po_sa <- design_po_sa$time[right_po_sa] -
      design_po_sa$time[left_po_sa]
    
    gL_po_sa <- design_po_sa$group[left_po_sa]
    gR_po_sa <- design_po_sa$group[right_po_sa]
    
    # Original parameter order:
    # treated time, control time, treatment effects by visit
    x_po_sa <- cbind(
      dt_po_sa * gR_po_sa * gL_po_sa,
      dt_po_sa * (1 - gR_po_sa) * (1 - gL_po_sa),
      sapply(seq_len(nr_of_times), function(t) {
        (gR_po_sa - gL_po_sa) *
          (design_po_sa$time[right_po_sa] == t) *
          (design_po_sa$time[left_po_sa] == t)
      })
    )
    
    stopifnot(
      ncol(x_po_sa) == ncol(L_trt),
      ncol(x_po_sa) == ncol(L_time),
      ncol(x_po_sa) == ncol(L_tukey)
    )
    
    base_po_sa <- data.frame(
      y = numeric(nrow(x_po_sa)),
      x_po_sa
    )
    
    names(base_po_sa) <- c(
      "y", paste0("x", seq_len(ncol(x_po_sa)))
    )
    
    base_po_sa$C1 <- design_po_sa$id[left_po_sa]
    base_po_sa$C2 <- design_po_sa$id[right_po_sa]
    base_po_sa$C3 <- paste(
      base_po_sa$C1, base_po_sa$C2, sep = "_"
    )
    
    orders_po_sa <- list(
      order(base_po_sa$C1),
      order(base_po_sa$C2),
      order(base_po_sa$C3)
    )
    
    for (l in seq_len(nr_of_iters)) {
      
      if (l == 1L || l %% 100L == 0L) {
        cat("PO GEE without FW: N =", N, "; run =", l, "\n")
      }
      
      # --------------------------------------------------------
      # EXACT original sampling sequence
      # --------------------------------------------------------
      
      set.seed(l)
      
      index0_po_sa <- sample.int(
        nrow(po_template$trajectories),
        N / 2,
        replace = TRUE,
        prob = po_template$joint_probabilities$control
      )
      
      index1_po_sa <- sample.int(
        nrow(po_template$trajectories),
        N / 2,
        replace = TRUE,
        prob = po_template$joint_probabilities$treatment
      )
      
      design_po_sa$Y <- as.vector(t(rbind(
        po_template$trajectories[
          index0_po_sa, , drop = FALSE
        ],
        po_template$trajectories[
          index1_po_sa, , drop = FALSE
        ]
      )))
      
      dat_po_sa <- base_po_sa
      
      dat_po_sa$y <-
        as.numeric(
          design_po_sa$Y[left_po_sa] <
            design_po_sa$Y[right_po_sa]
        ) +
        0.5 * as.numeric(
          design_po_sa$Y[left_po_sa] ==
            design_po_sa$Y[right_po_sa]
        )
      
      # --------------------------------------------------------
      # Three fits: ordinary GEE, unadjusted sandwich
      # --------------------------------------------------------
      
      warnings_po_sa <- character()
      
      fit_po_sa <- tryCatch(
        withCallingHandlers({
          
          d1_po_sa <- dat_po_sa[orders_po_sa[[1]], ]
          
          m1_po_sa <- geessbin_adjusted(
            y ~ . - 1 - C1 - C2 - C3,
            data = d1_po_sa,
            id = C1,
            corstr = "independence",
            beta.method = "GEE",
            SE.method = "SA"
          )
          
          d2_po_sa <- dat_po_sa[orders_po_sa[[2]], ]
          
          m2_po_sa <- geessbin_adjusted(
            y ~ . - 1 - C1 - C2 - C3,
            data = d2_po_sa,
            id = C2,
            corstr = "independence",
            beta.method = "GEE",
            SE.method = "SA",
            b = as.numeric(coef(m1_po_sa))
          )
          
          d3_po_sa <- dat_po_sa[orders_po_sa[[3]], ]
          
          m3_po_sa <- geessbin_adjusted(
            y ~ . - 1 - C1 - C2 - C3,
            data = d3_po_sa,
            id = C3,
            corstr = "independence",
            beta.method = "GEE",
            SE.method = "SA",
            b = as.numeric(coef(m1_po_sa))
          )
          
          b_po_sa <- coef(m1_po_sa)
          
          V_po_sa <- m1_po_sa$covb +
            m2_po_sa$covb - m3_po_sa$covb
          
          V_po_sa <- (V_po_sa + t(V_po_sa)) / 2
          
          stopifnot(
            length(b_po_sa) == ncol(x_po_sa),
            all(dim(V_po_sa) == ncol(x_po_sa)),
            all(is.finite(b_po_sa)),
            all(is.finite(V_po_sa))
          )
          
          B_po_sa <- rbind(
            coef(m1_po_sa),
            coef(m2_po_sa),
            coef(m3_po_sa)
          )
          
          list(
            b = b_po_sa,
            V = V_po_sa,
            max_coefficient_difference = max(
              apply(B_po_sa, 2, function(z) diff(range(z)))
            ),
            min_eigenvalue = min(
              eigen(
                V_po_sa,
                symmetric = TRUE,
                only.values = TRUE
              )$values
            ),
            status = "ok",
            error = NA_character_
          )
          
        }, warning = function(w) {
          warnings_po_sa <<- c(
            warnings_po_sa, conditionMessage(w)
          )
          invokeRestart("muffleWarning")
        }),
        
        error = function(e) {
          list(
            b = rep(NA_real_, ncol(x_po_sa)),
            V = matrix(
              NA_real_, ncol(x_po_sa), ncol(x_po_sa)
            ),
            max_coefficient_difference = NA_real_,
            min_eigenvalue = NA_real_,
            status = "failed",
            error = conditionMessage(e)
          )
        }
      )
      
      # --------------------------------------------------------
      # Store new model and diagnostics
      # --------------------------------------------------------
      
      model_index_po_sa <- model_index_po_sa + 1L
      
      model_store_PO_GEE_SA[[model_index_po_sa]] <- list(
        N = N,
        run = l,
        model = "GEE_SA",
        coefficients = fit_po_sa$b,
        vcov = fit_po_sa$V,
        status = fit_po_sa$status,
        error = fit_po_sa$error,
        warning = paste(
          unique(warnings_po_sa), collapse = "; "
        ),
        max_coefficient_difference =
          fit_po_sa$max_coefficient_difference,
        min_eigenvalue = fit_po_sa$min_eigenvalue
      )
      
      # --------------------------------------------------------
      # Global and Holm: treatment and time
      # --------------------------------------------------------
      
      for (hypothesis_po_sa in names(contrast_families_po_sa)) {
        
        test_index_po_sa <- test_index_po_sa + 1L
        L_po_sa <- contrast_families_po_sa[[hypothesis_po_sa]]
        
        result_po_sa <- tryCatch({
          
          if (fit_po_sa$status != "ok") {
            stop(fit_po_sa$error)
          }
          
          test_contrasts(
            b = fit_po_sa$b,
            V = fit_po_sa$V,
            L = L_po_sa,
            nu = N - 1,
            alpha = 0.05
          )$global
          
        }, error = function(e) {
          data.frame(
            procedure = c(
              "Asymptotic_Global", "Asymptotic_Holm",
              "Finite_Global", "Finite_Holm"
            ),
            statistic = NA_real_,
            df1 = NA_real_,
            df2 = NA_real_,
            p_global = NA_real_,
            reject = NA,
            status = "failed",
            error = conditionMessage(e)
          )
        })
        
        test_store_PO_GEE_SA[[test_index_po_sa]] <- data.frame(
          N = N,
          run = l,
          model = "GEE_SA",
          hypothesis = hypothesis_po_sa,
          result_po_sa,
          row.names = NULL
        )
        
        estimate_po_sa <- drop(L_po_sa %*% fit_po_sa$b)
        
        variance_po_sa <- diag(
          L_po_sa %*% fit_po_sa$V %*% t(L_po_sa)
        )
        
        contrast_store_PO_GEE_SA[[test_index_po_sa]] <-
          data.frame(
            N = N,
            run = l,
            model = "GEE_SA",
            hypothesis = hypothesis_po_sa,
            contrast = seq_len(nrow(L_po_sa)),
            estimate = estimate_po_sa,
            variance = variance_po_sa,
            status = ifelse(
              is.finite(estimate_po_sa) &
                is.finite(variance_po_sa) &
                variance_po_sa > 0,
              "ok",
              "invalid estimate or variance"
            ),
            error = fit_po_sa$error,
            row.names = NULL
          )
      }
      
      # --------------------------------------------------------
      # Single-step and Tukey: treatment contrasts
      # --------------------------------------------------------
      
      joint_env_po_sa$records <- list(
        model_store_PO_GEE_SA[[model_index_po_sa]]
      )
      
      joint_store_PO_GEE_SA[[model_index_po_sa]] <-
        test_joint_po_sa(1L)
    }
  }
  
  # ============================================================
  # 3. Extract new results
  # ============================================================
  
  test_runs_PO_GEE_SA <- bind_rows(test_store_PO_GEE_SA)
  
  contrast_runs_PO_GEE_SA <- bind_rows(
    contrast_store_PO_GEE_SA
  )
  
  joint_test_runs_PO_GEE_SA <- bind_rows(
    lapply(joint_store_PO_GEE_SA, `[[`, "family")
  )
  
  joint_elementary_runs_PO_GEE_SA <- bind_rows(
    lapply(joint_store_PO_GEE_SA, `[[`, "elementary")
  )
  
  # ============================================================
  # 4. Combined PO objects: originals remain unchanged
  # ============================================================
  
  model_store_PO_all <- c(
    Filter(
      function(z) {
        !is.null(z) && !identical(z$model, "GEE_SA")
      },
      model_store
    ),
    model_store_PO_GEE_SA
  )
  
  test_runs_PO_all <- bind_rows(test_store) |>
    filter(as.character(model) != "GEE_SA") |>
    bind_rows(test_runs_PO_GEE_SA)
  
  contrast_runs_PO_all <- bind_rows(contrast_store) |>
    filter(as.character(model) != "GEE_SA") |>
    bind_rows(contrast_runs_PO_GEE_SA)
  
  joint_test_runs_PO_all <- joint_test_runs |>
    filter(as.character(model) != "GEE_SA") |>
    bind_rows(joint_test_runs_PO_GEE_SA)
  
  joint_elementary_runs_PO_all <- joint_elementary_runs |>
    filter(as.character(model) != "GEE_SA") |>
    bind_rows(joint_elementary_runs_PO_GEE_SA)
  
  stopifnot(
    !anyDuplicated(
      test_runs_PO_all[
        c("N", "run", "model", "hypothesis", "procedure")
      ]
    ),
    !anyDuplicated(
      joint_test_runs_PO_all[
        c("N", "run", "model", "hypothesis", "procedure")
      ]
    )
  )
  
  # ============================================================
  # 5. Summary of the additional method
  # ============================================================
  
  summary_PO_GEE_SA <- bind_rows(
    test_runs_PO_GEE_SA |>
      select(N, run, model, hypothesis, procedure, reject),
    joint_test_runs_PO_GEE_SA |>
      select(N, run, model, hypothesis, procedure, reject)
  ) |>
    group_by(N, hypothesis, procedure) |>
    summarise(
      n_total = n(),
      n_valid = sum(!is.na(reject)),
      n_failed = sum(is.na(reject)),
      n_rejected = sum(reject, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      rejection_rate = n_rejected / na_if(n_valid, 0L)
    )
  
  print(summary_PO_GEE_SA, n = Inf)
  
  diagnostics_PO_GEE_SA <- bind_rows(
    lapply(model_store_PO_GEE_SA, function(z) {
      data.frame(
        N = z$N,
        run = z$run,
        status = z$status,
        warning = z$warning,
        max_coefficient_difference =
          z$max_coefficient_difference,
        min_eigenvalue = z$min_eigenvalue
      )
    })
  )
  
  diagnostics_PO_GEE_SA |>
    group_by(N) |>
    summarise(
      n_failed = sum(status == "failed"),
      n_with_warnings = sum(nzchar(warning)),
      n_indefinite = sum(
        min_eigenvalue < -1e-10, na.rm = TRUE
      ),
      n_coefficients_differ = sum(
        max_coefficient_difference > 1e-5,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    print(n = Inf)
  
  
  

  
  library(dplyr)
  
  # PO:
  exclusion_runs <- bind_rows(
    test_runs_PO_all |>
      select(N, run, model, hypothesis, procedure, reject),
    joint_test_runs_PO_all |>
      select(N, run, model, hypothesis, procedure, reject)
  )
  
  # For Gaussian simulations, use instead:
  # exclusion_runs <- bind_rows(
  #   test_runs_all |>
  #     select(N, run, model, hypothesis, procedure, reject),
  #   joint_test_runs_all |>
  #     select(N, run, model, hypothesis, procedure, reject)
  # )
  
  stopifnot(
    !anyDuplicated(
      exclusion_runs[
        c("N", "run", "model", "hypothesis", "procedure")
      ]
    )
  )
  
  exclusion_summary <- exclusion_runs |>
    group_by(N, model, hypothesis, procedure) |>
    summarise(
      n_total = n(),
      n_valid = sum(!is.na(reject)),
      n_excluded = sum(is.na(reject)),
      .groups = "drop"
    ) |>
    mutate(
      excluded_percent = 100 * n_excluded / n_total
    ) |>
    filter(n_excluded > 0) |>
    arrange(model, hypothesis, procedure, N)
  
  print(exclusion_summary, n = Inf)
  
  
  denominator_check <- exclusion_runs |>
    group_by(N, model, hypothesis, procedure) |>
    summarise(
      n_stored = n(),
      n_excluded = sum(is.na(reject)),
      denominator = sum(!is.na(reject)),
      n_rejected = sum(reject, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      rejection_rate = n_rejected / na_if(denominator, 0L)
    ) |>
    arrange(hypothesis, procedure, model, N)
  
  print(denominator_check, n = Inf)
  
  # Smallest denominator across all combinations
  min(denominator_check$denominator)
  
  # Print only combinations using fewer than 10,000 runs
  denominator_check |>
    filter(denominator < 10000& model%in%c("PIM_adj","GEE_SA"))  |>
    print(n = Inf)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  library(dplyr)
  library(ggplot2)
  
  # ============================================================
  # 1. Choose simulation
  # ============================================================
  
  simulation <- "PO"  # Change to "RI" for Gaussian simulations
  
  if (simulation == "PO") {
    
    # Rebuild from stores to include any additional sample sizes
    core_new <- bind_rows(test_store_PO_GEE_SA)
    
    joint_new <- bind_rows(
      lapply(joint_store_PO_GEE_SA, `[[`, "family")
    )
    
  } else if (simulation == "RI") {
    
    core_new <- bind_rows(test_store_GEE_SA)
    
    joint_new <- bind_rows(
      lapply(joint_store_GEE_SA, `[[`, "family")
    )
    
  } else {
    stop("simulation must be 'PO' or 'RI'.")
  }
  
  # Original objects must belong to the chosen simulation.
  # Preserve them; create new plotting objects.
  core_plot_runs <- bind_rows(test_store) |>
    filter(as.character(model) != "GEE_SA") |>
    bind_rows(core_new) |>
    mutate(model = as.character(model))
  
  joint_plot_runs <- joint_test_runs |>
    filter(as.character(model) != "GEE_SA") |>
    bind_rows(joint_new) |>
    mutate(model = as.character(model))
  
  procedure_labels <- c(
    Asymptotic_Global     = "Global",
    Finite_Global         = "Global",
    Asymptotic_Holm       = "Holm",
    Finite_Holm           = "Holm",
    Asymptotic_SingleStep = "Single-step",
    Finite_SingleStep     = "Single-step",
    Asymptotic_Tukey      = "Tukey",
    Finite_Tukey          = "Tukey"
  )
  
  all_plot_runs <- bind_rows(
    core_plot_runs |>
      filter(procedure %in% names(procedure_labels)) |>
      select(N, run, model, hypothesis, procedure, reject),
    
    joint_plot_runs |>
      filter(procedure %in% names(procedure_labels)) |>
      select(N, run, model, hypothesis, procedure, reject)
  )
  
  if (anyDuplicated(
    all_plot_runs[
      c("N", "run", "model", "hypothesis", "procedure")
    ]
  )) {
    stop("Duplicate simulation results detected.")
  }
  
  # ============================================================
  # 2. Rejection rates and failures
  # ============================================================
  
  alpha <- 0.05
  M <- nr_of_iters
  
  type1_summary <- all_plot_runs |>
    group_by(N, model, hypothesis, procedure) |>
    summarise(
      n_stored = n(),
      n_valid = sum(!is.na(reject)),
      n_failed = sum(is.na(reject)),
      n_rejected = sum(reject, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      rejection_rate = n_rejected / na_if(n_valid, 0L),
      adjustment = unname(procedure_labels[procedure]),
      reference = if_else(
        grepl("^Asymptotic_", procedure),
        "Asymptotic",
        "Finite approximation"
      )
    )
  
  print(
    type1_summary |>
      select(
        N, model, hypothesis, procedure,
        n_valid, n_failed, rejection_rate
      ),
    n = Inf
  )
  
  # ============================================================
  # 3. Lachin benchmark
  # ============================================================
  
  lachin_plot_summary <- NULL
  
  if (exists("lachin_store")) {
    
    lachin_plot_runs <- bind_rows(lachin_store)
    
    if (nrow(lachin_plot_runs) > 0L) {
      
      if (anyDuplicated(
        lachin_plot_runs[c("N", "run")]
      )) {
        stop("Duplicate Lachin results detected.")
      }
      
      lachin_plot_summary <- lachin_plot_runs |>
        group_by(N) |>
        summarise(
          n_stored = n(),
          n_valid = sum(!is.na(reject)),
          n_failed = sum(is.na(reject)),
          n_rejected = sum(reject, na.rm = TRUE),
          .groups = "drop"
        ) |>
        mutate(
          rejection_rate =
            n_rejected / na_if(n_valid, 0L),
          model = "Lachin"
        )
    }
  }
  
  # ============================================================
  # 4. Common plotting settings
  # ============================================================
  
  method_levels <- c(
    "GEE_SA", "GEE", "PGEE", "PIM", "PIM_adj", "Lachin"
  )
  
  method_labels <- c(
    GEE_SA = "GEE",
    GEE = "GEE+FW",
    PGEE = "PGEE+FW",
    PIM = "PIM",
    PIM_adj = "PIM_adj",
    Lachin = "Lachin"
  )
  
  method_shapes <- c(
    GEE_SA = 15,
    GEE = 16,
    PGEE = 17,
    PIM = 1,
    PIM_adj = 2,
    Lachin = 8
  )
  
  method_lines <- c(
    GEE_SA = "dashed",
    GEE = "solid",
    PGEE = "longdash",
    PIM = "dotted",
    PIM_adj = "dotdash",
    Lachin = "twodash"
  )
  
  # Pointwise 95% reference band under a true rejection rate of 5%.
  # Based on the planned number of replications.
  mc_band <- alpha + c(-1, 1) *
    qnorm(0.975) * sqrt(alpha * (1 - alpha) / M)
  
  mc_band <- pmax(0, pmin(1, mc_band))
  
  make_type1_plot <- function(plot_data, title, caption) {
    
    plot_data <- plot_data |>
      mutate(
        model = factor(model, levels = method_levels),
        reference = factor(
          reference,
          levels = c("Asymptotic", "Finite approximation")
        ),
        adjustment = factor(
          adjustment,
          levels = c("Global", "Holm", "Single-step", "Tukey")
        )
      )
    
    ggplot(
      plot_data,
      aes(
        x = N,
        y = rejection_rate,
        shape = model,
        linetype = model,
        group = model
      )
    ) +
      annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = mc_band[1], ymax = mc_band[2],
        fill = "grey85", alpha = 0.6
      ) +
      geom_hline(
        yintercept = alpha,
        linetype = "dashed",
        colour = "grey45"
      ) +
      geom_line(
        linewidth = 0.6,
        na.rm = TRUE
      ) +
      geom_point(
        size = 2.4,
        na.rm = TRUE
      ) +
      facet_grid(reference ~ adjustment) +
      scale_shape_manual(
        values = method_shapes,
        breaks = method_levels,
        labels = method_labels
      ) +
      scale_linetype_manual(
        values = method_lines,
        breaks = method_levels,
        labels = method_labels
      ) +
      scale_x_continuous(
        breaks = sort(unique(plot_data$N))
      ) +
      scale_y_continuous(
        labels = scales::label_percent(accuracy = 1)
      ) +
      labs(
        title = title,
        x = "Total number of subjects",
        y = "Empirical Type I error",
        shape = "Method",
        linetype = "Method",
        caption = caption
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        legend.key.width = grid::unit(1.2, "cm"),
        axis.text.x = element_text(
          angle = 45, hjust = 1
        ),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(
          hjust = 0, size = 9
        )
      )
  }
  
  # ============================================================
  # 5. Treatment-effect constancy: four procedures
  # ============================================================
  
  treatment_plot_data <- type1_summary |>
    filter(hypothesis == "treatment")
  
  # Repeat the SAME Lachin test across panels for comparison
  if (!is.null(lachin_plot_summary)) {
    
    panels <- expand.grid(
      adjustment = c("Global", "Holm", "Single-step", "Tukey"),
      reference = c("Asymptotic", "Finite approximation"),
      stringsAsFactors = FALSE
    )
    
    lachin_overlay <- bind_rows(
      lapply(seq_len(nrow(panels)), function(i) {
        lachin_plot_summary |>
          mutate(
            adjustment = panels$adjustment[i],
            reference = panels$reference[i]
          )
      })
    )
    
    treatment_plot_data <- bind_rows(
      treatment_plot_data,
      lachin_overlay
    )
  }
  
  p_type1_treatment <- make_type1_plot(
    treatment_plot_data,
    title = paste0(
      simulation, ": treatment-effect constancy"
    ),
    caption = paste0(
      "Dashed line: nominal 5%; shaded band: pointwise 95% ",
      "Monte Carlo reference range for M = ",
      format(M, big.mark = ","), ".\n",
      "Rates exclude failed tests. Lachin, when shown, is the ",
      "same asymptotic benchmark repeated across panels."
    )
  )
  
  print(p_type1_treatment)
  
  # ============================================================
  # 6. Equality of time trends: global and Holm
  # ============================================================
  
  time_plot_data <- type1_summary |>
    filter(
      hypothesis == "time",
      adjustment %in% c("Global", "Holm")
    )
  
  p_type1_time <- make_type1_plot(
    time_plot_data,
    title = paste0(
      simulation, ": equality of temporal effects between groups"
    ),
    caption = paste0(
      "Dashed line: nominal 5%; shaded band: pointwise 95% ",
      "Monte Carlo reference range for M = ",
      format(M, big.mark = ","), ".\n",
      "Rates exclude failed tests. For this single contrast, ",
      "global and Holm results should coincide within each ",
      "reference distribution."
    )
  )
  
  print(p_type1_time)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  
  