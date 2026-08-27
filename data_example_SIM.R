## -------------------------------------------------------------------------
## 0. Packages
## -------------------------------------------------------------------------

needed <- c("dplyr", "tidyr", "ggplot2", "multgee", "multcomp", "geessbin", "MASS", "aod")
to_install <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(tidyr)
library(ggplot2)
library(multgee)
library(multcomp)
library(geessbin)
library(MASS)
library(aod)

has_pim  <- requireNamespace("pim", quietly = TRUE)
has_npar <- requireNamespace("nparLD", quietly = TRUE)
has_lwo  <- requireNamespace("lwo", quietly = TRUE)

## -------------------------------------------------------------------------
## 1. Modified geessbin function allowing pseudo-outcomes in {0, 0.5, 1}
## -------------------------------------------------------------------------

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
  if (!is.numeric(y) || any(!unique(as.vector(y)) %in% c(0, 0.5, 1))) {
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

# Make the function behave as a geessbin object and allow hidden geessbin helpers.
environment(geessbin_adjusted) <- asNamespace("geessbin")

## -------------------------------------------------------------------------
## 2. Helper functions
## -------------------------------------------------------------------------

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

pseudo_score <- function(y_left, y_right, higher_is_better = TRUE) {
  # left = control, right = treatment for between-group pairs.
  # score = 1 means that the right observation wins.
  if (higher_is_better) {
    ifelse(y_right > y_left, 1, ifelse(y_right < y_left, 0, 0.5))
  } else {
    ifelse(y_right < y_left, 1, ifelse(y_right > y_left, 0, 0.5))
  }
}

make_small_arthritis <- function(n_per_group = 15,
                                 seed = 1988,
                                 id_col = "id",
                                 group_col = "trt",
                                 time_col = "time",
                                 outcome_col = "y") {
  data("arthritis", package = "multgee")
  dat <- as.data.frame(multgee::arthritis)
  dat <- na.omit(dat)
  names(dat) <- tolower(names(dat))
  id_col <- tolower(id_col)
  group_col <- tolower(group_col)
  time_col <- tolower(time_col)
  outcome_col <- tolower(outcome_col)

  required <- c(id_col, group_col, time_col, outcome_col)
  miss <- setdiff(required, names(dat))
  if (length(miss) > 0) {
    stop("Required columns not found in multgee::arthritis: ", paste(miss, collapse = ", "),
         "\nAvailable columns are: ", paste(names(dat), collapse = ", "))
  }

  dat <- dat %>%
    mutate(
      id_raw = .data[[id_col]],
      group_raw = .data[[group_col]],
      time_raw = .data[[time_col]],
      y_raw = .data[[outcome_col]]
    )

  group_levels <- sort(unique(dat$group_raw))
  if (length(group_levels) != 2) stop("Expected exactly two treatment groups.")

  dat <- dat %>%
    mutate(
      group = ifelse(group_raw == group_levels[1], 0L, 1L),
      id = as.integer(factor(id_raw)),
      time = as.integer(factor(time_raw, levels = sort(unique(time_raw)))),
      y = as.numeric(y_raw)
    )

  set.seed(seed)
  selected_ids <- dat %>%
    distinct(id, group) %>%
    group_by(group) %>%
    slice_sample(n = n_per_group) %>%
    ungroup()

  dat_small <- dat %>%
    semi_join(selected_ids, by = c("id", "group")) %>%
    mutate(id = as.integer(factor(id))) %>%
    arrange(group, id, time)

  dat_small
}

make_compare_and_dat_gee <- function(dat,
                                     higher_is_better = TRUE,
                                     include_within_subject_pairs = TRUE) {
  dat <- dat %>% arrange(id, time) %>% mutate(row_id = row_number())
  times <- sort(unique(dat$time))
  all_pairs <- list()

  # Between-treatment pairs within each visit: Var1 = control, Var2 = treatment.
  for (tt in times) {
    id.fac <- which(dat$group == 0 & dat$time == tt)
    id.nonfac <- which(dat$group == 1 & dat$time == tt)
    if (length(id.fac) > 0 && length(id.nonfac) > 0) {
      tmp <- expand.grid(Var1 = id.fac, Var2 = id.nonfac)
      tmp$pair_type <- "between"
      all_pairs[[length(all_pairs) + 1]] <- tmp
    }
  }

  # Within-subject pairs over time: Var1 = earlier visit, Var2 = later visit.
  if (include_within_subject_pairs) {
    for (ii in sort(unique(dat$id))) {
      idx <- which(dat$id == ii)
      idx <- idx[order(dat$time[idx])]
      if (length(idx) >= 2) {
        tmp <- t(utils::combn(idx, 2))
        tmp <- data.frame(Var1 = tmp[, 1], Var2 = tmp[, 2])
        tmp$pair_type <- "within"
        all_pairs[[length(all_pairs) + 1]] <- tmp
      }
    }
  }

  compare <- dplyr::bind_rows(all_pairs)
  L <- dat[compare$Var1, ]
  R <- dat[compare$Var2, ]

  y <- pseudo_score(L$y, R$y, higher_is_better = higher_is_better)

  # Model 2-type design:
  #   two within-subject time-trend parameters, one for treated and one for controls;
  #   one treatment effect parameter per visit.
  X <- data.frame(
    trend_treat = (R$time - L$time) * R$group * L$group,
    trend_ctrl  = (R$time - L$time) * (1 - R$group) * (1 - L$group)
  )
  for (tt in times) {
    X[[paste0("trt_visit", tt)]] <- (R$group - L$group) * (R$time == tt) * (L$time == tt)
  }

  # Remove accidental all-zero columns.
  keep <- vapply(X, function(z) any(abs(z) > 0), logical(1))
  X <- X[, keep, drop = FALSE]

  C1 <- dat[compare$Var1, "id"] %>% unlist(use.names = FALSE)
  C2 <- dat[compare$Var2, "id"] %>% unlist(use.names = FALSE)

  dat_GEE <- data.frame(y = y, X, C1 = C1, C2 = C2)
  dat_GEE$C3 <- paste(dat_GEE$C1, dat_GEE$C2, sep = "_")

  list(
    data_long = dat,
    compare = compare,
    dat_GEE = dat_GEE,
    x_names = names(X),
    treatment_terms = grep("^trt_visit", names(X), value = TRUE),
    times = times
  )
}

fit_three_geessbin_FW <- function(dat_GEE) {
  dat1 <- dat_GEE[order(dat_GEE$C1), ]
  mod1 <- geessbin_adjusted(y ~ . - 1 - C1 - C2 - C3,
                            data = dat1,
                            id = C1,
                            corstr = "independence",
                            beta.method = "PGEE",
                            SE.method = "FW")

  dat2 <- dat_GEE[order(dat_GEE$C2), ]
  mod2 <- geessbin_adjusted(y ~ . - 1 - C1 - C2 - C3,
                            data = dat2,
                            id = C2,
                            corstr = "independence",
                            beta.method = "PGEE",
                            SE.method = "FW")

  dat3 <- dat_GEE[order(dat_GEE$C3), ]
  mod3 <- geessbin_adjusted(y ~ . - 1 - C1 - C2 - C3,
                            data = dat3,
                            id = C3,
                            corstr = "independence",
                            beta.method = "PGEE",
                            SE.method = "FW")

  # Inclusion--exclusion covariance estimator. This is the estimator used in
  # the manuscript; no symmetrisation is applied here.
  V <- mod1$covb + mod2$covb - mod3$covb

  # Follow the simulation code: the point estimate is the average of the three
  # working-independence PGEE estimates, while V is V1 + V2 - V3.
  beta <- colMeans(rbind(coef(mod1), coef(mod2), coef(mod3)), na.rm = TRUE)

  mod_use <- mod1
  mod_use$coefficients <- beta
  mod_use$covb <- V

  list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod_use = mod_use,
       beta = beta, V = V)
}

nearest_psd <- function(V, eps = 1e-8) {
  V <- (V + t(V)) / 2
  ee <- eigen(V, symmetric = TRUE)
  vals <- pmax(ee$values, eps)
  out <- ee$vectors %*% diag(vals, length(vals)) %*% t(ee$vectors)
  dimnames(out) <- dimnames(V)
  out
}

result_table <- function(beta, V, treatment_terms) {
  est <- beta[treatment_terms]
  se <- sqrt(diag(V))[treatment_terms]
  z <- est / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  PI <- expit(est)
  PI_low <- expit(est - 1.96 * se)
  PI_high <- expit(est + 1.96 * se)
  data.frame(
    contrast = treatment_terms,
    estimate_logit_PI = est,
    se = se,
    z = z,
    p_unadjusted = p,
    p_holm = p.adjust(p, method = "holm"),
    PI = PI,
    PI_low = PI_low,
    PI_high = PI_high,
    win_odds = PI / (1 - PI),
    row.names = NULL
  )
}

full_model_result_table <- function(beta, V, treatment_terms) {
  # Full model output: both within-group time effects and visit-specific
  # treatment effects. Estimates are on the logit probabilistic-index scale.
  nm <- names(beta)
  est <- as.numeric(beta)
  se <- sqrt(diag(V))
  z <- est / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)

  term_type <- ifelse(nm %in% treatment_terms,
                      "Treatment effect",
                      "Time effect")
  term_label <- nm
  term_label[term_label == "trend_treat"] <- "Time trend, treatment group"
  term_label[term_label == "trend_ctrl"]  <- "Time trend, control group"
  term_label[nm %in% treatment_terms] <- paste("Treatment effect, visit",
                                                sub("trt_visit", "", nm[nm %in% treatment_terms]))

  PI <- expit(est)
  PI_low <- expit(est - 1.96 * se)
  PI_high <- expit(est + 1.96 * se)

  out <- data.frame(
    term = nm,
    label = term_label,
    type = term_type,
    estimate_logit_PI = est,
    se = se,
    z = z,
    p_unadjusted = p,
    PI = PI,
    PI_low = PI_low,
    PI_high = PI_high,
    win_odds = PI / (1 - PI),
    row.names = NULL
  )
  out
}

make_pairwise_interaction_L <- function(beta_names, treatment_terms) {
  # Diagnostic global Wald test: equality of visit-specific treatment effects by
  # comparing the first visit with all remaining visits.
  if (length(treatment_terms) < 2) stop("Need at least two treatment terms.")
  L <- matrix(0, nrow = length(treatment_terms) - 1, ncol = length(beta_names))
  colnames(L) <- beta_names
  rownames(L) <- paste0(treatment_terms[1], " - ", treatment_terms[-1])
  for (i in seq_along(treatment_terms[-1])) {
    L[i, treatment_terms[1]] <- 1
    L[i, treatment_terms[i + 1]] <- -1
  }
  L
}

make_deviation_from_mean_L <- function(beta_names, treatment_terms) {
  # Holm procedure used in the simulations and manuscript:
  # H0t: beta_At - mean_t(beta_At) = 0 for each visit t.
  K <- length(treatment_terms)
  if (K < 2) stop("Need at least two treatment terms.")
  L <- matrix(0, nrow = K, ncol = length(beta_names))
  colnames(L) <- beta_names
  rownames(L) <- paste0(treatment_terms, " - mean(treatment effects)")
  for (i in seq_len(K)) {
    L[i, treatment_terms] <- -1 / K
    L[i, treatment_terms[i]] <- 1 - 1 / K
  }
  L
}


write_latex_table <- function(df, path, digits = 3) {
  fmt <- function(x) if (is.numeric(x)) sprintf(paste0("%.", digits, "f"), x) else as.character(x)
  out <- df
  for (nm in names(out)) if (is.numeric(out[[nm]])) out[[nm]] <- fmt(out[[nm]])
  lines <- c("\\begin{tabular}{lrrrrrr}",
             "\\toprule",
             paste(names(out), collapse = " & "), "\\\\",
             "\\midrule")
  rows <- apply(out, 1, paste, collapse = " & ")
  lines <- c(lines, paste0(rows, " \\\\"), "\\bottomrule", "\\end{tabular}")
  writeLines(lines, path)
}

## -------------------------------------------------------------------------
## 3. Analysis settings
## -------------------------------------------------------------------------
set.seed(1988)
n_per_group <- 15
higher_is_better <- TRUE
out_dir <- "data_example_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

## -------------------------------------------------------------------------
## 4. Load data, construct pseudo-observations, fit three-GEE procedure
## -------------------------------------------------------------------------

arthritis_small <- make_small_arthritis(n_per_group = n_per_group, seed = 1988)
pairdat <- make_compare_and_dat_gee(arthritis_small,
                                    higher_is_better = higher_is_better,
                                    include_within_subject_pairs = TRUE)

cat("Small arthritis subset:\n")
print(arthritis_small %>% distinct(id, group) %>% count(group))
cat("\nPseudo-observations by type:\n")
print(table(pairdat$compare$pair_type))
cat("\nModel covariates:\n")
print(pairdat$x_names)

fit <- fit_three_geessbin_FW(pairdat$dat_GEE)

# Ensure that coefficient and covariance names match the design-matrix names.
# Some versions of geessbin rename columns generated by the formula interface.
# If this happens, downstream extraction of trt_visit* terms returns an empty
# data frame and the treatment-effect plots are blank.
if (length(fit$beta) == length(pairdat$x_names)) {
  names(fit$beta) <- pairdat$x_names
  dimnames(fit$V) <- list(pairdat$x_names, pairdat$x_names)
  fit$mod_use$coefficients <- fit$beta
  fit$mod_use$covb <- fit$V
}

# Recreate the treatment-term vector from the fitted coefficient names.
pairdat$treatment_terms <- grep("^trt_visit", names(fit$beta), value = TRUE)
if (length(pairdat$treatment_terms) == 0) {
  stop("No treatment-effect terms detected in the fitted model. Coefficient names are: ",
       paste(names(fit$beta), collapse = ", "))
}


# If the combined covariance has numerical negative eigenvalues, use a PSD
# projection only for glht; raw V is still saved.
V_raw <- fit$V
V_for_inference <- V_raw
asymmetry <- max(abs(V_raw - t(V_raw)), na.rm = TRUE)
cat("
Maximum asymmetry in V1 + V2 - V3: ", signif(asymmetry, 3), "
", sep = "")
# The estimator is V1 + V2 - V3. A nearest-PSD projection is only used as a
# numerical fallback if the raw matrix cannot be used for inference.
V_eig_check <- (V_raw + t(V_raw)) / 2
if (min(eigen(V_eig_check, symmetric = TRUE, only.values = TRUE)$values) < -1e-8 || any(diag(V_raw) <= 0)) {
  warning("Combined V has negative eigenvalues or non-positive variances; using nearest PSD matrix for numerical inference. Raw V is still saved.")
  V_for_inference <- nearest_psd(V_raw)
}
fit$mod_use$covb <- V_for_inference

# Full model output and descriptive visit-specific treatment effects.
full_param_results <- full_model_result_table(fit$beta, V_for_inference, pairdat$treatment_terms)
param_results <- result_table(fit$beta, V_for_inference, pairdat$treatment_terms)

# Primary Holm test for treatment-effect constancy:
# H0t: beta_At - mean_t(beta_At) = 0.
L_const <- make_deviation_from_mean_L(names(fit$beta), pairdat$treatment_terms)
glht_holm_const <- summary(multcomp::glht(fit$mod_use, linfct = L_const),
                           test = multcomp::adjusted("holm"))


holm_results <- data.frame(
  contrast = rownames(L_const),
  estimate = as.numeric(glht_holm_const$test$coefficients),
  se = as.numeric(glht_holm_const$test$sigma),
  z = as.numeric(glht_holm_const$test$tstat),
  p_holm = as.numeric(glht_holm_const$test$pvalues),
  row.names = NULL
)
holm_results$p_unadjusted <- 2 * stats::pnorm(abs(holm_results$z), lower.tail = FALSE)
holm_results <- holm_results[, c("contrast", "estimate", "se", "z", "p_unadjusted", "p_holm")]

# Diagnostic global Wald interaction test: equality of treatment effects across visits.
L_inter <- make_pairwise_interaction_L(names(fit$beta), pairdat$treatment_terms)
out_global <- aod::wald.test(V_for_inference, fit$beta, L = L_inter,
                             df = length(unique(arthritis_small$id)) - qr(L_inter)$rank)

global_stat <- as.numeric(out_global$result$chi2["chi2"])
global_p <- as.numeric(out_global$result$chi2["P"])
global_df <- nrow(L_inter)


# Black-and-white trajectory plot.
p <- ggplot(arthritis_small, aes(x = time, y = y, group = id, linetype = factor(group), shape = factor(group))) +
  geom_line(alpha = 0.55, colour = "black") +
  geom_point(alpha = 0.85, colour = "black", size = 2) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Treatment group") +
  scale_shape_manual(values = c(16, 17), name = "Treatment group") +
  labs(x = "Visit", y = "Ordinal response") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "arthritis_trajectories_bw.png"), p, width = 7, height = 4.5, dpi = 300)

## -------------------------------------------------------------------------
## Secondary visit-specific treatment tests: H0 beta_A(t) = 0
## Ordinary pointwise 95% CIs + Holm-adjusted p-values
##
## These p-values answer a different question from the constancy contrasts.
## They test whether the treatment effect differs from zero at each visit.
## The ordinary 95% CIs are shown for effect-size interpretation; the Holm
## p-values should be used for simultaneous statements across visits.
## -------------------------------------------------------------------------

visit_specific_tests <- param_results %>%
  dplyr::mutate(
    visit = dplyr::row_number(),
    low95 = estimate_logit_PI - qnorm(0.975) * se,
    high95 = estimate_logit_PI + qnorm(0.975) * se,
    p_raw = 2 * pnorm(-abs(z)),
    p_holm = p.adjust(p_raw, method = "holm"),
    PI_low95 = plogis(low95),
    PI_high95 = plogis(high95),
    p_holm_label = dplyr::case_when(
      is.na(p_holm) ~ "Holm p = NA",
      p_holm < 0.001 ~ "Holm p < 0.001",
      TRUE ~ paste0("Holm p = ", formatC(p_holm, format = "f", digits = 3))
    )
  )


## Figure 5: visit-specific effects with ordinary 95% CIs and Holm p-values.
## No file is read here; the figure is generated directly from param_results.

if (nrow(visit_specific_tests) == 0 || all(is.na(visit_specific_tests$estimate_logit_PI))) {
  stop("Figure 5 cannot be created: visit_specific_tests is empty or contains only NA estimates. Check treatment-term names and V matrix.")
}
print(visit_specific_tests)

y_range <- range(c(visit_specific_tests$low95, visit_specific_tests$high95), na.rm = TRUE)
y_pad <- 0.14 * diff(y_range)
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.10

## Horizontal reference line for the average visit-specific treatment effect.
mean_effect <- mean(expit(visit_specific_tests$estimate_logit_PI), na.rm = TRUE)

## Put all Holm p-value labels at a common height, above the CIs.
## The x-axis is deliberately widened so the label at Visit 3 is not clipped.
visit_specific_tests$label_y <- max(expit(visit_specific_tests$high95+ y_pad), na.rm = TRUE) 

p_eff <- ggplot(visit_specific_tests, aes(x = visit, y = expit(estimate_logit_PI))) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "black") +
  geom_hline(yintercept = mean_effect, linetype = "dashed", colour = "black") +
  annotate(
    "text",
    x = max(visit_specific_tests$visit) + 0.30,
    y = mean_effect,
    label = "Mean effect",
    hjust = 0,
    vjust = -0.35,
    size = 3.0,
    colour = "black"
  ) +
  geom_errorbar(aes(ymin = expit(low95), ymax = expit(high95)), width = 0.08, colour = "black") +
  geom_point(shape = 16, size = 2.5, colour = "black") +
  geom_text(
    aes(y = label_y, label = p_holm_label),
    size = 3.1,
    colour = "black",
    check_overlap = FALSE
  ) +
  scale_x_continuous(
    breaks = visit_specific_tests$visit,
    limits = c(min(visit_specific_tests$visit) - 0.35,
               max(visit_specific_tests$visit) + 0.75)
  ) +
  scale_y_continuous(
    breaks = c(0.5,0.7,1),
    limits = c(0.45,1)
  ) +
  labs(
    x = "Visit",
    y = "Treatment effect on PI scale"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.margin = margin(t = 8, r = 35, b = 8, l = 8)
  )


cat("\nFull parameter estimates for proposed longitudinal PIM:\n")
print(full_param_results)
cat("\nHolm-adjusted tests for treatment-effect constancy:\n")
print(holm_results)

ggsave(file.path(out_dir, "arthritis_figure5_ordinary95CI_HolmP_bw.png"),
       p_eff, width = 6.5, height = 4.5, dpi = 300)

cat("\nDone. Output written to: ", normalizePath(out_dir), "\n", sep = "")


