## Simulations related to LSPIM, as introduced in the paper "Longitudinal Probabilistic Index Models:
## Small-Sample Inference and Tests for Interaction"

## Two data generating models are considered, i.e. random intercept and proportional odds model, 
## where the latter is based on the model used by Long et al. (2026)

## Load packages and functions
library(dplyr)
library(aod)
library(geessbin)
source("functions.R")
source("lwo.R")          # Fuctions from Long et al. (2026)
library(pim)
library(tidyverse)
library(aod)
library(multcomp)
library(mvtnorm)
library(nparLD)
library(simstudy)
library(reshape2)

## We need to adapt the geessbin function to allow for ties

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
  if (!is.numeric(y) | !setequal(unique(y), c(0,0.5,1))) {
    stop("outcome vector must be numeric and take values in {0, 1}")
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

nr_of_iters = 10000

## RI data generating

# Contrast matrices considered for RI model

L <- rbind(
  c(1, -1,rep(0,4))
)


L_new <- rbind(
  c(0,0,1,-1,0,0),
  c(0,0,1,0,-1,0),
  c(0,0,1,0,0,-1)
)

L_new2 = L_test = matrix(c(c(0,0,c(3/4,rep(-1/4,3))),
                           c(0,0,c(-1/4,3/4,rep(-1/4,2))),
                           c(0,0,c(rep(-1/4,2),3/4,-1/4)),
                           c(0,0,c(rep(-1/4,3),3/4))),nrow=4,byrow = TRUE)


for(N in c(16,18,20,24,26,30,40,50,60)){
  nr_of_subjects=N
  use_times = nr_of_times = 4
  id = rep(1:nr_of_subjects,each=nr_of_times)
  time = rep(1:nr_of_times,nr_of_subjects)
  group = rep(rep(c(0,1),each=nr_of_subjects/2),each=nr_of_times)
  
  beta0=0
  beta1=sqrt(2)*0.1
  beta2=sqrt(2)*0.1
  beta3 =sqrt(2)*0
  

  design  = data.frame(id,time,group)

  
  store_parms_mod2 = matrix(nrow=nr_of_iters,ncol = 6)
  store_parms_mod2_adj = matrix(nrow=nr_of_iters,ncol = 6)
  store_se_mod2 = matrix(nrow=nr_of_iters,ncol = 6)
  store_se_mod2_adj = matrix(nrow=nr_of_iters,ncol = 6)
  
  store_parms_gee = matrix(nrow=nr_of_iters,ncol = 6)
  store_se_gee = matrix(nrow=nr_of_iters,ncol = 6)
  
  store_interaction_tests = matrix(nrow=nr_of_iters,ncol = 18)
  score_inter_group = c()
  score_inter_time = c()
  score_inter_group_global = c()
  score_inter_time_global = c()
  
  for(l in 1:nr_of_iters){
    if(l%%100==0){cat("Iteration: ",l,"\n")}
    set.seed(l)
    re = rmvnorm(nr_of_subjects,c(0,0),matrix(c(0.8^2,0,0,1.5^2),nrow=2))
    design$Y = beta0+
      (beta1)*design$time+
      beta2*(design$group==1)+
      beta3*design$time*(design$group==1)+
      rnorm(nrow(design),0,1)+
      rep(re[,2],each=nr_of_times)
    
    design$Y = round(design$Y)
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
    
    # Model 2
    assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")
    mod2_orig = pim(Y~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) +I((R(group)-L(group))*(R(time)==4)*(L(time)==4)),data=design,compare=compare,link="probit")
    store_parms_mod2[l,] = coef(mod2_orig)
    store_se_mod2[l,] = diag(vcov(mod2_orig))
    
    assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")
    
    mod2_adj = pim(Y~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) +I((R(group)-L(group))*(R(time)==4)*(L(time)==4)),data=design,compare=compare,link="probit")
    store_parms_mod2_adj[l,] = coef(mod2_adj)
    store_se_mod2_adj[l,] = diag(vcov(mod2_adj))
    
    C1 = individuals_var1
    C2 = individuals_var2
    
    
    x = model.matrix(mod2_adj)
    y = response(mod2_adj)
    
    
    dat_GEE = data.frame(y,x)
    
    names(dat_GEE) = c("y",c(paste0("x",1:6)))
    dat_GEE$C1=C1
    dat_GEE$C2=C2
    
    
    dat_GEE$C3 <- paste(dat_GEE$C1, dat_GEE$C2,sep="_")
    
    
    
    dat_GEE=dat_GEE[order(dat_GEE$C1),]
    mod1 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C1, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C2),]
    mod2 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C2, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C3),]
    mod3 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C3, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    
    store_se_gee[l,] = diag(mod1$covb+mod2$covb-mod3$covb)
    store_parms_gee[l,] = apply(rbind(coef(mod1),coef(mod2),coef(mod3)),2,mean)
    
    
    ### Interaction tests based on adjusted PIM
    
    V=vcov(mod2_adj)
    b = coef(mod2_adj)
    
    ### Standard degrees of freedom
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nrow(L_new), verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
    
    while(is.null(out)){
      V = V+rnorm(length(V),0.0001,0.00001)
      out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                               df = nrow(L_new), verbose = FALSE), error = function(e) {
                                 return(NULL)
                               })
    }
    
    store_interaction_tests[l,1] = out$result$chi2[3]
    store_interaction_tests[l,2] = out$result$F[4]
    
    V=vcov(mod2_adj)
    b = coef(mod2_adj)
    
    ### Between-Within degrees of freedom
    out = wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                    df = nr_of_subjects-qr(L_new)$rank, verbose = FALSE)
    
    # while(is.null(out)){
    #   V = V+rnorm(length(V),0.0001,0.00001)
    #   out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
    #                            df = nr_of_subjects-qr(L_new)$rank, verbose = FALSE), error = function(e) {
    #                              return(NULL)
    #                            })
    # }
    
    store_interaction_tests[l,3] = out$result$F[4]
    
    ### Containment degrees of freedom
    
    V=vcov(mod2_adj)
    b = coef(mod2_adj)
    
    
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
    
    while(is.null(out)){
      V = V+rnorm(length(V),0.0001,0.00001)
      out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                               df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE), error = function(e) {
                                 return(NULL)
                               })
    }
    
    store_interaction_tests[l,4] = out$result$F[4]
    
    ## Same tests, but for the time difference
    
    V=vcov(mod2_adj)
    b = coef(mod2_adj)
    
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nrow(L), verbose = FALSE)
    
    store_interaction_tests[l,5] = out$result$chi2[3]
    store_interaction_tests[l,6] = out$result$F[4]
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nr_of_subjects-qr(L)$rank, verbose = FALSE)
    
    store_interaction_tests[l,7] = out$result$F[4]
    
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE)
    
    store_interaction_tests[l,8] = out$result$F[4]
    
    
    ### Interaction tests based on GEE approach
    
    V=mod1$covb+mod2$covb-mod3$covb
    b = apply(rbind(coef(mod1),coef(mod2),coef(mod3)),2,mean)
    ### Standard degrees of freedom
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nrow(L_new), verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
    
    while(is.null(out)){
      V = V+rnorm(length(V),0.0001,0.00001)
      out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                               df = nrow(L_new), verbose = FALSE), error = function(e) {
                                 return(NULL)
                               })
    }
    
    store_interaction_tests[l,9] = out$result$chi2[3]
    store_interaction_tests[l,10] = out$result$F[4]
    
    
    ### Between-Within degrees of freedom
    
    V=mod1$covb+mod2$covb-mod3$covb
    b = apply(rbind(coef(mod1),coef(mod2),coef(mod3)),2,mean)
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nr_of_subjects-qr(L_new)$rank, verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
    
    store_interaction_tests[l,11] = out$result$F[4]
    
    ### Containment degrees of freedom
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
    
    while(is.null(out)){
      V = V+rnorm(length(V),0.0001,0.00001)
      out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                               df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE), error = function(e) {
                                 return(NULL)
                               })
    }
    
    store_interaction_tests[l,12] = out$result$F[4]
    
    
    
    ## Same tests, but for the time difference
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nrow(L), verbose = FALSE)
    
    store_interaction_tests[l,13] = out$result$chi2[3]
    store_interaction_tests[l,14] = out$result$F[4]
    
    V=mod1$covb+mod2$covb-mod3$covb
    b = apply(rbind(coef(mod1),coef(mod2),coef(mod3)),2,mean)
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nr_of_subjects-qr(L)$rank, verbose = FALSE)
    
    store_interaction_tests[l,15] = out$result$F[4]
    
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE)
    
    store_interaction_tests[l,16] = out$result$F[4]
    
    
    mod_narlD <- nparLD(Y ~ time * group, data = design,
                        subject = "id", description = FALSE)
    
    
    store_interaction_tests[l,17] = mod_narlD$Wald.test[3,3]
    store_interaction_tests[l,18] =  mod_narlD$ANOVA.test[3,3]
    
    mod_use = mod1
    mod_use$covb=V
    
    out = summary(glht(mod_use, linfct = L_test), test = adjusted("holm"))
    score_inter_group = c(score_inter_group,!prod(out$test$pvalues>0.05))
    mod_use$df.residual = nr_of_subjects-qr(L_test)$rank-1
    out = summary(glht(mod_use, linfct = L_test), test = Ftest())
    score_inter_group_global = c(score_inter_group_global,out$test$pvalue)
    mod_use$df.residual = nr_of_subjects-qr(L)$rank-1
    out = summary(glht(mod_use, linfct = L), test = adjusted("holm"))
    score_inter_time = c(score_inter_time,!prod(out$test$pvalues>0.05))
    out = summary(glht(mod_use, linfct = L), test = Ftest())
    score_inter_time_global = c(score_inter_time_global,out$test$pvalue)
  }  
  
  keep = which((abs(store_parms_mod2[,3])<5)&(abs(store_parms_mod2[,4])<5)&(abs(store_parms_mod2[,5])<5)&(abs(store_parms_mod2[,6])<5)&
                 (abs(store_parms_mod2[,7])<5)&(abs(store_parms_mod2[,8])<5))
  type1_interact = apply(store_interaction_tests[keep,],2,function(x)mean(x<0.05,na.rm=TRUE))
  
  save.image(sprintf("Simulation_output/small_sample_type1_pim_gee_%s_correct_PGEE_score_unrounded_testwithmeanhypothesis_4timepoints.Rdata",
      N))
}


## PO data generating

## Contrast matrices considered for PO model

L <- rbind(
  c(1, -1,rep(0,5))
)

L_new <- rbind(
  c(0,0,0,1,-1,0,0),
  c(0,0,0,1,0,-1,0),
  c(0,0,0,1,0,0,-1))

L_test <- rbind(
  c(0, 0,0, 3/4, -1/4,-1/4, -1/4),
  c(0, 0,0, -1/4, 3/4,-1/4, -1/4),
  c(0, 0,0, -1/4, -1/4, 3/4, -1/4),
  c(0, 0,0, -1/4, -1/4, -1/4, 3/4)
)

visits <- c(0,4,8,12,16) # baseline is visit = 0 
nr_of_times = length(visits)
# base ordinal category probabilities if all covariates are zero
baseprobs <- rev(c(0.06,0.11,0.12,0.50,0.21)) # reversed so larger outcomes are better
# covariate effects 
covs_effects_baseline <- c("trt"=0,"age"=-0.005,"pre_diarrhea"= 0.23)
corstr <- "ar1"
rho <- 0.6
corMatrix <- generate_corMatrix(n_visits = length(visits),rho=rho,corstr = corstr)
# simulation scenario, one of null, SID, pos and con
scenario <- "pos"
time_effect <- 0.15

trt_effect <-  0.2

time_trt_interaction <-  0

# working correlation used for modeling in LWO method
corstr_working <- "independence"

for(N in c(16,18,20,24,26,30,40,50,60)){
  nr_of_subjects = N
  store_parms_mod2 = matrix(nrow=nr_of_iters,ncol = 7)
  store_parms_mod2_adj = matrix(nrow=nr_of_iters,ncol = 7)
  store_se_mod2 = matrix(nrow=nr_of_iters,ncol = 7)
  store_se_mod2_adj = matrix(nrow=nr_of_iters,ncol = 7)
  
  store_parms_gee = matrix(nrow=nr_of_iters,ncol = 7)
  store_se_gee = matrix(nrow=nr_of_iters,ncol = 7)
  
  store_interaction_tests = matrix(nrow=nr_of_iters,ncol = 18)
  score_inter_group = c()
  score_inter_time = c()
  score_inter_group_global = c()
  score_inter_time_global = c()
  
  for(sims in 1:nr_of_iters){
    if(sims%%100==0){cat("Iteration: ",sims,"\n")}
    l = sims
    set.seed(sims)
    dat <- gen_data(N,baseprobs,covs_effects_baseline,
                    time_effect,trt_effect,time_trt_interaction,
                    visits,corMatrix)
    dat_wide <- dat$wide_format
    dat_long <- dat$long_format
    
    # create indicator for week 4 and week 8
    dat_long <- dat_long |>
      mutate(w4 = 1*(time==4),w8 = 1*(time==8),w12 = 1*(time==12),w16 = 1*(time==16))
    
    
    # Define compare for PIM  (both within and between subjects are compared)
    
    id.fac1 <- which((dat_long[,"trt"] == "0")&(dat_long[,"visit_cat"] == "visit1"))
    id.nonfac1 <- which((dat_long[,"trt"] == "1")&(dat_long[,"visit_cat"] == "visit1"))
    compare1 <- expand.grid(id.fac1,id.nonfac1)
    
    id.fac2 <- which((dat_long[,"trt"] == "0")&(dat_long[,"visit_cat"] == "visit2"))
    id.nonfac2 <- which((dat_long[,"trt"] == "1")&(dat_long[,"visit_cat"] == "visit2"))
    compare2 <- expand.grid(id.fac2,id.nonfac2)
    
    id.fac3 <- which((dat_long[,"trt"] == "0")&(dat_long[,"visit_cat"] == "visit3"))
    id.nonfac3 <- which((dat_long[,"trt"] == "1")&(dat_long[,"visit_cat"] == "visit3"))
    compare3 <- expand.grid(id.fac3,id.nonfac3)
    
    id.fac4 <- which((dat_long[,"trt"] == "0")&(dat_long[,"visit_cat"] == "visit4"))
    id.nonfac4 <- which((dat_long[,"trt"] == "1")&(dat_long[,"visit_cat"] == "visit4"))
    compare4 <- expand.grid(id.fac4,id.nonfac4)
    
    id.fac5 <- which((dat_long[,"trt"] == "0")&(dat_long[,"visit_cat"] == "visit5"))
    id.nonfac5 <- which((dat_long[,"trt"] == "1")&(dat_long[,"visit_cat"] == "visit5"))
    compare5 <- expand.grid(id.fac5,id.nonfac5)
    
    compare_between = rbind(compare1,compare2,compare3,compare4,compare5)
    compare_between = compare_between[order(compare_between[,"Var2"]),]
    
    start_timepoints <- c()
    later_timepoints <- c()
    k=1
    for (i in 1:(nrow(dat_long) - 1)) {
      if(i%%5 != 0){
        start <- i           # Current time point
        later <- (i+1) : (5*k)  # Later time points
        # Store the results in vectors
        start_timepoints <- c(start_timepoints, rep(start, length(later)))
        later_timepoints <- c(later_timepoints, later)}
      if(i%%5 == 0){
        k=k+1}
    }
    
    within1 = start_timepoints
    within2 = later_timepoints
    
    compare_within = cbind(within1,within2)
    colnames(compare_within) = c("Var1","Var2")
    compare = rbind(compare_between,compare_within)
    
    individuals_var1=dat_long[compare$Var1,"id"]%>%unlist()
    individuals_var2=dat_long[compare$Var2,"id"]%>%unlist()
    
    # Model 2
    assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")
    mod2_orig = pim(GBS_DS~I((R(time) - L(time))*R(trt)*L(trt))+I((R(time) - L(time))*(1-R(trt))*(1-L(trt)))+I((R(trt)-L(trt))*(R(visit_cat)=="visit1")*(L(visit_cat)=="visit1"))+I((R(trt)-L(trt))*(R(visit_cat)=="visit2")*(L(visit_cat)=="visit2"))+I((R(trt)-L(trt))*(R(visit_cat)=="visit3")*(L(visit_cat)=="visit3"))+
                      +I((R(trt)-L(trt))*(R(visit_cat)=="visit4")*(L(visit_cat)=="visit4"))+I((R(trt)-L(trt))*(R(visit_cat)=="visit5")*(L(visit_cat)=="visit5")),data=dat_long,compare=compare,link="probit")
    store_parms_mod2[l,] = coef(mod2_orig)
    store_se_mod2[l,] = diag(vcov(mod2_orig))
    
    assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")
    
    mod2_adj = pim(GBS_DS~I((R(time) - L(time))*R(trt)*L(trt))+I((R(time) - L(time))*(1-R(trt))*(1-L(trt)))+I((R(trt)-L(trt))*(R(visit_cat)=="visit1")*(L(visit_cat)=="visit1"))+I((R(trt)-L(trt))*(R(visit_cat)=="visit2")*(L(visit_cat)=="visit2"))+I((R(trt)-L(trt))*(R(visit_cat)=="visit3")*(L(visit_cat)=="visit3"))+
                     +I((R(trt)-L(trt))*(R(visit_cat)=="visit4")*(L(visit_cat)=="visit4"))+I((R(trt)-L(trt))*(R(visit_cat)=="visit5")*(L(visit_cat)=="visit5")),data=dat_long,compare=compare,link="probit")
    
    store_se_mod2_adj[l,] = diag(vcov(mod2_adj))
    
    C1 = individuals_var1
    C2 = individuals_var2
    
    
    x = model.matrix(mod2_adj)
    y = response(mod2_adj)
    
    
    dat_GEE = data.frame(y,x)
    
    names(dat_GEE) = c("y",c(paste0("x",1:7)))
    dat_GEE$C1=C1
    dat_GEE$C2=C2
    
    
    dat_GEE$C3 <- paste(dat_GEE$C1, dat_GEE$C2,sep="_")
    
    
    
    dat_GEE=dat_GEE[order(dat_GEE$C1),]
    mod1 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C1, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C2),]
    mod2 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C2, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C3),]
    mod3 = geessbin_adjusted(y~.-1-C1-C2-C3,data=dat_GEE,id=C3, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    
    store_se_gee[l,] = diag(mod1$covb+mod2$covb-mod3$covb)
    store_parms_gee[l,] = apply(rbind(coef(mod1),coef(mod2),coef(mod3)),2,mean)
    
    
    ### Interaction tests based on adjusted PIM
    
    V=vcov(mod2_adj)
    b = coef(mod2_adj)
    
    ### Standard degrees of freedom
    
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nrow(L_new), verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
  
    
    store_interaction_tests[l,1] = out$result$chi2[3]
    store_interaction_tests[l,2] = out$result$F[4]
    
    V=vcov(mod2_adj)
    b = coef(mod2_adj)
    
    ### Between-Within degrees of freedom
    out = wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                    df = nr_of_subjects-qr(L_new)$rank, verbose = FALSE)
    
    
    store_interaction_tests[l,3] = out$result$F[4]
    
    ### Containment degrees of freedom
    
    V=vcov(mod2_adj)
    b = coef(mod2_adj)
    
    
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
    
    
    store_interaction_tests[l,4] = out$result$F[4]
    
    
    
    ## Same tests, but for the time difference
    
    V=vcov(mod2_adj)
    b = coef(mod2_adj)
    
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nrow(L), verbose = FALSE)
    
    store_interaction_tests[l,5] = out$result$chi2[3]
    store_interaction_tests[l,6] = out$result$F[4]
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nr_of_subjects-qr(L)$rank, verbose = FALSE)
    
    store_interaction_tests[l,7] = out$result$F[4]
    
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE)
    
    store_interaction_tests[l,8] = out$result$F[4]
    
    
    
    
    ### Interaction tests based on GEE approach
    
    V=mod1$covb+mod2$covb-mod3$covb
    b = apply(rbind(coef(mod1),coef(mod2),coef(mod3)),2,mean)
    ### Standard degrees of freedom
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nrow(L_new), verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
    
    while(is.null(out)){
      V = V+rnorm(length(V),0.0001,0.00001)
      out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                               df = nrow(L_new), verbose = FALSE), error = function(e) {
                                 return(NULL)
                               })
    }
    
    store_interaction_tests[l,9] = out$result$chi2[3]
    store_interaction_tests[l,10] = out$result$F[4]
    
    
    ### Between-Within degrees of freedom
    
    V=mod1$covb+mod2$covb-mod3$covb
    b = apply(rbind(coef(mod1),coef(mod2),coef(mod3)),2,mean)
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nr_of_subjects-qr(L_new)$rank, verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
    
   
    store_interaction_tests[l,11] = out$result$F[4]
    
    ### Containment degrees of freedom
    
    out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                             df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE), error = function(e) {
                               return(NULL)
                             })
    
    while(is.null(out)){
      V = V+rnorm(length(V),0.0001,0.00001)
      out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
                               df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE), error = function(e) {
                                 return(NULL)
                               })
    }
    
    store_interaction_tests[l,12] = out$result$F[4]
    
    
    ## Same tests, but for the time difference
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nrow(L), verbose = FALSE)
    
    store_interaction_tests[l,13] = out$result$chi2[3]
    store_interaction_tests[l,14] = out$result$F[4]
    
    V=mod1$covb+mod2$covb-mod3$covb
    b = apply(rbind(coef(mod1),coef(mod2),coef(mod3)),2,mean)
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nr_of_subjects-qr(L)$rank, verbose = FALSE)
    
    store_interaction_tests[l,15] = out$result$F[4]
    
    
    out = wald.test(V, b, Terms = NULL, L = L, H0 = NULL,
                    df = nr_of_subjects*nr_of_times - nr_of_subjects, verbose = FALSE)
    
    store_interaction_tests[l,16] = out$result$F[4]
    
    
    mod_nparLD <- nparLD(GBS_DS ~ visit_cat * trt, data = dat_long,
                         subject = "id", description = FALSE)
    
    
    store_interaction_tests[l,17] = mod_nparLD$Wald.test[3,3]
    store_interaction_tests[l,18] =  mod_nparLD$ANOVA.test[3,3]
    
    mod_use = mod1
    mod_use$covb=V
    
    out = summary(glht(mod_use, linfct = L_test), test = adjusted("holm"))
    score_inter_group = c(score_inter_group,!prod(out$test$pvalues>0.05))
    mod_use$df.residual = nr_of_subjects-qr(L_test)$rank-1
    out = summary(glht(mod_use, linfct = L_test), test = Ftest())
    score_inter_group_global = c(score_inter_group_global,out$test$pvalue)
    mod_use$df.residual = nr_of_subjects-qr(L)$rank-1
    out = summary(glht(mod_use, linfct = L), test = adjusted("holm"))
    score_inter_time = c(score_inter_time,!prod(out$test$pvalues>0.05))
    out = summary(glht(mod_use, linfct = L), test = Ftest())
    score_inter_time_global = c(score_inter_time_global,out$test$pvalue)
  }  
  
  save.image(paste0("Simulation_output/small_sample_type1_pim_gee_",N,"_correct_PGEE_POmodel_testwithmeanhypothesis.Rdata"))

  }

