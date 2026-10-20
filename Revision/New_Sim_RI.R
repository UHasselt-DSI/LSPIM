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

# Gaussian RI model
use_times = nr_of_times = 4
beta0=0
beta1=sqrt(2)*0.1
beta2=sqrt(2)*0.1
beta3 =sqrt(2)*0

gaps <- 1:3
weights <- 4 - gaps

beta_T_true <- uniroot(
  function(b) {
    sum(weights * gaps *
          (pnorm(beta1 * gaps / sqrt(2)) - plogis(b * gaps)))
  },
  c(-10, 10)
)$root

beta_A_true <- qlogis(
  pnorm(beta2 / sqrt(2 * (1 + 1.5^2)))
)

truth <- c(
  rep(beta_T_true, 2),
  rep(beta_A_true, 4)
)

truth


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
    re = rmvnorm(nr_of_subjects,c(0,0),matrix(c(0.8^2,0,0,1.5^2),nrow=2))
    design$Y = beta0+
      (beta1)*design$time+
      beta2*(design$group==1)+
      beta3*design$time*(design$group==1)+
      rnorm(nrow(design),0,1)+
      rep(re[,2],each=nr_of_times)
    
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
    mod2_orig  <- capture_fit_warnings(
      pim(Y~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) +I((R(group)-L(group))*(R(time)==4)*(L(time)==4)),data=design,compare=compare,link="logit",start = truth),
      model = "Original PIM",
      run = l,
      N = nr_of_subjects
    )

    V_PIM = vcov(mod2_orig)
    b_PIM = coef(mod2_orig)
    
    # Fit sandwich adjusted PIM
    assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")
    mod2_adj  <- capture_fit_warnings(
      pim(Y~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) +I((R(group)-L(group))*(R(time)==4)*(L(time)==4)),data=design,compare=compare,link="logit",start=truth),
      model = "Adjusted PIM",
      run = l,
      N = nr_of_subjects
    )

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
    
    head(dat_GEE)
    
    # Fit models

    dat_GEE=dat_GEE[order(dat_GEE$C1),]
    mod1 = geessbin(y~.-1-C1-C2-C3,data=dat_GEE,id=C1, corstr = "independence",beta.method="GEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C2),]
    mod2 = geessbin(y~.-1-C1-C2-C3,data=dat_GEE,id=C2, corstr = "independence",beta.method="GEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C3),]
    mod3 = geessbin(y~.-1-C1-C2-C3,data=dat_GEE,id=C3, corstr = "independence",beta.method="GEE",SE.method = "FW")

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
    coef(mod1)
    dat_GEE=dat_GEE[order(dat_GEE$C1),]
    mod1 = geessbin(y~.-1-C1-C2-C3,data=dat_GEE,id=C1, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C2),]
    mod2 = geessbin(y~.-1-C1-C2-C3,data=dat_GEE,id=C2, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    dat_GEE=dat_GEE[order(dat_GEE$C3),]
    mod3 = geessbin(y~.-1-C1-C2-C3,data=dat_GEE,id=C3, corstr = "independence",beta.method="PGEE",SE.method = "FW")
    
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
    
    
    
    
    
    models$PGEE
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







old <- readRDS("/Users/lucp2490/Library/CloudStorage/GoogleDrive-stijn.jaspers@uhasselt.be/My Drive/FWO Postdoc copy/PIM/Longitudinal PIM/old_run.rds")

# Outcomes
identical(as.numeric(old$y), as.numeric(y))

mod1

for(l in 1:100){
pim_parms = parm(b_PIM_adj ,    V_PIM_adj , N - 1)
out_test = glht(pim_parms,linfct = L_trt)  
summary(out_test, test = adjusted("Westfall")) # Westfall
summary(out_test) # single setp
summary(out_test, test = adjusted("holm")) # Westfall
}



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



out_test = glht(pim_parms,linfct = L_tukey)  
summary(out_test) # single setp



glm(dat_GEE)

