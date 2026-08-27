##################################################
# data generating function
##################################################
# N: total sample size, 1:1 randomization
# baseprobs: marginal probabilities for ordinal categories if all prognostic effects are zero/reference group
# covs_effects_baseline: prognostic effects on the log odds ratio scale
# time_effect: effect of time per unit on the log odds ratio scale
# time_trt_interaction: time*treatment interaction effect in the log odds ratio scale
# visits: visit times
# corMatrix: correlation matrix for within-subject dependency
#' @import simstudy
require(simstudy)
gen_data <- function(N,baseprobs,covs_effects_baseline,
                     time_effect,trt_effect,time_trt_interaction,
                     visits,corMatrix)
{
  if (!requireNamespace("simstudy", quietly = TRUE)) {
    stop("Package 'simstudy' is required but not installed.")
  }
  # calculate marginal probabilities for trt/control at each visit, without covariate effects
  probs_control_allvisits <- lapply(visits,function(visit)
  {
    # for N individuals at each visit
    shift_control <- time_effect*visit
    p0 <- calculate_multinomial_probabilities(shift_control,baseprobs)
    return(p0)
  })
  probs_trt_allvisits <- lapply(visits,function(visit)
  {
    shift_trt <- trt_effect*(visit > 0) + time_effect*visit + time_trt_interaction*visit
    p1 <- calculate_multinomial_probabilities(shift_trt,baseprobs)
    return(p1)
  })
  # define covariate distribution
  # 2 prognostic covariates
  # age, mean 60 sd 10, from SID-GBS trial
  def <- simstudy::defData(varname = "age",dist = "normal",formula = 60,variance = 10^2)
  # 40% patients have preceding diarrhea, from SID-GBS trial
  def <- simstudy::defData(def,"pre_diarrhea",dist = "binomial",formula = 0.4,variance = 1)
  # 1:1 randomization
  def <- simstudy::defData(def,varname = "trt",dist = "trtAssign",formula = "1;1")
  # relationship specification, for continuous variables, mean-centered
  def <- simstudy::defData(def,varname = "z",
                           formula = paste0(
                             ifelse(covs_effects_baseline['age']>=0,"+",""),covs_effects_baseline['age'],'*age',
                             ifelse(covs_effects_baseline['pre_diarrhea']>=0,"+",""),covs_effects_baseline['pre_diarrhea'],'*pre_diarrhea'),
                           dist="nonrandom")
  data_cov <- simstudy::genData(N,def)
  data_cov$age <- round(data_cov$age,0)
  data_trt_cov <- data_cov[data_cov$trt==1,]
  data_trt <- simstudy::genOrdCat(data_trt_cov,adjVar = "z",
                                  baseprobs = matrix(unlist(probs_trt_allvisits),byrow = TRUE,ncol = length(baseprobs)),
                                  corMatrix = corMatrix,
                                  prefix = "visit",asFactor = FALSE)
  data_control_cov <- data_cov[data_cov$trt==0,]
  data_control <- simstudy::genOrdCat(data_control_cov,adjVar = "z",
                                      baseprobs = matrix(unlist(probs_control_allvisits),byrow = TRUE,ncol = length(baseprobs)),
                                      corMatrix = corMatrix,
                                      prefix = "visit",asFactor = FALSE)
  data_wide <- rbind(data_trt,data_control)
  data_wide <- data_wide[order(data_wide$id),]
  # combine and convert to long format
  data_long <- reshape2::melt(data_wide,id.vars = c("id",names(covs_effects_baseline),"z"),measure.vars = paste0("visit",1:length(visits)),
                              variable.name = "visit_cat",value.name = "GBS_DS")
  data_long <- data_long[order(data_long$id),]
  data_long$time <- rep(visits,N)
  data_sim <- list("wide_format"=data_wide,
                   "long_format"=data_long)
  return(data_sim)
}


##################################################
# calculate_multinomial_probabilities
##################################################
#' function to calculate the new multinomial probabilities 
#' lnORs: shifted log OR at each cut-point, if assume common OR then can supply
#' a single log OR
#' p0: the starting multinomial probabilities 
#' common : whether assume common OR, if so the first element of input lnORs
#' will be repeated for every cutpoint
calculate_multinomial_probabilities <- function(lnORs,p0,common=TRUE) {
  # number of ordinal categories
  k <- length(p0)
  if(common) lnORs = rep(lnORs[1],k-1)
  # cumulative logit of p0
  logit0 <- -qlogis(cumsum(p0))[-k]
  # cumulative logit of p1 is shifted by supplied log ORs
  logit1 <- logit0 + lnORs
  # get cumulative probabilities for p1
  cump1 <- plogis(-logit1)
  # get multinomial probabilities for all categories
  p1 <- c(cump1,1)-c(0,cump1)
  return(p1)
}

#############################################
# function to calculate true win odds
#############################################
calculate_win_odds_one_visit <- function(p0,p1)
{
  k <- length(p0)
  # P(Y1 > Y0)
  Y1_gt_Y0 <- sum(sapply(1:(k-1), function(x) sum(p1[(x+1):k])*p0[x]),na.rm=T)
  # P(Y1 < Y0)
  Y1_lt_Y0 <- sum(sapply(2:k, function(x) sum(p1[1:(x-1)])*p0[x]),na.rm=T)
  # P(Y1 = Y0)
  Y1_eq_Y0 <- sum(sapply(1:k, function(x) p1[x]*p0[x]),na.rm=T)
  # win probability
  WP_true <- Y1_gt_Y0+0.5*Y1_eq_Y0
  # win odds
  WO_true <- (Y1_gt_Y0+0.5*Y1_eq_Y0)/(Y1_lt_Y0+0.5*Y1_eq_Y0)
  return(c("true win probability"=WP_true,"true win odds"=WO_true))
}
calculate_win_odds_multiple_visits <-  function(baseprobs,conditional_effect,
                                                time_effect,trt_effect,time_trt_interaction,visits)
{
  res <- sapply(visits,function(visit)
  {
    # for N individuals at each visit
    shift_control <- time_effect*visit + conditional_effect
    p0 <- calculate_multinomial_probabilities(shift_control,baseprobs)
    shift_trt <- trt_effect * (visit>0) + time_effect*visit + time_trt_interaction*visit +  conditional_effect
    p1 <- calculate_multinomial_probabilities(shift_trt,baseprobs)
    win_odds <- calculate_win_odds_one_visit(p0,p1)["true win odds"]
    return(win_odds)
  })
  return(res)
}
# use large simulated sample to marginal out the covariate distribution
calculate_adjusted_win_odds_multiple_visits <-  function(N_approx,baseprobs,covs_effects_baseline,
                                                         time_effect,trt_effect,time_trt_interaction,visits)
{
  dat_full <- gen_data(N_approx,baseprobs,covs_effects_baseline,
                       time_effect,trt_effect,time_trt_interaction,
                       visits,corMatrix=diag(length(visits)))[["long_format"]]
  dat_full$z_rounded <- round(dat_full$z,1)
  table_z <- data.frame(table(dat_full$z_rounded))
  colnames(table_z) <- c("conditional_effect","count")
  table_z$conditional_effect <- as.numeric(levels(table_z$conditional_effect))
  table_z$freq <- table_z$count/sum(table_z$count)
  wos_conditional=sapply(table_z$conditional_effect, function(i)
    calculate_win_odds_multiple_visits(baseprobs,conditional_effect=i,
                                       time_effect,trt_effect,time_trt_interaction,visits)
  )
  # weigh each conditional win odds w.r.t the strata frequency to get marginal win odds
  wos_marginal <- apply(wos_conditional, 1, function(x) t(c(x))%*%c(table_z$freq))
  return(wos_marginal)
}

##########################################
# function to generate correlation matrix
##########################################
# n_visits: number of visits
# rho: correlation coefficient
# corstr: correlation structure, either "independent","exchangeable", or "ar1"
generate_corMatrix <- function(n_visits,rho,corstr)
{
  if(corstr %in% c("independence","ind"))
  {
    diag(n_visits)
  } else if (corstr == "ar1")
  {
    exponent <- abs(matrix(1:n_visits - 1, nrow = n_visits, ncol = n_visits, byrow = TRUE) - 
                      (1:n_visits - 1))
    rho^exponent
  } else if (corstr %in% c("exchangeable","cs","compound symmetry","exch"))
  {
    mat <- matrix(rho, nrow = n_visits, ncol = n_visits)
    diag(mat) <- 1
    mat
  } else
  {
    warning("unknown correlation structure")
    return(NULL)
  }
}


###########################################################
# function to assess simulation results
###########################################################


assess_results <- function(type=NULL,N=NULL,rho=NULL,hypothesis=NULL,results=NULL)
{
  load("estimands_win_odds.RData")
  if(is.null(results))
  {
    load(paste0("C:/Users/ylong/OneDrive - LUMC/Projects/Longitudinal_win_odds/codes/results/results-",type,"-N",N,"-rho",rho,"-",hypothesis,".RData"))
  }
  if(hypothesis=="null")
  {
    true_wos <- true_wos_null
  } else if (hypothesis=="SID")
  {
    true_wos <- true_wos_SID
  } else if (hypothesis=="pos")
  {
    true_wos <- true_wos_pos
  } else
  {
    true_wos <- true_wos_con
  }
  # assess coverage of 95% CI for log win odds
  coverage <- rowMeans(sapply(results, function(i) 
  {
    i[[4]]
  }),na.rm=T)
  # power/type I error
  power_each_visit <- rowMeans(sapply(results, function(i) 
  {
    CI <- i[[3]]
    # CI does not contain 0
    CI[,1] > 0 | CI[,2] < 0
  }),na.rm=T)
  # assess bias of log win odds estimator
  est_log_wos <- sapply(results, function(i) return(i[[1]]))
  bias_log_wos <- rowMeans(est_log_wos,na.rm = TRUE) - log(true_wos)
  # convert to odds scale
  est_wos <- exp(est_log_wos)
  bias_wos <- rowMeans(est_wos,na.rm = TRUE) - true_wos
  # on the probability scale it should be unbiased
  bias_MWS <- rowMeans(est_wos/(1+est_wos),na.rm = TRUE) - true_wos/(1+true_wos)
  
  # assess estimated variance
  # "true" variance by directly calculating the variance of MC estimates
  var_log_wos_MC <- apply(est_log_wos, 1, function(x) var(x,na.rm = TRUE))
  # average estimated variance
  var_log_wos_avg <- rowMeans(sapply(results, function(i) return(i[[2]])),na.rm = TRUE)
  # global power
  power_global <- mean(sapply(results, function(i) i[[5]]<0.05))
  print(paste0("global power is: ",power_global))
  res_assessed <- data.frame(rbind(coverage,power_each_visit,true_wos,bias_log_wos,bias_wos,bias_MWS,
                                   var_log_wos_MC,var_log_wos_avg))
  colnames(res_assessed) <- paste0("visit",1:length(power_each_visit))
  return(res_assessed)
  
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
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
  utu_test = matrix(0,nrow=ncol(U),ncol=ncol(U))
  phi = ((outer(individuals_var1, individuals_var1, FUN = "==")+outer(individuals_var2, individuals_var2, FUN = "==")+outer(individuals_var1, individuals_var2, FUN = "==")+outer(individuals_var2, individuals_var1, FUN = "=="))>0)*1
  for(j in 1:nrow(U)){
    for(k in 1:nrow(U)){
      utu_test = utu_test+ (phi[j,k]*U[j,]%*%t(U[k,]))
    }
  }
  rownames(utu)=rownames(utu_test)=colnames(utu_test) = colnames(utu) = c()

  
  #skip this last part if U.diff is identity... (note: NULL here will mean identity)
  if(is.null(U.diff)) return (utu)
  
  U.diff.inv<-solve(U.diff)
  
  return(U.diff.inv%*%utu%*%U.diff.inv)
}

