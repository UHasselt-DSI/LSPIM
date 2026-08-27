# Type 1 error for treatment effect in smaller samples
# Compare PIM with 3GEE approach as in our first paper

# We fit model 3 (gives same parameter estimates as compared to model 2)

library(multcomp)
library(nparLD)


## Sim
head(design)
str(design)
design$visit_cat= paste0("visit",time)

design <- design |>
  mutate(w2 = 1*(time==2),w3 = 1*(time==3),w4 = 1*(time==4),w5 = 1*(time==5),w6 = 1*(time==6))

mod_lwo <- lwo(Y ~ group + 
                 group:w2 + group:w3 + group:w4 + group:w5 +group:w6,
               data = design,
               id = "id",
               visit = "visit_cat",
               time.varname = paste0("w",2:6),
               corstr = "independence")

summary(mod_lwo)
store_interaction_tests[l,]
L_inter <- matrix(c(0,1,0,0,0,0,
                    0,0,1,0,0,0,
                    0,0,0,1,0,0,
                    0,0,0,0,1,0,
                    0,0,0,0,0,1),nrow = 5,byrow = T)
out
W_2_inter <- t(L_inter%*%mod_lwo$coefficients)%*%solve(L_inter%*%mod_lwo$var%*%t(L_inter))%*%(L_inter%*%mod_lwo$coefficients)
p_val_inter_null <- pchisq(q=W_2_inter,df=5,lower.tail = FALSE)
pval_inter_LWO=c(pval_inter_LWO,p_val_inter_null)

nrow(L_new)
trans_matrix <- matrix(c(1,0,0,0,0,0,
                         1,1,0,0,0,0,
                         1,0,1,0,0,0,
                         1,0,0,1,0,0,
                         1,0,0,0,1,0,
                         1,0,0,0,0,1),nrow = 6,byrow = T)

# estimated log win odds at each visit
est_log_wos <- trans_matrix%*%c(mod_lwo$coefficients)
# estimated variance for the log win odds
var_log_wos <- diag(trans_matrix%*%mod_lwo$var%*%t(trans_matrix))

cbind(plogis(est_log_wos),
      pnorm(coef(mod2_adj))[3:8])


###
###

library(pim)
library(tidyverse)
library(aod)
library(geessbin)

nr_of_subjects=16
use_times = nr_of_times = 4
id = rep(1:nr_of_subjects,each=nr_of_times)
time = rep(1:nr_of_times,nr_of_subjects)
group = rep(rep(c(0,1),each=nr_of_subjects/2),each=nr_of_times)

beta0=0
beta1=sqrt(2)*0.1
beta2=sqrt(2)*0.1
beta3 =sqrt(2)*0

library(mvtnorm)

design  = data.frame(id,time,group)

nr_of_iters = 10000

L <- rbind(
  c(1, -1,rep(0,4))
)

L_new <- rbind(
  c(0,0,1,-1,0,0,0,0),
  c(0,0,1,0,-1,0,0,0),
  c(0,0,1,0,0,-1,0,0),
  c(0,0,1,0,0,0,-1,0),
  c(0,0,1,0,0,0,0,-1)
)

L_new <- rbind(
  c(0,0,1,-1,0,0),
  c(0,0,1,0,-1,0),
  c(0,0,1,0,0,-1)
)


L_new2 <- rbind(
  c(0, 0, 5/6, -1/6,-1/6, -1/6, -1/6, -1/6),
  c(0, 0, -1/6, 5/6,-1/6, -1/6, -1/6, -1/6),
  c(0, 0, -1/6, -1/6, 5/6, -1/6, -1/6, -1/6),
  c(0, 0, -1/6, -1/6, -1/6, 5/6, -1/6, -1/6),
  c(0, 0, -1/6, -1/6, -1/6, -1/6, 5/6, -1/6)
)

L_new2 = L_test = matrix(c(c(0,0,c(3/4,rep(-1/4,3))),
                  c(0,0,c(-1/4,3/4,rep(-1/4,2))),
                  c(0,0,c(rep(-1/4,2),3/4,-1/4)),
                  c(0,0,c(rep(-1/4,3),3/4))),nrow=4,byrow = TRUE)



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

for(l in 7514:nr_of_iters){
  if(l%%100==0){cat("Iteration: ",l,"\n")}
  set.seed(l)
  re = rmvnorm(nr_of_subjects,c(0,0),matrix(c(0.8^2,0,0,1.5^2),nrow=2))
  design$Y = beta0+
    (beta1)*design$time+
    beta2*(design$group==1)+
    beta3*design$time*(design$group==1)+
    rnorm(nrow(design),0,1)+
    rep(re[,2],each=nr_of_times)
  
  #design$Y = round(design$Y)
  ## Model poset, defined via compare.

  # id.fac1 <- which((design[,"group"] == 0)&(design[,"time"] == 1))
  # id.nonfac1 <- which((design[,"group"] == 1)&(design[,"time"] == 1))
  # compare1 <- expand.grid(id.fac1,id.nonfac1)
  # 
  # id.fac2 <- which((design[,"group"] == 0)&(design[,"time"] == 2))
  # id.nonfac2 <- which((design[,"group"] == 1)&(design[,"time"] == 2))
  # compare2 <- expand.grid(id.fac2,id.nonfac2)
  # 
  # id.fac3 <- which((design[,"group"] == 0)&(design[,"time"] == 3))
  # id.nonfac3 <- which((design[,"group"] == 1)&(design[,"time"] == 3))
  # compare3 <- expand.grid(id.fac3,id.nonfac3)
  # 
  # id.fac4 <- which((design[,"group"] == 0)&(design[,"time"] == 4))
  # id.nonfac4 <- which((design[,"group"] == 1)&(design[,"time"] == 4))
  # compare4 <- expand.grid(id.fac4,id.nonfac4)
  #
  # id.fac5 <- which((design[,"group"] == 0)&(design[,"time"] == 5))
  # id.nonfac5 <- which((design[,"group"] == 1)&(design[,"time"] == 5))
  # compare5 <- expand.grid(id.fac5,id.nonfac5)
  #
  # id.fac6 <- which((design[,"group"] == 0)&(design[,"time"] == 6))
  # id.nonfac6 <- which((design[,"group"] == 1)&(design[,"time"] == 6))
  # compare6 <- expand.grid(id.fac6,id.nonfac6)
  #
  # compare_between = rbind(compare1,compare2,compare3,compare4,compare5,compare6)
  # compare_between = compare_between[order(compare_between[,"Var2"]),]
  # 
  # start_timepoints <- c()
  # later_timepoints <- c()
  # k=1
  # for (i in 1:(nrow(design) - 1)) {
  #   if(i%%4 != 0){
  #     start <- i           # Current time point
  #     later <- (i+1) : (4*k)  # Later time points
  #     # Store the results in vectors
  #     start_timepoints <- c(start_timepoints, rep(start, length(later)))
  #     later_timepoints <- c(later_timepoints, later)}
  #   if(i%%4 == 0){
  #     k=k+1}
  # }
  # within1 = start_timepoints
  # within2 = later_timepoints
  # 
  # compare_within = cbind(within1,within2)
  # colnames(compare_within) = c("Var1","Var2")
  # compare_within = compare_within[order(compare_within[,2]),]
  # compare = rbind(compare_between,compare_within)
  
  
  individuals_var1=design[compare$Var1,"id"]
  individuals_var2=design[compare$Var2,"id"]

  # Model 2
  assignInNamespace("sandwich.estimator", sandwich.estimator, ns = "pim")
  mod2_orig = pim(Y~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) +I((R(group)-L(group))*(R(time)==4)*(L(time)==4)),data=design,compare=compare,link="probit")
  store_parms_mod2[l,] = coef(mod2_orig)
  store_se_mod2[l,] = diag(vcov(mod2_orig))
  response(mod2_orig)
  model.matrix(mod2_orig)
  assignInNamespace("sandwich.estimator", sandwich.estimator_clustered_mod1_mod2, ns = "pim")
  
  mod2_adj = pim(Y~I((R(time) - L(time))*R(group)*L(group))+I((R(time) - L(time))*(1-R(group))*(1-L(group)))+I((R(group)-L(group))*(R(time)==1)*(L(time)==1))+I((R(group)-L(group))*(R(time)==2)*(L(time)==2))+I((R(group)-L(group))*(R(time)==3)*(L(time)==3)) +I((R(group)-L(group))*(R(time)==4)*(L(time)==4)),data=design,compare=compare,link="probit")
  store_parms_mod2_adj[l,] = coef(mod2_adj)
  store_se_mod2_adj[l,] = diag(vcov(mod2_adj))
  
  rownames(check) = colnames(check)=NULL
  
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
  mod1 = geessbin(y~.-1-C1-C2-C3,data=dat_GEE,id=C1, corstr = "independence",beta.method="PGEE",SE.method = "FW")
  dat_GEE=dat_GEE[order(dat_GEE$C2),]
  mod2 = geessbin(y~.-1-C1-C2-C3,data=dat_GEE,id=C2, corstr = "independence",beta.method="PGEE",SE.method = "FW")
  dat_GEE=dat_GEE[order(dat_GEE$C3),]
  mod3 = geessbin(y~.-1-C1-C2-C3,data=dat_GEE,id=C3, corstr = "independence",beta.method="PGEE",SE.method = "FW")
  
  store_se_gee[l,] = diag(mod1$covb+mod2$covb-mod3$covb)
  store_parms_gee[l,] = apply(rbind(coef(mod1),coef(mod2),coef(mod3)),2,mean)
  
  
  ### Interaction tests based on adjusted PIM
  
  V=vcov(mod2_adj)
  #V = V+rnorm(length(V),0.0001,0.00001)
  b = coef(mod2_adj)
  
  ### Standard degrees of freedom
  
  # diag(L_new%*%vcov(mod2_adj)%*%t(L_new))/diag(L_new%*%vcov(mod2_orig)%*%t(L_new))
  # 
  # L_new%*%coef(mod2_adj)/sqrt(diag(L_new%*%vcov(mod2_adj)%*%t(L_new)))
  
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
  
  # while(is.null(out)){
  #   V = V+rnorm(length(V),0.0001,0.00001)
  #   out = tryCatch(wald.test(V, b, Terms = NULL, L = L_new, H0 = NULL,
  #                            df = nr_of_subjects-8, verbose = FALSE), error = function(e) {
  #                              return(NULL)
  #                            })
  # }
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


save.image("LMM_4_26.Rdata")


length(keep)
keep = which((abs(store_parms_mod2[,3])<5)&(abs(store_parms_mod2[,4])<5)&(abs(store_parms_mod2[,5])<5)&(abs(store_parms_mod2[,6])<5)&
               (abs(store_parms_mod2[,7])<5)&(abs(store_parms_mod2[,8])<5))
boxplot(store_parms_mod2[keep,])


type1_interact = apply(store_interaction_tests[keep,],2,function(x)mean(x<0.05,na.rm=TRUE))

apply(store_parms_mod2[keep,],2,mean,na.rm=TRUE)

cbind(apply(store_se_gee[keep,],2,mean,na.rm=TRUE),
      apply(store_se_mod2[keep,],2,mean,na.rm=TRUE),
      apply(store_se_mod2_adj[keep,],2,mean,na.rm=TRUE),
      apply(store_parms_mod2[keep,] ,2,var,na.rm=TRUE),
      apply(store_parms_gee[keep,] ,2,var,na.rm=TRUE))


save.image("small_sample_type1_pim_gee_26_correct_PGEE_score_unrounded_testwithmeanhypothesis.Rdata")
load("small_sample_type1_pim_gee_20_correct_PGEE_score_unrounded_testwithmeanhypothesis.Rdata")

rbind(type1_interact[1:4],
      type1_interact[5:8],
      type1_interact[9:12],
      type1_interact[13:16])


plogis(apply(store_parms_gee[keep,] ,2,mean,na.rm=TRUE))
pnorm(apply(store_parms_mod2[keep,] ,2,mean,na.rm=TRUE))
getwd()
qlogis(pnorm(truth))


out_interact = data.frame(rbind(type1_interact[1:4],
                                type1_interact[5:8],
                                type1_interact[9:12],
                                type1_interact[13:16]))
names(out_interact) = c('chisquare','ddf_L',"ddf_BW","ddf_containment")
out_interact$contrast = c("Group","Time","Group","Time")
out_interact$method = c('PIM','PIM',"GEE","GEE")



sum(store_interaction_tests[,11]<0.05,na.rm = TRUE)/sum(!is.na(store_interaction_tests[,11]))

truth = beta2/sqrt(2*(1.5^2+1))

apply(store_parms_gee,2,mean,na.rm=TRUE)[3:8] - qlogis(pnorm(truth))





sub1 = subset(design,group==0)
sub2 = subset(design,group==1)


sub1_1 = subset(sub1,time==1)
sub2_1 = subset(sub2,time==1)

use1 = sample(1:nrow(sub1_1),1000000,replace = TRUE)
use1 = sample(1:nrow(sub1_1),10000000000,replace = TRUE)


for(samplesizes_fin1 in c(16,18,20,24,30,40,50,60)){
  load( paste0("small_sample_type1_pim_gee_",samplesizes_fin1,"_correct_PGEE_score_unrounded_testwithmeanhypothesis.Rdata"))
  truth = c(rep(beta2/sqrt(2*(1)),2),rep(beta2/sqrt(2*(1+1.5^2)),6))
  bias_PIM = apply(t(t(store_parms_mod2) - truth),2,mean)
  bias_GEE = apply(t(t(qnorm(plogis(store_parms_gee))) - truth),2,mean,na.rm=TRUE)
  
  emp_PIM = apply(store_parms_mod2,2,var,na.rm=TRUE)
  emp_PIM_adjusted = apply(store_parms_mod2_adj,2,var,na.rm=TRUE)
  var_PIM_adjusted = apply(store_se_mod2_adj,2,mean,na.rm=TRUE)
  var_PIM = apply(store_se_mod2,2,mean,na.rm=TRUE)
  emp_GEE = apply(store_parms_gee,2,var,na.rm=TRUE)
  var_GEE = apply(store_se_gee,2,mean,na.rm=TRUE)
  
  type1_interact = apply(store_interaction_tests,2,function(x)mean(x<0.05,na.rm=TRUE))
  keep_bw = type1_interact[seq(3,15,4)]
  
  score_inter = c(mean(score_inter_group,na.rm=TRUE),mean(score_inter_time))
  LDA_inter = apply(store_interaction_tests[,c(17,18)],2,function(x)mean(x<0.05,na.rm=TRUE))
  if(samplesizes_fin1==20){
    bias_out_fin1 = data.frame("bias" = c(bias_PIM,bias_GEE),"method"=c(rep("PIM",8),rep("GEE",8)))
    bias_out_fin1$sampsize=samplesizes_fin1
    bias_out_fin1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),2)
    interact_out_fin1 = data.frame("TypeI" = c(keep_bw,score_inter,LDA_inter),"sampsize"=samplesizes_fin1,"method"=c("PIM","PIM","GEE","GEE","GEE_score","GEE_score","LDA","LDA"),"contrast" = c(rep(c("Group","Time"),3),"W", "A"))
    var_out_fin1 = data.frame("model_variance" = c(var_PIM,var_PIM_adjusted,var_GEE),"empirical_variance" = c(emp_PIM ,emp_PIM_adjusted ,emp_GEE),"method"=c(rep("PIM",8),rep("PIM_adjusted",8),rep("GEE",8)))
    var_out_fin1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),3)
    var_out_fin1$sampsize=samplesizes_fin1
  }
  if(samplesizes_fin1>20){
    bias_out_sub1 = data.frame("bias" = c(bias_PIM,bias_GEE),"method"=c(rep("PIM",8),rep("GEE",8)))
    bias_out_sub1$sampsize=samplesizes_fin1
    bias_out_sub1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),2)
    bias_out_fin1 = rbind(bias_out_fin1,bias_out_sub1)
    interact_out_sub1 =  data.frame("TypeI" = c(keep_bw,score_inter,LDA_inter),"sampsize"=samplesizes_fin1,"method"=c("PIM","PIM","GEE","GEE","GEE_score","GEE_score","LDA","LDA"),"contrast" = c(rep(c("Group","Time"),3),"W", "A"))
    interact_out_fin1 = rbind(interact_out_fin1,interact_out_sub1 )
    var_out1_sub1 = data.frame("model_variance" = c(var_PIM,var_PIM_adjusted,var_GEE),"empirical_variance" = c(emp_PIM ,emp_PIM_adjusted ,emp_GEE),"method"=c(rep("PIM",8),rep("PIM_adjusted",8),rep("GEE",8)))
    var_out1_sub1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),3)
    var_out1_sub1$sampsize=samplesizes_fin1
    var_out_fin1 = rbind(var_out_fin1,var_out1_sub1)
  }
}


ggplot(bias_out_fin1%>%filter(par%in%paste0("T",1:2)),aes(x=sampsize,y=bias,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0,slope=0))+ylim(min(bias_out_fin1$bias),max(bias_out_fin1$bias))+ggtitle("Bias of time coefficients")
ggplot(bias_out_fin1%>%filter(par%in%paste0("A",1:6)),aes(x=sampsize,y=bias,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0,slope=0))+ggtitle("Bias of group coefficients")
ggplot(interact_out_fin1%>%filter(method!="LDA"),aes(x=sampsize,y=TypeI,color=contrast,linetype=method))+geom_line()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Type I error for test interaction (df = nr_subjects - rank of L)")

ggplot(var_out_fin1,aes(x=sampsize,y=model_variance/empirical_variance,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Average model variance divided by empirical variance")



ggplot(interact_out_fin%>%filter(method%in%c("GEE_score","LDA")),aes(x=sampsize,y=TypeI,color=contrast,linetype=method))+geom_line()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Type I error for test interaction (df = nr_subjects - rank of L)")




var(store_parms_mod2)



pnorm(coef(mod2_orig)[3:8])





ggplot(var_out1%>%filter(par=="T1"),aes(x=sampsize,y=model_variance/empirical_variance,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Average model variance divided by empirical variance")



mean(score_inter_time)


load("small_sample_type1_pim_gee_20_correct_PGEE_unrounded.Rdata")


0.1919 0.0570 0.1075 0.0385
0.3380 0.0400 0.1075 0.0385
boxplot(store_parms_mod2)

store_interaction_tests[l,seq(3,15,4)]
out


which(score_inter_group==TRUE)




mean(store_interaction_tests[,11])
mod_use$df.residual = nr_of_subjects-qr(L_new)$rank-1
out = summary(glht(mod_use, linfct = L_new), test = Ftest())
out$test$pvalue

mean( (score_inter_group_global<0.05) & score_inter_group,na.rm=TRUE )
mean(score_inter_group,na.rm = TRUE)
mean(score_inter_group_global<0.05)

out = summary(glht(mod_use, linfct = L_new), test = Chisqtest())



mean((score_inter_time_global<0.05) & score_inter_time,na.rm=TRUE )
mean(score_inter_time,na.rm = TRUE)
mean(score_inter_time_global<0.05)

mean((score_inter_group_global<0.05) & score_inter_group,na.rm=TRUE )
mean(score_inter_group,na.rm = TRUE)
mean(score_inter_group_global<0.05)




L_test = matrix(c(c(0,0,c(5/6,rep(-1/6,5))),
                  c(0,0,c(-1/6,5/6,rep(-1/6,4))),
                  c(0,0,c(rep(-1/6,2),5/6,rep(-1/6,3))),
                  c(0,0,c(rep(-1/6,3),5/6,rep(-1/6,2))),
                  c(0,0,c(rep(-1/6,4),5/6,-1/6)),
                  c(0,0,c(rep(-1/6,5),5/6))),nrow=6,byrow = TRUE)


out = summary(glht(mod_use, linfct = L_new), test = adjusted("holm"))
out = summary(glht(mod_use, linfct = L_test), test = adjusted("holm"))

score_inter_group = c(score_inter_group,!prod(out$test$pvalues>0.05))

which(score_inter_group)
l=28
``






for(samplesizes_fin1 in c(16,18,20,24,26,30,40,50,60)){
  load( paste0("small_sample_type1_pim_gee_",samplesizes_fin1,"_correct_PGEE_score_unrounded_testwithmeanhypothesis.Rdata"))
  truth = c(rep(beta2/sqrt(2*(1)),2),rep(beta2/sqrt(2*(1+1.5^2)),6))
  bias_PIM = apply(t(t(store_parms_mod2) - truth),2,mean)
  bias_GEE = apply(t(t(qnorm(plogis(store_parms_gee))) - truth),2,mean,na.rm=TRUE)
  
  emp_PIM = apply(store_parms_mod2,2,var,na.rm=TRUE)
  emp_PIM_adjusted = apply(store_parms_mod2_adj,2,var,na.rm=TRUE)
  var_PIM_adjusted = apply(store_se_mod2_adj,2,mean,na.rm=TRUE)
  var_PIM = apply(store_se_mod2,2,mean,na.rm=TRUE)
  emp_GEE = apply(store_parms_gee,2,var,na.rm=TRUE)
  var_GEE = apply(store_se_gee,2,mean,na.rm=TRUE)
  
  type1_interact = apply(store_interaction_tests,2,function(x)mean(x<0.05,na.rm=TRUE))
  keep_bw = type1_interact[seq(3,15,4)]
  
  score_inter = c(mean(score_inter_group,na.rm=TRUE),mean(score_inter_time))
  LDA_inter = apply(store_interaction_tests[,c(17,18)],2,function(x)mean(x<0.05,na.rm=TRUE))
  if(samplesizes_fin1==20){
    bias_out_fin1 = data.frame("bias" = c(bias_PIM,bias_GEE),"method"=c(rep("PIM",8),rep("GEE",8)))
    bias_out_fin1$sampsize=samplesizes_fin1
    bias_out_fin1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),2)
    interact_out_fin1 = data.frame("TypeI" = c(keep_bw,score_inter,LDA_inter),"sampsize"=samplesizes_fin1,"method"=c("PIM_wald.test","PIM_wald.test","GEE_wald.test","GEE_wald.test","GEE_glht","GEE_glht","LDA_W","LDA_A"),"contrast" = c(rep(c("Group","Time"),3),"-", "-"))
    var_out_fin1 = data.frame("model_variance" = c(var_PIM,var_PIM_adjusted,var_GEE),"empirical_variance" = c(emp_PIM ,emp_PIM_adjusted ,emp_GEE),"method"=c(rep("PIM",8),rep("PIM_adjusted",8),rep("GEE",8)))
    var_out_fin1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),3)
    var_out_fin1$sampsize=samplesizes_fin1
  }
  if(samplesizes_fin1>20){
    bias_out_sub1 = data.frame("bias" = c(bias_PIM,bias_GEE),"method"=c(rep("PIM",8),rep("GEE",8)))
    bias_out_sub1$sampsize=samplesizes_fin1
    bias_out_sub1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),2)
    bias_out_fin1 = rbind(bias_out_fin1,bias_out_sub1)
    interact_out_sub1 =  data.frame("TypeI" = c(keep_bw,score_inter,LDA_inter),"sampsize"=samplesizes_fin1,"method"=c("PIM_wald.test","PIM_wald.test","GEE_wald.test","GEE_wald.test","GEE_glht","GEE_glht","LDA_W","LDA_A"),"contrast" = c(rep(c("Group","Time"),3),"-", "-"))
    interact_out_fin1 = rbind(interact_out_fin1,interact_out_sub1 )
    var_out1_sub1 = data.frame("model_variance" = c(var_PIM,var_PIM_adjusted,var_GEE),"empirical_variance" = c(emp_PIM ,emp_PIM_adjusted ,emp_GEE),"method"=c(rep("PIM",8),rep("PIM_adjusted",8),rep("GEE",8)))
    var_out1_sub1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),3)
    var_out1_sub1$sampsize=samplesizes_fin1
    var_out_fin1 = rbind(var_out_fin1,var_out1_sub1)
  }
}


ggplot(bias_out_fin1%>%filter(par%in%paste0("T",1:2)),aes(x=sampsize,y=bias,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0,slope=0))+ylim(min(bias_out_fin1$bias),max(bias_out_fin1$bias))+ggtitle("Bias of time coefficients")
ggplot(bias_out_fin1%>%filter(par%in%paste0("A",1:6)),aes(x=sampsize,y=bias,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0,slope=0))+ggtitle("Bias of group coefficients")
ggplot(interact_out_fin1,aes(x=sampsize,y=TypeI,color=contrast,linetype=method))+geom_line()+geom_point()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Type I error for test interaction (6 timepoints)")+ geom_hline(yintercept=c(0.05-sqrt(0.05*0.95/10000)*1.96,0.05+sqrt(0.05*0.95/10000)*1.96), linetype="dashed", color = "darkgreen")+ 
  xlab("Total sample size (evenly divided in two groups)") +
  ylab("Type I error") 

ggplot(var_out_fin1,aes(x=sampsize,y=model_variance/empirical_variance,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Average model variance divided by empirical variance")


library(ggplot2)
library(tidyverse)

for(sampsize_final_use in c(16,18,20,24,26,30,40,50,60)){
  load( paste0("small_sample_type1_pim_gee_",sampsize_final_use,"_correct_PGEE_score_unrounded_testwithmeanhypothesis_4timepoints.Rdata"))
  truth = c(rep(beta2/sqrt(2*(1)),2),rep(beta2/sqrt(2*(1+1.5^2)),4))
  bias_PIM = apply(t(t(store_parms_mod2) - truth),2,mean)
  bias_GEE = apply(t(t(qnorm(plogis(store_parms_gee))) - truth),2,mean,na.rm=TRUE)
  
  emp_PIM = apply(store_parms_mod2,2,var,na.rm=TRUE)
  emp_PIM_adjusted = apply(store_parms_mod2_adj,2,var,na.rm=TRUE)
  var_PIM_adjusted = apply(store_se_mod2_adj,2,mean,na.rm=TRUE)
  var_PIM = apply(store_se_mod2,2,mean,na.rm=TRUE)
  emp_GEE = apply(store_parms_gee,2,var,na.rm=TRUE)
  var_GEE = apply(store_se_gee,2,mean,na.rm=TRUE)
  
  type1_interact = apply(store_interaction_tests,2,function(x)mean(x<0.05,na.rm=TRUE))
  keep_bw = type1_interact[seq(3,15,4)]
  
  score_inter = c(mean(score_inter_group,na.rm=TRUE),mean(score_inter_time))
  LDA_inter = apply(store_interaction_tests[,c(17,18)],2,function(x)mean(x<0.05,na.rm=TRUE))
  if(sampsize_final_use==16){
    bias_out_use = data.frame("bias" = c(bias_PIM,bias_GEE),"method"=c(rep("PIM",6),rep("GEE",6)))
    bias_out_use$sampsize=sampsize_final_use
    bias_out_use$par = rep(c(paste0("T",1:2),paste0("A",1:4)),2)
    interact_out_use = data.frame("TypeI" = c(keep_bw,score_inter,LDA_inter),"sampsize"=sampsize_final_use,"method"=c("PIM_wald.test","PIM_wald.test","GEE_wald.test","GEE_wald.test","GEE_glht","GEE_glht","LDA_W","LDA_A"),"contrast" = c(rep(c("Group","Time"),3),"-", "-"))
    var_out_use= data.frame("model_variance" = c(var_PIM,var_PIM_adjusted,var_GEE),"empirical_variance" = c(emp_PIM ,emp_PIM_adjusted ,emp_GEE),"method"=c(rep("PIM",6),rep("PIM_adjusted",6),rep("GEE",6)))
    var_out_use$par = rep(c(paste0("T",1:2),paste0("A",1:4)),3)
    var_out_use$sampsize=sampsize_final_use
  }
  if(sampsize_final_use>16){
    bias_out_sub1 = data.frame("bias" = c(bias_PIM,bias_GEE),"method"=c(rep("PIM",6),rep("GEE",6)))
    bias_out_sub1$sampsize=sampsize_final_use
    bias_out_sub1$par = rep(c(paste0("T",1:2),paste0("A",1:4)),2)
    bias_out_use = rbind(bias_out_use,bias_out_sub1)
    interact_out_sub1 =  data.frame("TypeI" = c(keep_bw,score_inter,LDA_inter),"sampsize"=sampsize_final_use,"method"=c("PIM_wald.test","PIM_wald.test","GEE_wald.test","GEE_wald.test","GEE_glht","GEE_glht","LDA_W","LDA_A"),"contrast" = c(rep(c("Group","Time"),3),"-", "-"))
    interact_out_use= rbind(interact_out_use,interact_out_sub1 )
    var_out1_sub1 = data.frame("model_variance" = c(var_PIM,var_PIM_adjusted,var_GEE),"empirical_variance" = c(emp_PIM ,emp_PIM_adjusted ,emp_GEE),"method"=c(rep("PIM",6),rep("PIM_adjusted",6),rep("GEE",6)))
    var_out1_sub1$par = rep(c(paste0("T",1:2),paste0("A",1:4)),3)
    var_out1_sub1$sampsize=sampsize_final_use
    var_out_use= rbind(var_out_use,var_out1_sub1)
  }
}


ggplot(bias_out_use%>%filter(par%in%paste0("T",1:2)),aes(x=sampsize,y=bias,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0,slope=0))+ylim(min(bias_out_use$bias),max(bias_out_use$bias))+ggtitle("Bias of time coefficients")
ggplot(bias_out_use%>%filter(par%in%paste0("A",1:6)),aes(x=sampsize,y=bias,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0,slope=0))+ggtitle("Bias of group coefficients")
ggplot(interact_out_use,aes(x=sampsize,y=TypeI,color=contrast,linetype=method))+geom_line()+geom_point()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Type I error for test interaction (4 timepoints)")+ geom_hline(yintercept=c(0.05-sqrt(0.05*0.95/10000)*1.96,0.05+sqrt(0.05*0.95/10000)*1.96), linetype="dashed", color = "darkgreen")+ 
  xlab("Total sample size (evenly divided in two groups)") +
  ylab("Type I error") 

p4a=ggplot(var_out_use,aes(x=sampsize,y=model_variance/empirical_variance,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Average model variance divided by empirical variance (4 timepoints)")








for(sampsize_final_use in c(16,18,20,24,26,30,40,50,60)){
  load( paste0("small_sample_type1_pim_gee_",sampsize_final_use,"_correct_PGEE_score_unrounded_testwithmeanhypothesis.Rdata"))
  truth = c(rep(beta2/sqrt(2*(1)),2),rep(beta2/sqrt(2*(1+1.5^2)),6))
  bias_PIM = apply(t(t(store_parms_mod2) - truth),2,mean)
  bias_GEE = apply(t(t(qnorm(plogis(store_parms_gee))) - truth),2,mean,na.rm=TRUE)
  
  emp_PIM = apply(store_parms_mod2,2,var,na.rm=TRUE)
  emp_PIM_adjusted = apply(store_parms_mod2_adj,2,var,na.rm=TRUE)
  var_PIM_adjusted = apply(store_se_mod2_adj,2,mean,na.rm=TRUE)
  var_PIM = apply(store_se_mod2,2,mean,na.rm=TRUE)
  emp_GEE = apply(store_parms_gee,2,var,na.rm=TRUE)
  var_GEE = apply(store_se_gee,2,mean,na.rm=TRUE)
  
  type1_interact = apply(store_interaction_tests,2,function(x)mean(x<0.05,na.rm=TRUE))
  keep_bw = type1_interact[seq(3,15,4)]
  
  score_inter = c(mean(score_inter_group,na.rm=TRUE),mean(score_inter_time))
  LDA_inter = apply(store_interaction_tests[,c(17,18)],2,function(x)mean(x<0.05,na.rm=TRUE))
  if(sampsize_final_use==16){
    bias_out_fin1 = data.frame("bias" = c(bias_PIM,bias_GEE),"method"=c(rep("PIM",8),rep("GEE",8)))
    bias_out_fin1$sampsize=sampsize_final_use
    bias_out_fin1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),2)
    interact_out_fin1 = data.frame("TypeI" = c(keep_bw,score_inter,LDA_inter),"sampsize"=sampsize_final_use,"method"=c("PIM","PIM","GEE","GEE","GEE_score","GEE_score","LDA","LDA"),"contrast" = c(rep(c("Group","Time"),3),"W", "A"))
    var_out_fin1 = data.frame("model_variance" = c(var_PIM,var_PIM_adjusted,var_GEE),"empirical_variance" = c(emp_PIM ,emp_PIM_adjusted ,emp_GEE),"method"=c(rep("PIM",8),rep("PIM_adjusted",8),rep("GEE",8)))
    var_out_fin1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),3)
    var_out_fin1$sampsize=sampsize_final_use
  }
  if(sampsize_final_use>16){
    bias_out_sub1 = data.frame("bias" = c(bias_PIM,bias_GEE),"method"=c(rep("PIM",8),rep("GEE",8)))
    bias_out_sub1$sampsize=sampsize_final_use
    bias_out_sub1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),2)
    bias_out_fin1 = rbind(bias_out_fin1,bias_out_sub1)
    interact_out_sub1 =  data.frame("TypeI" = c(keep_bw,score_inter,LDA_inter),"sampsize"=sampsize_final_use,"method"=c("PIM","PIM","GEE","GEE","GEE_score","GEE_score","LDA","LDA"),"contrast" = c(rep(c("Group","Time"),3),"W", "A"))
    interact_out_fin1 = rbind(interact_out_fin1,interact_out_sub1 )
    var_out1_sub1 = data.frame("model_variance" = c(var_PIM,var_PIM_adjusted,var_GEE),"empirical_variance" = c(emp_PIM ,emp_PIM_adjusted ,emp_GEE),"method"=c(rep("PIM",8),rep("PIM_adjusted",8),rep("GEE",8)))
    var_out1_sub1$par = rep(c(paste0("T",1:2),paste0("A",1:6)),3)
    var_out1_sub1$sampsize=sampsize_final_use
    var_out_fin1 = rbind(var_out_fin1,var_out1_sub1)
  }
}


ggplot(bias_out_fin1%>%filter(par%in%paste0("T",1:2)),aes(x=sampsize,y=bias,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0,slope=0))+ylim(min(bias_out_fin1$bias),max(bias_out_fin1$bias))+ggtitle("Bias of time coefficients")
ggplot(bias_out_fin1%>%filter(par%in%paste0("A",1:6)),aes(x=sampsize,y=bias,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0,slope=0))+ggtitle("Bias of group coefficients")
p3b = ggplot(interact_out_fin1,aes(x=sampsize,y=TypeI,color=contrast,linetype=method))+geom_line()+geom_point()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Type I error for test interaction (6 timepoints)")+ geom_hline(yintercept=c(0.05-sqrt(0.05*0.95/10000)*1.96,0.05+sqrt(0.05*0.95/10000)*1.96), linetype="dashed", color = "darkgreen")

p4b = ggplot(var_out_use%>%filter(method%in%c("GEE","PIM_adjusted")&!(par%in%c("T1","T2"))),aes(x=sampsize,y=model_variance/empirical_variance,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Average model variance divided by empirical variance (4 timepoints)")

ggplot(var_out_use%>%filter(method%in%c("GEE","PIM_adjusted")&(par%in%c("T1","T2"))),aes(x=sampsize,y=model_variance/empirical_variance,color=par,linetype=method))+geom_line()+geom_abline(aes(intercept=0.05,slope=0))+ggtitle("Average model variance divided by empirical variance (4 timepoints)")

library(gridExtra)
grid.arrange(p3a,p3b)
grid.arrange(p4a,p4b)
