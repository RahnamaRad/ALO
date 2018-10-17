rm(list = ls())    #delete objects
cat("\014")
library(class)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(glmnet)
library(alocv)
library(rmutil)
library(tictoc)
library(scales)
set.seed(1)
# generate data
p_               =     seq(500, 5500, 1000)
delta            =     0.1
rho              =     0.1
alpha_elnet      =     0.5
MCMCsamples      =     5
m                =     10
error.plot       =     list()

time.lo.mean     =     rep(0, length(p_))
time.alo.mean    =     rep(0, length(p_))
time.fit.mean    =     rep(0, length(p_))

time.lo.se       =     rep(0, length(p_))
time.alo.se      =     rep(0, length(p_))
time.fit.se      =     rep(0, length(p_))


spikeCov         =     1

for (i in 1:length(p_)){
  time.lo.smp    =     rep(0, MCMCsamples)
  time.alo.smp   =     rep(0, MCMCsamples)
  time.fit.smp   =     rep(0, MCMCsamples)
  p              =     p_[i] 
  n              =     p * delta
  k              =     rho * n
  #dfmax_         =     floor(0.5* min(n,p) )
  if (spikeCov){
    # spike covariance
    a_             =    0.5 
    row_           =    c(1:p, 1:p)
    column_        =    c(1:p, rep(p+1, p))
    elmts          =    c(rep(sqrt(1-a_), p), rep(sqrt(a_), p))
    F              =    sparseMatrix(i = row_, j = column_, x = elmts)
  } else {
    a_               =     0.9 # AR(1) covariance
    F                =    matrix(rep(0, p*p), nrow = p, ncol = p)
    for (row_ in 1:p){
      for (column_ in 1:row_){
        F[row_, column_]   =   a_^abs(row_ - column_)
      }
    }
    F              =    t(F) 
  }
  for (s in 1:MCMCsamples){
    beta.star      =     rep(0, p)
    iS             =     sample(1:p, k)
    beta.star[iS]  =     rlaplace(k, m=0, s=1/sqrt(2))
    X              =     F %*% matrix(rnorm( n*ncol(F), mean = 0, sd = 1 ), nrow = ncol(F), ncol = n)
    X              =     t(X)
    py             =    exp(X %*% beta.star) / (1 + exp(X %*% beta.star))
    y              =    rep(0, n)
    for (ss in 1:n){
      y[ss]=rbinom(1,1,py[ss])
    }
    if (s==1){
      #lo             =     cv.glmnet(X, y, family = "binomial", alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, dfmax = dfmax_, nfolds = n)        
      #fit            =     glmnet(X, y, family = "binomial", alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, dfmax = dfmax_)
      #lambda.min     =     max(min(fit$lambda), min(lo$lambda))
      #lambda.max     =     min(max(fit$lambda), max(lo$lambda))
      #lambdaS        =     rev(exp(seq(log(lambda.min), log(lambda.max), length.out = m)))
    }
    lambdaS        =    exp(seq(log(1/n), log(100/n), length.out = m))
    
    ptm             =     proc.time()      
    lo              =     cv.glmnet(X, y, family = "binomial", alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, lambda = lambdaS, nfolds = n, type.measure = "deviance")
    ptm             =     proc.time() - ptm
    time.lo.smp[s]  =     ptm["elapsed"]
    
    ptm             =     proc.time() 
    fit             =     glmnet(X, y, family = "binomial", alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, lambda = lambdaS)
    ptm             =     proc.time() - ptm
    time.fit.smp[s] =     ptm["elapsed"]
    
    ptm             =     proc.time()
    alo             =     glmnetALO(X, y, glm_obj = fit, alpha = alpha_elnet, standardize = FALSE, type.measure = "deviance")
    ptm             =     proc.time() - ptm
    time.alo.smp[s] =     time.fit.smp[s] + ptm["elapsed"]
  }
  time.fit.mean[i]    =  mean(time.fit.smp)
  time.lo.mean[i]     =  mean(time.lo.smp)
  time.alo.mean[i]    =  mean(time.alo.smp)
  
  time.fit.se[i]    =  sd(time.fit.smp)/sqrt(MCMCsamples)
  time.lo.se[i]     =  sd(time.lo.smp)/sqrt(MCMCsamples)
  time.alo.se[i]    =  sd(time.alo.smp)/sqrt(MCMCsamples)
  
  cat(sprintf("n = %s| p = %s \n", n, p))
  cat(sprintf("TIME: lo = %.2f ± %.2f| alo = %.2f± %.2f| fit =%.2f± %.2f \n", time.lo.mean[i], time.lo.se[i], time.alo.mean[i], time.alo.se[i], time.fit.mean[i], time.fit.se[i]))
  cat(sprintf("-------------------------------------- \n"))
}

time           =     data.frame(c(rep("LO", length(p_)), rep("ALO", length(p_)), rep("FIT", length(p_)))
                                , c(p_, p_, p_)
                                , c(time.lo.mean, time.alo.mean, time.fit.mean)
                                , c(time.lo.se, time.alo.se, time.fit.se) )
colnames(time) =     c("method", "p", "time", "se")
time.plot      =     ggplot(time, aes(x=p, y = time, color=method)) +   geom_line(size=0.5) 
time.plot      =     time.plot  + theme(legend.text = element_text(colour="black", size=12, face="bold", family = "Courier")) 
time.plot      =     time.plot  + geom_pointrange(aes(ymin=time-se, ymax=time+se),  size=0.4,  shape=15)
time.plot      =     time.plot  + theme(legend.title=element_blank()) 
time.plot      =     time.plot  + scale_color_discrete(breaks=c("LO","ALO","FIT"))
time.plot      =     time.plot  + theme(axis.title.x = element_text(size=16, family = "Courier"),
                                        axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
                                        axis.text.y  = element_text(angle=0, vjust=0.5, size=16))
time.plot      =     time.plot  + theme(axis.title.y = element_text(size=16,  family = "Courier")) 
time.plot      =     time.plot  + xlab("p=number of predictors") + ylab("time(sec)")
#time.plot      =     time.plot  + ggtitle("computational complexity")
time.plot      =     time.plot  + theme(plot.title = element_text(hjust = 0.5, vjust = -10,size=14, face="bold",  family = "Courier"))
time.plot      =     time.plot  + annotation_logticks() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                                        labels = trans_format("log10", math_format(10^.x)))
time.plot      =     time.plot  +scale_x_log10(breaks = p_) 
print(time.plot)



