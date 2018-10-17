rm(list = ls())    #delete objects
cat("\014")
library(class)
library(ggplot2)
library(dplyr)
library(glmnet)
set.seed(0)

# generate data
p              =     1000
n              =     250
k              =     50
beta.star      =     rep(0, p)
beta.star[1:k] =     sqrt(5*10/9)/sqrt(k)#(2*rbinom(k,1,0.5)-1) 
o              =     sqrt(2)#0.7 * sqrt(k/n)
MCMCsamples    =     500
m              =     20


for (i in 1:MCMCsamples){
  X              =     matrix(rnorm(n*p, mean = 0, sd = 1), ncol = p, nrow = n)
  e              =     rnorm(n, mean = 0, sd = o)
  y              =     X %*% beta.star + e
  if (i==1){
    cv.lasso.3f    =     cv.glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, dfmax = floor(0.99 * n), nfolds = 3)
    cv.lasso.5f    =     cv.glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, dfmax = floor(0.99 * n), nfolds = 5)
    cv.lasso.10f   =     cv.glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, dfmax = floor(0.99 * n), nfolds = 10)
    cv.lasso.lo    =     cv.glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, dfmax = floor(0.99 * n), nfolds = n)
    lasso.fit      =     glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, dfmax = floor(0.99 * n))
    lambda.min     =     max(min(lasso.fit$lambda), min(cv.lasso.3f$lambda), min(cv.lasso.5f$lambda), min(cv.lasso.10f$lambda), min(cv.lasso.lo$lambda))
    lambda.max     =     min(max(lasso.fit$lambda), max(cv.lasso.3f$lambda), max(cv.lasso.5f$lambda), max(cv.lasso.10f$lambda), max(cv.lasso.lo$lambda))
    lambda         =     rev(exp(seq(log(lambda.min), log(lambda.max), length.out = m)))
    cv.3f          =     matrix(rep(0, m*MCMCsamples), nrow = m, ncol = MCMCsamples)
    cv.5f          =     matrix(rep(0, m*MCMCsamples), nrow = m, ncol = MCMCsamples)
    cv.10f         =     matrix(rep(0, m*MCMCsamples), nrow = m, ncol = MCMCsamples)
    cv.lo          =     matrix(rep(0, m*MCMCsamples), nrow = m, ncol = MCMCsamples)
    extraErr       =     matrix(rep(0, m*MCMCsamples), nrow = m, ncol = MCMCsamples)
  }
  cv.lasso.3f     =     cv.glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, lambda = lambda, nfolds = 3)
  cv.lasso.5f     =     cv.glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, lambda = lambda, nfolds = 5)
  cv.lasso.10f    =     cv.glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, lambda = lambda, nfolds = 10)
  cv.lasso.lo     =     cv.glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, lambda = lambda, nfolds = n)
  lasso.fit       =     glmnet(X, y, alpha = 1,  intercept = FALSE, standardize = FALSE, lambda = lambda)
  
  # note that the cv.glmnet doesn't necessarily calculate CV over the values of lambdas given to the function.
  # to test this claim compare the set lambdas feeded to the cv.glmnet function and the lambdas of the object provided by cv.glmnet.
  # this creates a problem because in each iteration we have the CV evaluate over a slightly different set of lambdas
  # in order to fix this problem we use a very flexible smoothing spline to evaluate the CV over a fixed set of lambdas
  # note that since the df=number of lambdas -1, the amount of smoothing is minimal and won't change our results in a meaningful manner
  cv.3f.sspline   =     smooth.spline(cv.lasso.3f$lambda, cv.lasso.3f$cvm, df = length(cv.lasso.3f$lambda)-1) 
  cv.3f.i         =     predict(cv.3f.sspline, lambda)
  cv.3f[ ,i]      =     cv.3f.i$y
    
  cv.5f.sspline   =     smooth.spline(cv.lasso.5f$lambda, cv.lasso.5f$cvm, df = length(cv.lasso.5f$lambda)-1) 
  cv.5f.i         =     predict(cv.5f.sspline, lambda)
  cv.5f[ ,i]      =     cv.5f.i$y
  
  cv.10f.sspline  =     smooth.spline(cv.lasso.10f$lambda, cv.lasso.10f$cvm, df = length(cv.lasso.10f$lambda)-1) 
  cv.10f.i        =     predict(cv.10f.sspline, lambda)
  cv.10f[ ,i]     =     cv.10f.i$y
  
  cv.lo.sspline   =     smooth.spline(cv.lasso.lo$lambda, cv.lasso.lo$cvm, df = length(cv.lasso.lo$lambda)-1)
  cv.lo.i         =     predict(cv.lo.sspline, lambda)
  cv.lo[ ,i]      =     cv.lo.i$y
  
  extraErr[ ,i]   =     colSums((lasso.fit$beta - beta.star)^2) + o^2
  print(i)
}

cv.3f.mean     =     rowMeans(cv.3f)
cv.3f.se       =     sqrt(apply(cv.3f, 1, var))/sqrt(MCMCsamples)
cv.5f.mean     =     rowMeans(cv.5f)
cv.5f.se       =     sqrt(apply(cv.5f, 1, var))/sqrt(MCMCsamples)
cv.10f.mean    =     rowMeans(cv.10f)
cv.10f.se      =     sqrt(apply(cv.10f, 1, var))/sqrt(MCMCsamples)
cv.lo.mean     =     rowMeans(cv.lo)
cv.lo.se       =     sqrt(apply(cv.lo, 1, var))/sqrt(MCMCsamples)
extraErr.mean  =     rowMeans(extraErr)
extraErr.se    =     sqrt(apply(extraErr, 1, var))/sqrt(MCMCsamples)

error          =     data.frame(c( rep("5 fold CV", m), rep("LO", m),   rep("extraErr", m) ), n*c(lambda, lambda, lambda) ,c(cv.5f.mean, cv.lo.mean, extraErr.mean),c(cv.5f.se, cv.lo.se, extraErr.se)*sqrt(MCMCsamples))
colnames(error)=     c("method", "lambda", "err", "sd")
error.plot     =     ggplot(error, aes(x=lambda, y = err, color=method)) +   geom_line(size=0.5) + scale_x_log10() + theme(legend.text = element_text(colour="black", size=10, face="bold")) + geom_pointrange(aes(ymin=err-sd, ymax=err+sd),  size=0.4,  shape=20) +  theme(legend.title=element_blank())


eror           =     data.frame(c( rep("3 fold CV", m), rep("5 fold CV", m), rep("10 fold CV", m), rep("LO", m),   rep("extraErr", m) ), n*c(lambda, lambda, lambda, lambda, lambda) ,c(cv.3f.mean, cv.5f.mean, cv.10f.mean, cv.lo.mean, extraErr.mean),c(cv.3f.se, cv.5f.se, cv.10f.se, cv.lo.se, extraErr.se))
colnames(eror) =     c("method", "lambda", "err", "se")
#eror           =     read.table("eror.txt")

eror.plot      =     ggplot(eror, aes(x=lambda, y = err, color=method)) +   geom_line(size=0.5) 
eror.plot      =     eror.plot  + scale_x_log10(breaks = c(seq(10,200,10)))   
eror.plot      =     eror.plot  + theme(legend.text = element_text(colour="black", size=12, face="bold", family = "Courier")) 
eror.plot      =     eror.plot  + geom_pointrange(aes(ymin=err-se, ymax=err+se),  size=0.4,  shape=15)
eror.plot      =     eror.plot  + theme(legend.title=element_blank()) 
eror.plot      =     eror.plot  + scale_color_discrete(breaks=c("3 fold CV","5 fold CV","10 fold CV","LO","extraErr"))
eror.plot      =     eror.plot  + theme(axis.title.x = element_text(size=18)) 
eror.plot      =     eror.plot  + theme(axis.title.y = element_text(size=16, face="bold", family = "Courier")) 
eror.plot      =     eror.plot  + xlab( expression(paste( lambda))) + ylab("error")
eror.plot      =     eror.plot  + ggtitle("Out-of-sample error \n versus \n k-fold cross validation")
eror.plot      =     eror.plot  + theme(plot.title = element_text(hjust = 0.5, vjust = -10,size=14, face="bold",  family = "Courier"))
eror.plot


write.table(eror, "eror.txt", sep="\t")
