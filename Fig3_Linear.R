rm(list = ls())    #delete objects
cat("\014")
library(class)
library(ggplot2)
library(dplyr)
library(glmnet)
library(alocv)
library(rmutil)
library(tictoc)
library(latex2exp)
set.seed(1)
p_values         =     c(200, 1000, 10000)
for (p in p_values){
n                =     1000
k                =     100
alpha_elnet      =     0.5
o                =     1  
m                =     30 #number of lambdas
dfmax_           =     floor(min(p,n) * 0.7)
spikeCov         =     1 # Spiked(=1) or Toeplitz(=0) covariance

if (spikeCov){
  # spike covariance
  a_             =    0.5 # correlation coeeficient
  row_           =    c(1:p, 1:p)
  column_        =    c(1:p, rep(p+1, p))
  elmts          =    c(rep(sqrt(1-a_), p), rep(sqrt(a_), p))
  F              =    sparseMatrix(i = row_, j = column_, x = elmts)
} else {
  a_               =     0.9 # Toeplitz AR(1) covariance
  F                =    matrix(rep(0, p*p), nrow = p, ncol = p)
  for (row_ in 1:p){
    for (column_ in 1:row_){
      F[row_, column_]   =   a_^abs(row_ - column_)
    }
  }
  F              =    t(F) 
}



# to make sure the var(x^T * beta.star) = 1
F              =    F / sqrt(sum(F[1, ]^2) * k) # C   =    F %*% t(F)
beta.star      =    rep(0, p)
iS             =    sample(1:p, k)
beta.star[iS]  =    rlaplace(k, m=0, s=1/sqrt(2))
X              =    F %*% matrix(rnorm( n*ncol(F), mean = 0, sd = 1 ), nrow = ncol(F), ncol = n)
X              =    t(X)
e              =    rnorm(n, mean = 0, sd = o)
y              =    X %*% beta.star + e
#fit            =    glmnet(X, y, alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, dfmax = dfmax_)
#lambdaS        =    exp(seq(log(min(fit$lambda)), log(max(fit$lambda)), length.out = m))
lambdaS        =    exp(seq(log(1/n), log(100/n), length.out = m))

ptm            =     proc.time()      
lo             =     cv.glmnet(X, y, alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, lambda = lambdaS, nfolds = n)
ptm            =     proc.time() - ptm
time.lo        =     ptm["elapsed"] 

ptm            =     proc.time() 
fit            =     glmnet(X, y, alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, lambda = lambdaS)
ptm            =     proc.time() - ptm
time.fit       =     ptm["elapsed"] 

ptm            =     proc.time()
alo            =     glmnetALO(X, y, glm_obj = fit, alpha = alpha_elnet, standardize = FALSE, type.measure = "mse")
ptm            =     proc.time() - ptm
time.alo       =     time.fit + ptm["elapsed"] 

cat(sprintf("n = %s| p = %s \n", n, p))
cat(sprintf("TIME: lo = %.2f| alo = %.2f| fit =%.2f \n", time.lo, time.alo, time.fit))
cat(sprintf(" df_max/p = %.3f \n", max(fit$df/p)))
cat(sprintf("-------------------------------------- \n"))

eror           =     data.frame(c(rep("LO", length(lo$lambda)),  rep("ALO", length(alo$lambda)) ), 
                                    n*c(lo$lambda, alo$lambda) ,
                                    c(lo$cvm, alo$alom),
                                    c(lo$cvsd, rep(0, length(alo$lambda))))
colnames(eror) =     c("method", "lambda", "err", "se")
eror.plot      =     ggplot(eror, aes(x=lambda, y = err, color=method)) +   geom_line(size=1) 
eror.plot      =     eror.plot  + scale_x_log10()#(breaks = c(seq(0.1,2.4,0.2)))   
eror.plot      =     eror.plot  + theme(legend.text = element_text(colour="black", size=16, face="bold", family = "Courier")) 
eror.plot      =     eror.plot  + geom_pointrange(aes(ymin=err-se, ymax=err+se),  size=0.8,  shape=15)
eror.plot      =     eror.plot  + theme(legend.title=element_blank()) 
eror.plot      =     eror.plot  + scale_color_discrete(breaks=c("LO", "ALO"))
eror.plot      =     eror.plot  + theme(axis.title.x = element_text(size=24),
                                        axis.text.x  = element_text(angle=0, vjust=0.5, size=14),
                                        axis.text.y  = element_text(angle=0, vjust=0.5, size=14)) 
#eror.plot      =     eror.plot  + theme(axis.title.y = element_text(size=16, face="bold", family = "Courier")) 
eror.plot      =     eror.plot  + xlab( expression(paste( lambda))) + ylab("")
eror.plot      =     eror.plot  + theme(plot.title = element_text(hjust = 0.5, vjust = -32, size=20, family = "Courier"))
#eror.plot      =     eror.plot  + ggtitle(TeX(sprintf("$n$=%s,$p$=%s,$t_{LO}$=%s,$t_{ALO}$=%0.3f,$t_{FIT}$=%.3f",n,p,time.lo,time.alo,time.fit))) 
eror.plot      =     eror.plot  + ggtitle((sprintf("n=%s, p=%s \n\n LO:%0.3f(sec) \n ALO:%0.3f(sec) \n FIT:%.3f(sec)",n,p,time.lo,time.alo,time.fit))) 


if (spikeCov){
  if (alpha_elnet == 1) {
    pdf(sprintf("/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/gaussian_lasso_spiked_n_%s_p_%s_k_%s.pdf", n,p,k))
  } else {
    pdf(sprintf("/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/gaussian_elnet_spiked_n_%s_p_%s_k_%s.pdf", n,p,k))
  }
} else {
  if (alpha_elnet == 1) {
    pdf(sprintf("/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/gaussian_lasso_toeplitz_n_%s_p_%s_k_%s.pdf", n,p,k))
  } else {
    pdf(sprintf("/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/gaussian_elnet_toeplitz_n_%s_p_%s_k_%s.pdf", n,p,k))
  }
}
print(eror.plot)
dev.off()

}
