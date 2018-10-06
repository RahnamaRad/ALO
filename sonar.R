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
set.seed(0)
# generate data
p                =     1000
delta            =     1   
n                =     delta * p
rho              =     0.1
k                =     rho * n
alpha_elnet      =     1
m                =     30

#sig_var        =   as.numeric(t(beta.star) %*% C %*% beta.star)
#C              =    C /  sig_var
#F              =    F /  sqrt(sig_var)
dfmax_         =     floor(min(n,p) * 0.7)
spikeCov       =     0

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



# to make sure the var(x^T * beta.star) = 1
F              =    F / sqrt(sum(F[1, ]^2) * k) # C   =    F %*% t(F)
beta.star      =    rep(0, p)
iS             =    sample(1:p, k)
beta.star[iS]  =    rlaplace(k, m=0, s=1/sqrt(2))
X              =    F %*% matrix(rnorm( n*ncol(F), mean = 0, sd = 1 ), nrow = ncol(F), ncol = n)
X              =    t(X)
py             =    exp(X %*% beta.star) / (1 + exp(X %*% beta.star))
y              =    rep(0, n)
for (i in 1:n){
  y[i]=rbinom(1,1,py[i])
}

fit            =    glmnet(X, y, family = "binomial", alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, dfmax = dfmax_)
lambdaS        =    exp(seq(log(min(fit$lambda)), log(max(fit$lambda)*0.9), length.out = m))

ptm            =     proc.time()      
lo             =     cv.glmnet(X, y, family = "binomial", alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, lambda = lambdaS, nfolds = n, type.measure="class")
ptm            =     proc.time() - ptm
time.lo        =     ptm["elapsed"] / (length(lo$lambda))

ptm            =     proc.time() 
fit            =     glmnet(X, y, family = "binomial", alpha = alpha_elnet,  intercept = FALSE, standardize = FALSE, lambda = lambdaS)
ptm            =     proc.time() - ptm
time.fit       =     ptm["elapsed"] / (length(fit$lambda))


ptm            =     proc.time()
alo_raw        =     rep(0, length(fit$lambda))
for (i in 1:length(fit$lambda)) {
  z_alo  = rep(0,n)
  if (fit$df[i] > 0) {
    S        =    which(fit$beta[ ,i] != 0)
    XS       =    as.matrix(X[ ,S])
    beta.hat =    fit$beta[, i]
    ez       =    exp(XS %*% beta.hat[S])
    ld       =    ez/(1+ez) - y
    ldd      =    ez/(1+ez)^2
    diag_ldd =    sparseMatrix(i = 1:n , j = 1:n, x = as.numeric(ldd))
    J        =    t(XS) %*% diag_ldd %*% XS + fit$lambda[i] * n * (1-alpha_elnet) * sparseMatrix(i = 1:fit$df[i] , j = 1:fit$df[i], x = rep(1,fit$df[i]))
    H        =    XS %*% solve(J, t(XS)) %*% diag_ldd
    z_alo    =    XS %*% beta.hat[S] + (ld/ldd) * diag(H) / (1-diag(H)) 
  }
  alo_raw[i] =  mean(as.numeric(exp(z_alo)/(1+exp(z_alo))>0.5) != y)  #2*mean(-z_alo*y + log(1+exp(z_alo)))
    #mean(as.numeric(exp(z_alo)/(1+exp(z_alo))>0.5) != y)      
}
ptm            =     proc.time()   -    ptm
time.alo_raw   =     time.fit      +    ptm["elapsed"] / (length(fit$lambda))



ptm            =     proc.time()
alo            =     glmnetALO(X, y, glm_obj = fit, alpha = alpha_elnet, standardize = FALSE, type.measure = "class")
ptm            =     proc.time() - ptm
time.alo       =     time.fit + ptm["elapsed"] / (length(alo$lambda))

cat(sprintf("n = %s| p = %s \n", n, p))
cat(sprintf("TIME: lo = %.2f| alo = %.2f| fit =%.2f \n", time.lo, time.alo, time.fit))
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
eror.plot      =     eror.plot  + theme(axis.title.x = element_text(size=24)) 
#eror.plot      =     eror.plot  + theme(axis.title.y = element_text(size=16, face="bold", family = "Courier")) 
eror.plot      =     eror.plot  + xlab( expression(paste( lambda))) + ylab("")
eror.plot      =     eror.plot  + theme(plot.title = element_text(hjust = 0.5, vjust = -10, size=20, face="bold",  family = "Courier"))
eror.plot      =     eror.plot  + ggtitle(TeX(sprintf("$n$=%s, $p$=%s, $k=%s$",n,p,k)) ) 

if (spikeCov){
  if (alpha_elnet == 1) {
    pdf(sprintf("/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/logistic_lasso_spiked_n_%s_p_%s_k_%s.pdf", n,p,k))
  } else {
    pdf(sprintf("/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/logistic_elnet_spiked_n_%s_p_%s_k_%s.pdf", n,p,k))
  }
} else {
  if (alpha_elnet == 1) {
    pdf(sprintf("/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/logistic_lasso_toeplitz_n_%s_p_%s_k_%s.pdf", n,p,k))
  } else {
    pdf(sprintf("/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/logistic_elnet_toeplitz_n_%s_p_%s_k_%s.pdf", n,p,k))
  }
}
print(eror.plot)
dev.off()



cat(sprintf("n = %s| p = %s \n", n, p))
cat(sprintf("TIME: lo = %.2f| alo = %.2f| alo_raw = %.2f| fit =%.2f \n", time.lo, time.alo, time.alo_raw, time.fit))
cat(sprintf("-------------------------------------- \n"))

eror           =     data.frame(c(rep("LO", length(lo$lambda)),  rep("ALO", length(alo$lambda)),  rep("ALO(raw)", length(fit$lambda)) ), 
                                n*c(lo$lambda, alo$lambda, fit$lambda) ,
                                c(lo$cvm, alo$alom, alo_raw),
                                c(lo$cvsd, rep(0, length(alo$lambda)), rep(0, length(fit$lambda))))
colnames(eror) =     c("method", "lambda", "err", "se")
eror.plot      =     ggplot(eror, aes(x=lambda, y = err, color=method)) +   geom_line(size=1) 
eror.plot      =     eror.plot  + scale_x_log10()#(breaks = c(seq(0.1,2.4,0.2)))   
eror.plot      =     eror.plot  + theme(legend.text = element_text(colour="black", size=16, face="bold", family = "Courier")) 
eror.plot      =     eror.plot  + geom_pointrange(aes(ymin=err-se, ymax=err+se),  size=0.8,  shape=15)
eror.plot      =     eror.plot  + theme(legend.title=element_blank()) 
eror.plot      =     eror.plot  + scale_color_discrete(breaks=c("LO", "ALO", "ALO(raw)"))
eror.plot      =     eror.plot  + theme(axis.title.x = element_text(size=24)) 
#eror.plot      =     eror.plot  + theme(axis.title.y = element_text(size=16, face="bold", family = "Courier")) 
eror.plot      =     eror.plot  + xlab( expression(paste( lambda))) + ylab("")
eror.plot      =     eror.plot  + theme(plot.title = element_text(hjust = 0.5, vjust = -10, size=20, face="bold",  family = "Courier"))
eror.plot      =     eror.plot  + ggtitle(TeX(sprintf("$n$=%s, $p$=%s$",n,p)) ) 

print(eror.plot)

