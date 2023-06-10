rm(list=ls())

# install.packages("devtools")
# install_github("qingyuanzhao/mr.raps")

library(glmnet)
library(parallel)
library( mvtnorm)
library(AER)
library(devtools)
library(mr.raps)




n <- 4000
#Number of selected SNPs
p <- 300
#Number of SNPs affecting the exposure
np_act <- 5 #sample at 0.3 MAF HWE
sub_split <- 10 #number of split

set.seed(1)


simu_cross_fit_IV_est <- function(n,beta0,h2 )
{
  
  #Genotype data
  G <- matrix(sample(c(0,1,2),
                     size= (n*p),
                     prob = c(0.7^2, 2*0.7*0.3, 0.3^2),
                     replace = TRUE),
              nrow= n)
  
  #Hidden confounding
  H <- rnorm(n)
  #correlated noise 
  noise <-rmvnorm(n, 
                  mean = rep(0, 2), 
                  sigma = matrix(nrow = 2,
                                 byrow = TRUE,
                                 c(1,0.2,0.2,1)
                  )
  )
  
  
  
  V <- noise[,2]
  
  U <- noise[,1]
  
  
  #rescaling for insuring h2 = h2 on average
  tt <- apply( G[, (1:np_act)],1,sum)#genotype effect
  if( h2==0)
  {
    resc <- 0
  }else{resc <- ((var(H)+var(V))*(h2/(1-h2)))/var(tt)
  
  
  }
  
  ZPI <- sqrt(resc)*tt 
  
  
  #X is compound of ZPI , H hidden confounding, and a noise that correlate with the noise of Y
  
  X <- ZPI + H + V
  
  #Y the exposure alos affected by H
  
  Y <- beta0 * X + H + U
  
  build_IV_sub_sample <- function(i) #i in 0:(sub_split-1)
  {
    indx <- (1+i*(n/sub_split)):((i+1)*(n/sub_split))
    
    
    Gtemp <- G[-indx,]
    Xtemp <- X[-indx]
    Ytemp <- Y[indx]
    
    lasso <-  cv.glmnet(x=Gtemp, y=Xtemp,alpha=1)
    return(predict(lasso, G[indx,], s = "lambda.min"))
  }
  #Cross fitting LASSO
  #res  <- mclapply(0:(sub_split-1), build_IV_sub_sample, mc.cores = 6)
  res  <-  lapply(0:(sub_split-1), build_IV_sub_sample )
  
  #CFMR1
  IV_each <- function(i){
    indx <- (1+i*(n/sub_split)):((i+1)*(n/sub_split))
    Ytemp <- Y[indx]
    return(summary(lm(Ytemp~as.vector(res[[i+1]])))$coef[2,])
  }
  CFMR1_each <- lapply(0:(sub_split-1), IV_each )
  res1 <- do.call(rbind,CFMR1_each)#cross fitted instrument
  
  # CFMR2, use the prediction of X as single IV
  IV   <- do.call(c,res)#cross fitted instrument
  res2 <- summary(ivreg(Y~X|IV))$coef[2,]
  
  # regression of Y on prediction of X directly, compare to CFMR2
  res3 <- summary(lm(Y~IV))$coef[2,]
  
  out <- c( n, beta0, h2, mean(res1[,"Estimate"]), sum((res1[,"Std. Error"])^2)/(length(res1[,"Std. Error"]))^2, res2,res3)
  return(out)
}

# start.time <- Sys.time()
# simu_cross_fit_IV_est(n=5000, beta0=1,h2=0.15)
# Sys.time()-start.time
  
beta0 <- c(0, 0.5,-0.5)
h2 <- c(0.2)
n <- c( 2000, 4000, 6000, 8000, 10000)
N         <- rep(n, each=10000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N))
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 6)
save(res, file = "simulation_CFMR_Comparison.RData")

start.time <- Sys.time()
beta0 <- c(0)
h2 <- c(0.2)
n <- c(4000)
N         <- rep(n, each=1000)  
B         <- rep( beta0, length(N)/length(beta0))
H2        <- rep(h2, length(N))
res <- mcmapply(simu_cross_fit_IV_est,n=N, beta0=B,h2=H2, mc.cores = 6)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

