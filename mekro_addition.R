# mekro p=50 ,real data 实验
##数据的来源包
data("Prostate")
p=ncol(prostate)-1

iter=100
mekro_lambda_real = matrix(0, iter,p)
mekro_loss_real = c(rep(0, iter))
tau.grid = c(seq(0.5, 5, 0.5))
source("./model/mekro/mekro.R")
for(i in 1:iter){ 
  ## Split dataset
  # return: X;y
  gam = c(rep(1, p))
  
  set.seed(i+6666)
  prostate <- subset( prostate, train=TRUE )[,1:9] # X 
  test_sub = sample(nrow(prostate),0.1*nrow(prostate))
  train = prostate[-test_sub,]
  X <- as.matrix(train[,1:8])
  y <- train[,9]
  Xy = cbind(X,y)
  #lasso
  out = mekro.path(X, y, tau.grid = tau.grid, newload = FALSE)
  which.min(out$gcv.mat[,1]) -> idx
  tau = tau.grid[idx]
  
  n = length(y)
  
  
  out = mekro.est(gam = gam,X, y, tau=tau, newload = FALSE)
  mekro_lambda_real[i, ] = out$lam.final
  mekro_loss_real[i] = out$sse / n
  
}

thresold<-function(x){
  for(i in 1:length(x)){
    ifelse(x[i]>0.0005, x[i]<-1, x[i]<- 0)
  }
  return(x)
}
mekro_lambda_real_thresold = apply(mekro_lambda_real, MARGIN = 2, thresold)
apply(mekro_lambda_real_thresold, MARGIN = 2, sum)/100
mean(mekro_loss_real)
NumericalGeneration <- function(n, p, k=0, model_select){
  # set.seed(i+6666)
  # ModelError = rnorm(n, 0, 0.25)
  X <- matrix(0,nrow=n,ncol=p)
  U_star <- c(runif(min=0, max=1, n=n))
  #k <- 0 # 控制协变量的相关性，K=0时相互独立，K=1时，相关性为0.5
  for(i in 1:p){
    X[,i] <- (runif(min=0, max=1, n=n)+k*U_star) / (1+k)
  }
  # Calculate var(epsilon), i.e. measurement error
  # According to :
  # set σ_epsilon^2 so that the theoreticalmodel
  # R2 = Var{g(X)}/[Var{g(X)} + σ_epsilon^2] is 0.75 (a signalto-noise ratio of 3)
  ##  calculate Var(g(x))
  if (model_select == 1){
    tmp = c(length(X[,2]))
    for (j in 1:length(X[,2])){
      if (X[j,2] > 0.55){
        tmp[j]=0
      }else{
        tmp[j] = X[j,2]
      }
    }
    y = 3.5*X[,1]+1.5*pi*tmp + 2*X[, 3]  # Model 1
    
    
  }else{
    y = sin(2*pi*(X[,1]+X[,2])/(1+X[,3]))  # Model 2
  }
  R = var(y)# var(g(x))
  ModelError = rnorm(n, 0, sqrt(R/3))
  
  
  if (model_select == 1){
    tmp = c(length(X[,2]))
    for (j in 1:length(X[,2])){
      if (X[j,2] > 0.55){
        tmp[j]=0
      }else{
        tmp[j] = X[j,2]
      }
    }
    y = 3.5*X[,1]+1.5*pi*tmp + 2*X[, 3] + 0.5*ModelError # Model 1
    
    
  }else{
    y = sin(2*pi*(X[,1]+X[,2])/(1+X[,3]))  + 0.5*ModelError # Model 2
  }
  # XY = zeroone(cbind(X,y))
  return(cbind(X, y))
  
  
  
}

#################################################################
#############p = 50
k=0
##################################################################################
# Part 1: get tau
# description: get tau respond to n like taui_j to ni and model_select = j
# out put: tau1_1, tau1_2, tau2_1, tau2_2, tau3_1, tau3_2, tau4_1, tau4_2
###################################################################################
# tau1_1
iter=100
p=50
k=0
n=100
XY <- zeroone(NumericalGeneration(n,p, k, model_select=1))

X <- XY[,1:p]
Y <- XY[,p+1]

out = mekro.path(X, Y, tau.grid = tau.grid, newload = FALSE)
which.min(out$gcv.mat[,1]) -> idx
tau1 = tau.grid[idx]

# tau3
n=100
XY <- zeroone(NumericalGeneration(n,p, k, model_select=2))

X <- XY[,1:p]
Y <- XY[,p+1]

out = mekro.path(X, Y, tau.grid = tau.grid, newload = FALSE)
which.min(out$gcv.mat[,1]) -> idx
tau2 = tau.grid[idx]

gam = c(rep(1/tau1,p))
resultMekroN2ms1<- matrix(0, iter, p)
LossMekroN2ms1<-c(rep(0,iter))
for(i in 1:iter){
  
  XY <- zeroone(NumericalGeneration(n,p, k, model_select=1))
  X <- XY[,1:p]
  Y <- XY[,p+1]
  out = mekro.est(gam = gam,X, Y, tau=tau1, newload = FALSE)
  resultMekroN2ms1[i, ] = out$lam.final
  LossMekroN2ms1[i] = out$sse / n
}



gam = c(rep(1/tau2,p))
resultMekroN2ms2<- matrix(0, iter, p)
LossMekroN2ms2<-c(rep(0,iter))
for(i in 1:iter){
  
  XY <- zeroone(NumericalGeneration(n,p, k, model_select=2))
  X <- XY[,1:p]
  Y <- XY[,p+1]
  out = mekro.est(gam = gam,X, Y, tau=tau1, newload = FALSE)
  resultMekroN2ms2[i, ] = out$lam.final
  LossMekroN2ms2[i] = out$sse / n
}

thresold<-function(x){
  for(i in 1:length(x)){
    ifelse(x[i]>0.0001, x[i]<-1, x[i]<- 0)
  }
  return(x)
}
resultMekroN2ms1_thresold<-thresold(resultMekroN2ms1)
resultMekroN2ms2_thresold<-thresold(resultMekroN2ms2)


################################################################
####################### Pi C a c u l a t e#######################
################################################################
PI <- function(X){
  tmp = apply(X, MARGIN = 2, sum)
  count = 0
  for(i in 1:nrow(X)){
    if(X[i,1]==1&X[i,2]==1&X[i,3]==1){
      count = count+1}
  }
  return(c(tmp[1]/100, tmp[2]/100, tmp[3]/100, count/100 ))
}

PI_MekroN2ms1 = PI(resultMekroN2ms1_thresold)
PI_MekroN2ms2 = PI(resultMekroN2ms2_thresold)

CP<-function(X){
  CP2 = c(rep(0, 100))
  CP1 = c(rep(0,100))
  for(i in 1:nrow(X)){
    tmp = X[i, ]
    CP2[i] = CP2[i] = sum(tmp[4:p])/7
    CP1[i] = (3-sum(tmp[1:3]))/3
  }
  return(c(sum(CP1)/100,sum(CP2)/100))
}

CP_MekroN2ms1 = CP(resultMekroN2ms1_thresold)

CP_MekroN2ms2 = CP(resultMekroN2ms2_thresold)





