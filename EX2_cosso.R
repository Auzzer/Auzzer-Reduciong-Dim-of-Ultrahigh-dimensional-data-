# 结果分别保存至result：lambda；cvm：loss
library(cosso)
n1 = 50; n2=100; n3=200; n4=400
p = 10

# Numerical Generation
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
zeroone <- function(x) (x-min(x)) / (max(x) - min(x))


iter = 100
k=1
###############################################################
###############################################################
model_select=2
# tau = 0.5
k1LambdaCOSSO50 = matrix(0, iter, p)
k1LambdaCOSSO100 = matrix(0, iter, p)
k1LambdaCOSSO200 = matrix(0, iter, p)
k1LambdaCOSSO400 = matrix(0, iter, p)


k1LossCosso50 = c(rep(0, iter))
k1LossCosso100 = c(rep(0, iter))
k1LossCosso200 = c(rep(0, iter))
k1LossCosso400 = c(rep(0, iter))


n=50
n=n1
for (i in 1:iter){
  XY <-NumericalGeneration(n,p, k, model_select)
  XY = zeroone(XY)
  X <- XY[,1:p]
  y <- XY[,p+1]
  
  
  out <- cosso(X, y, tau = 0.5, family = "Quan", cpus = 11)
  tune = tune.cosso(out, plot.it = FALSE)
  cvm = min(tune$cvm)
  result = tune$M[2:11]
  k1LambdaCOSSO50[i,] = result
  k1LossCosso50[i] = cvm 
}


n=100
n=100
for (i in 1:iter){
  XY <-NumericalGeneration(n,p, k, model_select)
  XY = zeroone(XY)
  X <- XY[,1:p]
  y <- XY[,p+1]
  
  
  out <- cosso(X, y, tau = 0.5, family = "Quantile", cpus = 11)
  tune = tune.cosso(out, plot.it = FALSE)
  cvm = min(tune$cvm)
  result = tune$M[2:11]
  k1LambdaCOSSO100[i,] = result
  k1LossCosso100[i] = cvm 
}
n=200
for (i in 1:iter){
  XY <-NumericalGeneration(n,p, k, model_select)
  XY = zeroone(XY)
  X <- XY[,1:p]
  y <- XY[,p+1]
  
  out <- cosso(X, y, tau = 0.5, family = "Quantile", cpus = 11)
  tune = tune.cosso(out, plot.it = FALSE)
  cvm = min(tune$cvm)
  result = tune$M[2:11]
  k1LambdaCOSSO200[i,] = result
  k1LossCosso200[i] = cvm 
}




n=400
for (i in 1:iter){
  XY <-NumericalGeneration(n,p, k, model_select)
  XY = zeroone(XY)
  X <- XY[,1:p]
  y <- XY[,p+1]
  
  out <- cosso(X, y, tau = 0.5, family = "Quantile", cpus = 11)
  tune = tune.cosso(out, plot.it = FALSE)
  cvm = min(tune$cvm)
  result = tune$M[2:11]
  k1LambdaCOSSO400[i,] = result
  k1LossCosso400[i] = cvm 
}
email %>%
  smtp_send(
    from = "haozhe_pang@nuist.edu.cn", # 修改发件人
    to = "haozhe_pang@nuist.edu.cn", # 收件人
    subject = "cosso_model1_tau=0.5完成",
    credentials = creds_file(file = "haozhe_nuist")
  )
