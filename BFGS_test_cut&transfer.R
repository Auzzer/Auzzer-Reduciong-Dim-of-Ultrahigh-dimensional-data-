#导入分位数回归的包
require(quantreg)
library(parallel)
library(snowfall)

OptimCostFunction <- function(eta, gamma1, X, Y, h, tau,eps, coe){
  n <- length(Y)
  lambda <- transfrom_eta_to_lambda(eta, gamma1)
  
  g <- gHatCalculate(X, Y, n, p, lambda, coe)
   
  # 计算方括号前的部分
  part1 <- 1/n * (Y-g)
  # 计算方括号后的部分
  part2 <- abs(rep(tau, n) - pnorm((g-Y)/h))
  # 计算代价函数
  L <- abs(part1 %*% part2)
  return(L)
}
# 给定optim所需梯度函数,公式10第二部分
OptimGradient <- function(eta, gamma1, X, Y, h, tau, eps, coe){
  lambda <- transfrom_eta_to_lambda(eta, gamma1)
  g <- gHatCalculate(X, Y, n, p, lambda, coe)
  deri_g <- deri_g_lambda(X, Y, g, lambda, h, tau, eps)
  deri_cost <- deri_cost_eta(X, Y, gamma1, tau, h, g, eta, deri_g)
  return(deri_cost)
}


# S计算，针对所有的k构造成一个矩阵,k=1,2,...,n(公式8下方）（用以计算g函数）
S_Express <- function(x, p, X, lambda){
  n=nrow(X)
  x_matrix <- rep(x, each=n)
  dim(x_matrix) <- c(n, p)
  # SUM <- lambda^2 %*% t((x_matrix-X)^2)
  # S <- exp(-1/2*SUM)  # S为1*n的矩阵，其中每个分量对应一个S_k
  tmp = (x_matrix-X)
  S = c(rep(0, n))
  for ( j in 1:n){
    S[j] = sum(exp(-(lambda^2 * tmp[j, ]^2)/2))
  }
  return(S)
  
}


# 将eta转化为lambda（公式9至公式10中间的部分）
transfrom_eta_to_lambda <- function(eta, gamma1){
  eta.square <- eta^2
  Sum_of_eta.square <- sum(eta.square)
  lambda <- gamma1 / Sum_of_eta.square * eta.square
  return(lambda)
} 


g_Express<-function(g_last, S, y){
  part1 = sum( (y*dnorm(-(y-g_last))/h)/h+ tau-pnorm(-(y-g_last)/h) *S)
  part2 = sum(dnorm(-(y-g_last)/h)/h*S)
  if(part2<0.01){part2 = 0.01}
  return(part1/part2)
}

get_coe<-function(qw){
  qw = data.frame(qw)
  rq(y~., data = qw)->res
  res$coefficients->coe
  return(coe)
}

gHatCalculate<- function(X, Y, n, p, lambda, coe){
  
  n = length(Y)
  g_final = c(rep(0,n))
  
  temp = coe[-1]*X
  apply(temp, MARGIN = 1, sum)->g0_former
  g0_former = -abs(g0_former+coe[1])
  for (i in 1:n){
    S <- S_Express(X[i,],p, X, lambda )
    g0 = -g0_former[i]
    g1 = g0_former[i]
    g2 = 0
    while(abs(g1-g0)>1e-5){
      
      g2 = g0-(g_Express(g0, S, y=Y[i])-g0)/( (g_Express(g1, S, y=Y[i])-g1) - (g_Express(g0, S, y=Y[i])-g0) )*(g1-g0)
      g1=g2
      g0=g1
    }
    g_final[i]=g2
  }
  return(g_final)
}

# 计算g对lambda的导数(公式11)
deri_g_lambda <- function(X, Y, g, lambda, h, tau, eps){
  # 得到的结果为n*p维矩阵
  n <- length(Y)
  p <- ncol(X)
  deri_g <- rep(1, n*p)
  dim(deri_g) <- c(n, p)
  for(t in 1:p){
    for(i in 1:n){
      x <- X[i,]
      S <- S_Express(x, p, X, lambda)
      S <- as.vector(S)
      # 计算phi与PHI（公式12）
      phi <- dnorm((g-Y)/h)
      # 因为精度无法达到，近似处理
      # phi[phi<min(phi[phi!=0])] <- min(phi[phi!=0])
      PHI <- pnorm((g-Y)/h)
      # 计算deri_phi，phi函数的导数
      deri_phi <- phi * (Y-g)/h
      # 计算A，A是第一个花括号中的公式
      A <- ((phi/h)%*%S) 
      # 计算分母
      # 计算推导的公式中分母的第一项
      part1_1 <- A^2 - A * sum((Y/h*deri_phi-phi)* S/h)  #若取h过小，会导致减号为Inf*NaN
      # 计算推导的公式中分母的第一项
      part1_2 <- ((Y*(phi/h)+tau-PHI) %*% S) * ((deri_phi/h) %*% (S/h))
      # 分母
      part1 <- part1_1 + part1_2
      
      # 计算分子
      # 计算分子的第一项
      part2_1 <-  A * sum((Y/h*phi+tau-PHI) * S * (-lambda[t]) * (rep(x[t], n)-X[,t])^2)
      # 计算分子的第二项
      part2_2 <- ((Y*(phi/h)+tau-PHI) %*% S) * (as.vector(phi*((-lambda[t])*(rep(x[t], n)-X[,t])^2)) %*% (S/h))
      # 分子
      part2 <- part2_1 - part2_2
      
      # 求根
      root <- part2/part1
      
      # 赋值
      deri_g[i, t] <- root
    }
  }
  return(deri_g)
}

# 计算costfunction对eta的导数（公式10第二部分）
deri_cost_eta <- function(X, Y, gamma1, tau, h, g, eta, deri_g){
  n <- length(Y)
  p <- ncol(X)
  # 计算Γ值
  gamma2 <- c(eta) %*% c(eta)
  # 计算求和号前的系数，记为part1
  part1 <- c(2*gamma1*eta) / c(n*gamma2^2)
  # 计算内部花括号中的内容，记作part2_1
  part2_1 <- tau - pnorm((g-Y)/h) + (Y-g)/h*dnorm((g-Y)/h)
  # 根据part2_1计算外部花括号内容，记作part2 (1*p维)
  part2 <- part2_1 %*% deri_g
  # 计算小括号中的内容，记作part3
  part3_1 <- c(gamma2) * diag(p)
  part3 <- c(eta)^2 - part3_1
  
  # 计算求导结果
  deri_cost <- part1*part2 %*% part3
  return(deri_cost)
}


thresold<-function(X){
  if(class(X)== "matrix"){
    for(i in 1:nrow(X)){
      for(j in 1:length(X[i,])){
        if(X[i,j]<1e-3){X[i,j] = 0}else{X[i,j] = X[i,j]}
      }
    }
  }else{for(i in 1:length(X)){
    if(X[i]<1e-3){X[i]=0}else{X[i]=X[i]}
  }}
  return(X)
}



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
  # for(i in 1:3){R = sum(var(X[,i]))}
  # According to :
  # set σ_epsilon^2 so that the theoreticalmodel
  #  R2 = Var{g(X)}/[Var{g(X)} + σ_epsilon^2] is 0.75 (a signalto-noise ratio of 3)
  if (model_select == 1){
    tmp = c(length(X[,2]))
    for (j in 1:length(X[,2])){
      if (X[j,2] > 0.55){
        tmp[j]=0
      }else{
        tmp[j] = X[j,2]
      }
    }
    y = 3.5*X[,1]+1.5*pi*tmp + 2*X[, 3] # + 0.5*ModelError # Model 1
    
    
  }else{
    y = sin(2*pi*(X[,1]+X[,2])/(1+X[,3])) #  + 0.5*ModelError # Model 2
  }
  R = var(y)
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



Result = list()
data = list()
Lambda_result = matrix(0, 20, p)
Loss_result = c(rep(0, 20))

for ( iter in 1:20){
  # Compose and Split Dataset 
  XY = NumericalGeneration(n, p=10, k, model_select=2)
  data[[iter]]=XY
  # XY <- zeroone(XY)
  X <- XY[,1:p]
  Y <- XY[,p+1]
  idx = sample(c(1:n), size = n)
  spilt = 0.8*n
  idx_train = idx[1:spilt]
  idx_test = idx[(spilt+1):n]
  X_train = X[idx_train, ]
  X_test = X[idx_test, ]
  Y_train = Y[idx_train]
  Y_test = Y[idx_test]
  Y = Y_train; X = X_train
  
  # iter 
  eta = c(rep(0.5, p))
  for ( i in 1:p){eta[i] = abs(cor(Y,X[,i]))}
  
  H = list()
  Grad = matrix(0, (maxit+1), p)
  # eta0=eta
  Eta = matrix(0, (maxit+1), p)
  S_k = matrix(0, maxit, p)
  coe = get_coe(XY)
  y_k = matrix(0, maxit, p)
  for ( i in 1:maxit){
    if(i==1){
      # when i = 1, define Hessian = I(p)
      Eta[i, ]= eta
      H[[i]] = diag(p)
      
      
      Grad[i, ] = OptimGradient(Eta[i, ], gamma1, X_train ,Y_train, tau = tau,eps=eps, h=h, coe)#gradient
      Eta[(i+1), ] = Eta[i, ]-t(solve(H[[i]])%*%Grad[i,])
      Grad[(i+1), ] =  OptimGradient(Eta[(i+1), ], gamma1, X_train ,Y_train, tau = tau,eps=eps, h=h, coe)
      
      S_k[i, ]  = Eta[(i+1), ]-Eta[i,]
      y_k[i, ] = Grad[(i+1), ]-Grad[i, ]
      tmp1 = as.numeric(1/(S_k[i, ]%*%H[[i]]%*%S_k[i, ]))
      tmp2 = as.numeric((y_k[i, ]%*%S_k[i, ]))
      H[[(i+1)]] = H[[i]] - (H[[i]]%*%S_k[i, ]%*%t(S_k[i, ])%*%H[[i]])/tmp1 + (y_k[i, ]%*%t(y_k[i, ]))/tmp2
    }else{
      
      # when i>1 
      Grad[i, ] = OptimGradient(Eta[i, ], gamma1, X_train ,Y_train, tau = tau,eps=eps, h=h, coe)#gradient
      Eta[(i+1), ] = Eta[i,]-t(solve(H[[i]])%*%Grad[i,])
      Grad[(i+1), ] =  OptimGradient(Eta[(i+1), ], gamma1, X_train ,Y_train, tau = tau,eps=eps, h=h, coe)
      diff = OptimCostFunction(Eta[(i+1), ], gamma1, X_train, Y_train, h, tau, eps, coe)-OptimCostFunction(Eta[i,], gamma1, X, Y, h, tau, eps, coe)
      if(abs(diff)>1e-32){
        # solve Hessian Matrix
        S_k[i, ]  = Eta[(i+1), ]-Eta[i,]
        y_k[i, ] = Grad[(i+1), ]-Grad[i, ]
        tmp1 = as.numeric(1/(S_k[i, ]%*%H[[i]]%*%S_k[i, ]))
        tmp2 = as.numeric((y_k[i, ]%*%S_k[i, ]))
        H[[(i+1)]] = H[[i]] - (H[[i]]%*%S_k[i, ]%*%t(S_k[i, ])%*%H[[i]])/tmp1 + (y_k[i, ]%*%t(y_k[i, ]))/tmp2
        
      }else{break}
      
    }
    # cat(i, "Completed. ", "\n")  
  }
  
  # calculate final iteration number: itNum
  itNum = maxit
  for(item in 1:maxit){
    if(sum(Eta[item, ]== c(rep(0, p))) >0){
      itNum = item-1
      break
    }
  }
  
thresold<-function(x){
  idx = order(Lambda[119,],decreasing = TRUE)[4:10]
  for(item in idx){
    x[idx]=0
  }
  return(x)
}
  # transform each eta to lambda
  lambda = matrix(0, itNum, p)
  loss = matrix(0, itNum)
  for (j in 1:itNum){
    lambda[j,] = transfrom_eta_to_lambda(Eta[j,], gamma1) 
    loss[j]= OptimCostFunction(lambda[j,], gamma1, X_test, Y_test, h, tau,eps, coe)/n
    
  }
  for(j in 2:itNum){
    if(abs(loss[j]-loss[j-1])<1e-10){
      itNum_final = j
      break
    }
  }
  # return(lambda[itNum_final,], loss[itNum_final]) 
  Lambda_result[iter, ] = lambda[itNum_final,]
  Loss_result[iter] = loss[itNum_final]
  Result[[iter]] = c(loss,lambda)
  cat("Completed: ", iter, "\n")
}