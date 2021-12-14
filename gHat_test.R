
#导入分位数回归的包
require(quantreg)
library(quantreg)


OptimCostFunction <- function(eta, gamma1, X, Y, h, tau,eps, coe){
  m=n
  lambda <- transfrom_eta_to_lambda(eta, gamma1)
  
  g <- gHatCalculate(X, Y, m, p, lambda, coe)
  m <- length(Y) 
  # 计算方括号前的部分
  part1 <- 1/m * (Y-g)
  # 计算方括号后的部分
  part2 <- abs(rep(tau, m) - pnorm((g-Y)/h))
  # 计算代价函数
  L <- abs(part1 %*% part2)
  return(L)
}
# 给定optim所需梯度函数,公式10第二部分
OptimGradient <- function(eta, gamma1, X, Y, h, tau, eps, coe){
  lambda <- transfrom_eta_to_lambda(eta, gamma1)
  g <- gHatCalculate(X, Y, m, p, lambda, coe)
  deri_g <- deri_g_lambda(X, Y, g, lambda, h, tau, eps)
  deri_cost <- deri_cost_eta(X, Y, gamma1, tau, h, g, eta, deri_g)
  return(deri_cost)
}


# S计算，针对所有的k构造成一个矩阵,k=1,2,...,n(公式8下方）（用以计算g函数）
S_Express <- function(x, m, p, X, lambda){
  x_matrix <- rep(x, each=m)
  dim(x_matrix) <- c(m, p)
  # SUM <- lambda^2 %*% t((x_matrix-X)^2)
  # S <- exp(-1/2*SUM)  # S为1*n的矩阵，其中每个分量对应一个S_k
  tmp = (x_matrix-X)
  S = c(rep(0, m))
  for ( j in 1:m){
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

gHatCalculate<- function(X, Y, m, p, lambda, coe){
  g_final = c(rep(0,m))
  
  temp = coe[-1]*X
  apply(temp, MARGIN = 1, sum)->g0_former
  g0_former = -abs(g0_former+coe[1])
  for (i in 1:m){
    S <- S_Express(X[i,], m,p, X, lambda )
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
  m <- length(Y)
  p <- ncol(X)
  deri_g <- rep(1, m*p)
  dim(deri_g) <- c(m, p)
  for(t in 1:p){
    for(i in 1:m){
      x <- X[i,]
      S <- S_Express(x, m, p, X, lambda)
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
      part2_1 <-  A * sum((Y/h*phi+tau-PHI) * S * (-lambda[t]) * (rep(x[t], m)-X[,t])^2)
      # 计算分子的第二项
      part2_2 <- ((Y*(phi/h)+tau-PHI) %*% S) * (as.vector(phi*((-lambda[t])*(rep(x[t], m)-X[,t])^2)) %*% (S/h))
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
  # m <- length(Y)
  # p <- ncol(X)
  # 计算Γ值
  gamma2 <- c(eta) %*% c(eta)
  # 计算求和号前的系数，记为part1
  part1 <- c(2*gamma1*eta) / c(m*gamma2^2)
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






# Numercial Generation
NumercialGeneration <- function(n, p, k=0, model_select){
  # set.seed(i+6666)
  ModelError = rnorm(n, 0, 0.25)
  X <- matrix(0,nrow=n,ncol=p)
  U_star <- c(runif(min=0, max=1, n=n))
  #k <- 0 # 控制协变量的相关性，K=0时相互独立，K=1时，相关性为0.5
  for(i in 1:p){
    X[,i] <- (runif(min=0, max=1, n=n)+k*U_star) / (1+k)
  }
  
  
  if (model_select == 1){
    tmp = c(length(X[,2]))
    for (j in 1:length(X[,2])){
      if (X[j,2] > 0){
        tmp[j]=0
      }else{
        tmp[j] = X[j,2]
      }
    }
    y = 3.5*X[,1]+1.5*pi*tmp + 2*X[, 3] + ModelError # Model 1
    
    
  }else{
    y = sin((X[,1]+X[,2])/(1+X[,3])) + ModelError # Model 2
  }
  return(cbind(X, y))
  
  
  
}

episode<- function(i){
  set.seed(i+6666)
  
  idx = sample(c(1:n) ,n)
  group = list()
  for (i in 1:(n/m)){
    low = m*i-(m-1)
    upper = m*i
    group[[i]] = idx[low:upper]
  }
  
  
  for (i in 1:(n/m)){
    tmp = group[[i]]
    X = q[tmp, ]
    Y = w[tmp]
    
    grad = OptimGradient(eta, gamma1, X,Y, h, eps=eps, tau = tau, coe)
    eta = eta - learning_rate*(1/m)*grad
    
  }
  return(eta)
}



gamma1 = 5
h=0.2
m = 20 # batch size
n=400
p=100
learning_rate = 0.5
tau = 0.5 


eps = 1e-16
lambda_list = list
for (i in 1:200){
  qw = NumercialGeneration(n=n, p=p, k=0, model_select=2)
  q <- qw[,1:p]
  w <- qw[,p+1]
  
  eta = c(rep(0, p))
  for ( i in 1:p){eta[i] = abs(cor(w,q[,i]))}
  
  coe = get_coe(qw)
  
  
  
  times = 400-n
  lambda = matrix(rep(0, times*p), times, p)
  loss = c(rep(0, times))
  for (i in 1:times){
    eta = episode(i)
    lambda[i, ] = transfrom_eta_to_lambda(eta, gamma1)
    loss[i] = OptimCostFunction(eta, gamma1, X=q, Y=w, h, tau, eps, coe)
    if(i%%2==0){
      cat("\n lambda: ",lambda[i, ],"\n Loss:", loss[i] , 
          "Completed: ", i/times*100, "% \n")
    }
  }
  
  
}


