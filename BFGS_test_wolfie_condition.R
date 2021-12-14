#导入分位数回归的包
require(quantreg)
library(parallel)
library(snowfall)

OptimCostFunction <- function(lambda, gamma1, X, Y, h, tau,eps, coe){
  # m=n
  # lambda <- transfrom_eta_to_lambda(eta, gamma1)
  
  g <- gHatCalculate(X, Y, n, p, lambda, coe)
  n <- length(Y) 
  # 计算方括号前的部分
  part1 <-  (Y-g)
  # 计算方括号后的部分
  part2 <- rep(tau, n) - pnorm((g-Y)/h)
  # 计算代价函数
  L <- abs(part1 %*% part2)/n
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
S_Express <- function(x, n, p, X, lambda){
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
  g_final = c(rep(0,n))
  
  temp = coe[-1]*X
  apply(temp, MARGIN = 1, sum)->g0_former
  g0_former = -abs(g0_former+coe[1])
  for (i in 1:n){
    S <- S_Express(X[i,], n,p, X, lambda )
    g0 = -g0_former[i]
    g1 = g0_former[i]
    g2 = 0
    while(abs(g1-g0)>1e-10){
      
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
      S <- S_Express(x, n, p, X, lambda)
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
  # p <- ncol(X)
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





NumericalGeneration <- function(n, p, k=0, model_select){
  # set.seed(i+6666)
  # ModelError = rnorm(n, 0, 0.25)
  X <- matrix(0,nrow=n,ncol=p)
  U_star <- c(runif(min=0, max=1, n=n))
  #k <- 0 # 控制协变量的相关性，K=0时相互独立，K=1时，相关性为0.5
  for(i in 1:p){
    X[,i] <- (runif(min=0, max=1, n=n)+k*U_star) / (1+k)
  }
  X = zeroone(X)
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
  return(zeroone(cbind(X, y)))
  
  
  
}



thresold<-function(x,l){
  # l:选择的变量个数
  if(l<p){
    idx = order(x,decreasing = TRUE)[(l+1):10]
    for(item in idx){
      x[idx]=0
    }
    return(x)
  }else if(l==p){
    return(x)
  }else {
    cat("Warning: the number of chosen variable should be smaller than the length of variable")
  }
  
}


getRes<-function(p,Lambda_final, gamma1, X, Y, h, tau,eps, coe, itNum){
  
  res<-list()
  for(l in 3:p){
    loss2 = c(rep(0, itNum))
    BIC_Quantile = c(rep(0, itNum))
    BIC = c(rep(0, itNum))
    AIC = c(rep(0, itNum))
    
    Lambda_thresold = t(apply(Lambda_final, MARGIN=1, thresold, l=l))
    loss2 = apply(Lambda_thresold, MARGIN = 1,  OptimCostFunction,gamma1, X, Y, h, tau, eps, coe)
    for (z in 1:length(loss2)){
      BIC_Quantile[z] = log(loss2[z]*n)+ l*log(n)/(2*sqrt(n))
      BIC[z] = l*log(n)-2*log(loss2[z]*n)
      AIC[z] = 2*l-2*log(loss2[z]*n)
    }
    res[[l]] = cbind(loss2, BIC_Quantile, BIC, AIC)
  }
  return(res)
}

zeroone <- function(x) {(x-min(x)) / (max(x) - min(x))}



MEM_BFGS<-function(times){
  
  XY = NumericalGeneration(n, p=10, k, model_select)
  # XY <- zeroone(XY)
  X <- XY[,1:p]
  Y <- XY[,p+1]
  
  eta = c(rep(0.5, p))
  for ( i in 1:p){eta[i] = abs(cor(Y,X[,i]))}
  for (i in 1:3){if(eta[i]<0.5){eta[i] = 0.5}}
  H = list()
  Grad = matrix(0, (maxit+1), p)
  # eta0=eta
  Eta = matrix(0, (maxit+1), p)
  S_k = matrix(0, maxit, p)
  coe = get_coe(XY)
  y_k = matrix(0, maxit, p)
  Lambda =  matrix(0, (maxit+1), p)
  for ( i in 1:maxit){
    if(i==1){
      # when i = 1, define Hessian = I(p)
      Eta[i, ]= eta
      Lambda[i, ] = transfrom_eta_to_lambda(Eta[i, ], gamma1)
      H[[i]] = diag(p)
      Grad[i, ] = OptimGradient(Eta[i, ], gamma1, X ,Y, tau = tau,eps=eps, h=h, coe)#gradient
      S_k[i,] = -t(solve(H[[i]])%*%Grad[i,])
      
      ## step choose with backstriking linge search
      t = 1
        
      Eta[(i+1), ] = Eta[i, ]+t*S_k[i,]
      
      Lambda[(i+1), ] = transfrom_eta_to_lambda(Eta[(i+1), ], gamma1)
      f1 = OptimCostFunction(Lambda[i, ], gamma1, X, Y, h, tau, eps, coe)
      f2 = OptimCostFunction(Lambda[(i+1), ], gamma1, X, Y, h, tau, eps, coe)
      alpha = 0.2
      beta = 0.5
      differ = f2-f1
      
      while(f2-f1+alpha*t*as.numeric(t(Grad[i,])%*%S_k[i,])>0&t>1e-5){
        t = beta*t
        Eta[(i+1), ] = Eta[i, ]+t*S_k[i,]
        
        Lambda[(i+1), ] = transfrom_eta_to_lambda(Eta[(i+1), ], gamma1)
        f1 = OptimCostFunction(Lambda[i, ], gamma1, X, Y, h, tau, eps, coe)
        f2 = OptimCostFunction(Lambda[(i+1), ], gamma1, X, Y, h, tau, eps, coe)
        differ = f2-f1
      }
        
        
      Grad[(i+1), ] =  OptimGradient(Eta[(i+1), ], gamma1, X ,Y, tau = tau,eps=eps, h=h, coe)
      y_k[i, ] = Grad[(i+1), ]-Grad[i, ]
      tmp1 = as.numeric(1/(S_k[i, ]%*%H[[i]]%*%S_k[i, ]))
      tmp2 = as.numeric((y_k[i, ]%*%S_k[i, ]))
      
      H[[(i+1)]] = H[[i]] - (H[[i]]%*%S_k[i, ]%*%t(S_k[i, ])%*%H[[i]])/tmp1 + (y_k[i, ]%*%t(y_k[i, ]))/tmp2
      
      
    }else if(i>1&differ>1e-10){
      
      # when i>1 
      Grad[i, ] = OptimGradient(Eta[i, ], gamma1, X ,Y, tau = tau,eps=eps, h=h, coe)#gradient
      
      
      S_k[i,] = -t(solve(H[[i]])%*%Grad[i,])
      
      ## step choose with backstriking linge search
      t = 1
      
      Eta[(i+1), ] = Eta[i, ]+t*S_k[i,]
      
      Lambda[(i+1), ] = transfrom_eta_to_lambda(Eta[(i+1), ], gamma1)
      f1 = OptimCostFunction(Lambda[i, ], gamma1, X, Y, h, tau, eps, coe)
      f2 = OptimCostFunction(Lambda[(i+1), ], gamma1, X, Y, h, tau, eps, coe)
      alpha = 0.08
      beta = 0.5
      
      
      while(f2-f1+alpha*t*as.numeric(t(Grad[i,])%*%S_k[i,])>0&t>1e-5){
        t = beta*t
        Eta[(i+1), ] = Eta[i, ]+t*S_k[i,]
        
        Lambda[(i+1), ] = transfrom_eta_to_lambda(Eta[(i+1), ], gamma1)
        f1 = OptimCostFunction(Lambda[i, ], gamma1, X, Y, h, tau, eps, coe)
        f2 = OptimCostFunction(Lambda[(i+1), ], gamma1, X, Y, h, tau, eps, coe)
        
      }
      
      
      Grad[(i+1), ] =  OptimGradient(Eta[(i+1), ], gamma1, X ,Y, tau = tau,eps=eps, h=h, coe)
      y_k[i, ] = Grad[(i+1), ]-Grad[i, ]
      tmp1 = as.numeric(1/(S_k[i, ]%*%H[[i]]%*%S_k[i, ]))
      tmp2 = as.numeric((y_k[i, ]%*%S_k[i, ]))
      H[[(i+1)]] = H[[i]] - (H[[i]]%*%S_k[i, ]%*%t(S_k[i, ])%*%H[[i]])/tmp1 + (y_k[i, ]%*%t(y_k[i, ]))/tmp2
      
    
      
    }
    cat(i, "Completed. ", "\n")  
    cat(Lambda[i,])
  }
  
  
  itNum = maxit #计算实际迭代次数
  for(item in 1:maxit){
    if(sum(Eta[item, ]== c(rep(0, p))) >0){
      itNum = item-1
      break
    }
  }
  Eta  = Eta[1:itNum,]
  Lambda_final = t(apply(Eta, MARGIN=1, transfrom_eta_to_lambda, gamma1))
  
  
  
  # judge部分
  out<-getRes(p,Lambda_final, gamma1, X, Y, h,tau, eps,coe,itNum)
  
  Variable_Select = matrix(0, length(c(3:p)), 5) # 合计五列，分别为：迭代次数，loss，分位数BIC，BIC和AIC
  for(i in 3:p){
    idx = which.min(out[[i]][,2])
    Variable_Select[i, ] = c(idx,out[[i]][idx,])
  }
  VarNum = which.min(Variable_Select[,3]) # 变量个数
  return(c(Variable_Select[VarNum,2],
           thresold(Lambda_final[Variable_Select[VarNum,1],], VarNum)))# 第一列为loss，后面为Lambda结果果

         
}





# MEM_BFGS(1)




