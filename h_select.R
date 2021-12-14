n=100
p=10
XY <- NumercialGeneration(n,p, k, model_select=2)
X <- XY[,1:p]
Y <- XY[,p+1]

# 以下为二分法中用到的一些参数
g0 <- rep(-max(Y)-5, times=n)
g1 <- rep(max(Y)+5, times=n) 

library(ncvreg)
out.lasso = ncvreg(X, Y, family = "gaussian", penalty = "lasso" )
eta = out.lasso$beta[, which.min(out.lasso$lambda)][-1]

h=0.01
optim.out1 <- optim(eta, OptimCostFunction, OptimGradient, method = "L-BFGS-B",control = list(factr = 1e-16,maxit = 66),
                    gamma1=gamma1, X=X, Y=Y, h=h, tau=tau, g0=g0, g1=g1, eps=eps) 




h=0.1
optim.out2 <- optim(eta, OptimCostFunction, OptimGradient, method = "L-BFGS-B",control = list(factr = 1e-16,maxit = 66),
                    gamma1=gamma1, X=X, Y=Y, h=h, tau=tau, g0=g0, g1=g1, eps=eps) 




h=1
optim.out3 <- optim(eta, OptimCostFunction, OptimGradient, method = "L-BFGS-B",control = list(factr = 1e-16,maxit = 66),
                    gamma1=gamma1, X=X, Y=Y, h=h, tau=tau, g0=g0, g1=g1, eps=eps) 




h=10
optim.out4 <- optim(eta, OptimCostFunction, OptimGradient, method = "L-BFGS-B",control = list(factr = 1e-16,maxit = 66),
                    gamma1=gamma1, X=X, Y=Y, h=h, tau=tau, g0=g0, g1=g1, eps=eps) 



optim <- rbind(eta, optim.out1$par,optim.out2$par,optim.out3$par,optim.out4$par)

cbind(optim.out1$value,optim.out2$value,optim.out3$value, optim.out4$value)


# 结果显示 确实有



