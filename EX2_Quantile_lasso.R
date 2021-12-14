source("./model/Quantile.R")


# require Package: rqPen
# install.packages("rqPen")
zeroone <- function(x) {(x-min(x)) / (max(x) - min(x))}
library("rqPen")
iter = 200 # 重复试验次数 
nlist =c(50, 100, 200, 400) # 样本量
iter = 150 
p = 10

k=0
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
  
  return(cbind(X, y))
  
  
  
}
main_function_Quantile <- function(i){
  # library("rqPen")
  XY <- NumericalGeneration(n,p, k, model_select)
  # XY = zeroone(XY)
  X <- XY[,1:p]
  Y <- XY[,p+1]
  out = cv.rq.pen(X, Y, tau = .5, penalty=penalty, criteria="BIC")
  loss = out$cv[which.min(out$cv[,2]),1]
  lambda = (out$models[[which.min(out$cv[,2])]])$coefficients[-1]
  return(c(loss, lambda))
}



#######################################################################
############################ k = 1 ################################
##########################################################################################
k=1
model_select = 2; penalty = "LASSO"
for (n in nlist){
  assign(paste0("k1Lambda", model_select, penalty, n), matrix(0, iter, p))
}
for (n in nlist){
  assign(paste0("k1Loss", model_select, penalty, n), c(rep(0, iter)))
}
##########
n=50
# main_function_Quantile(n, p, penalty, model_select)
library(parallel)
library(snowfall)
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(rqPen)
sfExport("n","p","penalty", "model_select","k") 
sfExport("NumericalGeneration")
results3_50 <- sfLapply(1:iter, main_function_Quantile)
sfStop()
for (i in 1:iter){
  k1Lambda2LASSO50[i, ] = results3_50[[i]][-1]
  k1Loss2LASSO50[i] = results3_50[[i]][1]
}

##########
n=100
# main_function_Quantile(n, p, penalty, model_select)
library(parallel)
library(snowfall)
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(rqPen)
sfExport("n","p","penalty", "model_select","k") 
sfExport("NumericalGeneration")
results3_100 <- sfLapply(1:iter, main_function_Quantile)
sfStop()
for (i in 1:iter){
  k1Lambda2LASSO100[i, ] = results3_100[[i]][-1]
  k1Loss2LASSO100[i] = results3_100[[i]][1]
}

##########
n=200
# main_function_Quantile(n, p, penalty, model_select)
library(parallel)
library(snowfall)
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(rqPen)
sfExport("n","p","penalty", "model_select","k") 
sfExport("NumericalGeneration")
results3_200 <- sfLapply(1:iter, main_function_Quantile)
sfStop()
for (i in 1:iter){
  k1Lambda2LASSO200[i, ] = results3_200[[i]][-1]
  k1Loss2LASSO200[i] = results3_200[[i]][1]
}

##########
n=400
# main_function_Quantile(n, p, penalty, model_select)
library(parallel)
library(snowfall)
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(rqPen)
sfExport("n","p","penalty", "model_select","k") 
sfExport("NumericalGeneration")
results3_400 <- sfLapply(1:iter, main_function_Quantile)
sfStop()
for (i in 1:iter){
  k1Lambda2LASSO400[i, ] = results3_400[[i]][-1]
  k1Loss2LASSO400[i] = results3_400[[i]][1]
}

email <-compose_email()
email %>%
  smtp_send(
    from = "haozhe_pang@nuist.edu.cn", # 修改发件人
    to = "haozhe_pang@nuist.edu.cn", # 收件人
    subject = "k1 model_select = 2; penalty = LASSO is completed",
    credentials = creds_file(file = "haozhe_nuist")
  )




##########################################################################################
model_select = 2; penalty = "SCAD"
for (n in nlist){
  assign(paste0("k1Lambda", model_select, penalty, n), matrix(0, iter, p))
}
for (n in nlist){
  assign(paste0("k1Loss", model_select, penalty, n), c(rep(0, iter)))
}
##########
n=50
# main_function_Quantile(n, p, penalty, model_select)
library(parallel)
library(snowfall)
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(rqPen)
sfExport("n","p","penalty", "model_select","k") 
sfExport("NumericalGeneration")
results4_50 <- sfLapply(1:iter, main_function_Quantile)
sfStop()
for (i in 1:iter){
  k1Lambda2SCAD50[i, ] = results4_50[[i]][-1]
  k1Loss2SCAD50[i] = results4_50[[i]][1]
}

##########
n=100
# main_function_Quantile(n, p, penalty, model_select)
library(parallel)
library(snowfall)
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(rqPen)
sfExport("n","p","penalty", "model_select","k") 
sfExport("NumericalGeneration")
results4_100 <- sfLapply(1:iter, main_function_Quantile)
sfStop()
for (i in 1:iter){
  k1Lambda2SCAD100[i, ] = results4_100[[i]][-1]
  k1Loss2SCAD100[i] = results4_100[[i]][1]
}

##########
n=200
# main_function_Quantile(n, p, penalty, model_select)
library(parallel)
library(snowfall)
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(rqPen)
sfExport("n","p","penalty", "model_select","k") 
sfExport("NumericalGeneration")
results4_200 <- sfLapply(1:iter, main_function_Quantile)
sfStop()
for (i in 1:iter){
  k1Lambda2SCAD200[i, ] = results4_200[[i]][-1]
  k1Loss2SCAD200[i] = results4_200[[i]][1]
}

##########
n=400
# main_function_Quantile(n, p, penalty, model_select)
library(parallel)
library(snowfall)
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(rqPen)
sfExport("n","p","penalty", "model_select","k") 
sfExport("NumericalGeneration")
results4_400 <- sfLapply(1:iter, main_function_Quantile)
sfStop()
for (i in 1:iter){
  k1Lambda2SCAD400[i, ] = results4_400[[i]][-1]
  k1Loss2SCAD400[i] = results4_400[[i]][1]
}








email <-compose_email()
email %>%
  smtp_send(
    from = "haozhe_pang@nuist.edu.cn", # 修改发件人
    to = "haozhe_pang@nuist.edu.cn", # 收件人
    subject = "k1 model_select = 2; penalty = SCAD is completed",
    credentials = creds_file(file = "haozhe_nuist")
  )





thresold<-function(x){
  for(i in 1:length(x)){
    if(x[i]!=0){x[i]=1}
  }
  return(x)
}
k1Lambda2LASSO100_thresold = t(apply(k1Lambda2LASSO100,MARGIN = 1, FUN=thresold))
k1Lambda2LASSO200_thresold = t(apply(k1Lambda2LASSO200,MARGIN = 1, FUN=thresold))
k1Lambda2LASSO400_thresold = t(apply(k1Lambda2LASSO400,MARGIN = 1, FUN=thresold))
k1Lambda2LASSO50_thresold = t(apply(k1Lambda2LASSO50,MARGIN = 1, FUN=thresold))
########################## M o d e l S e l e c t = 1 ###############
tmp = apply(k1Lambda2LASSO100_thresold, MARGIN = 2, sum)
count=0
for(i in 1:nrow(k1Lambda2LASSO100_thresold)){
  if(k1Lambda2LASSO100_thresold[i,1]==1&k1Lambda2LASSO100_thresold[i,2]==1&k1Lambda2LASSO100_thresold[i,3]==1){count = count+1}
}
PI_k1Lambda2LASSO100 = c(tmp[1]/150, tmp[2]/150, tmp[3]/150, count/150)
###################################################################
tmp = apply(k1Lambda2LASSO200_thresold, MARGIN = 2, sum)
count=0
for(i in 1:nrow(k1Lambda2LASSO200_thresold)){
  if(k1Lambda2LASSO200_thresold[i,1]==1&k1Lambda2LASSO200_thresold[i,2]==1&k1Lambda2LASSO200_thresold[i,3]==1){count = count+1}
}
PI_k1Lambda2LASSO200 = c(tmp[1]/150, tmp[2]/150, tmp[3]/150, count/150)
###################################################################
tmp = apply(k1Lambda2LASSO400_thresold, MARGIN = 2, sum)
count=0
for(i in 1:nrow(k1Lambda2LASSO400_thresold)){
  if(k1Lambda2LASSO400_thresold[i,1]==1&k1Lambda2LASSO400_thresold[i,2]==1&k1Lambda2LASSO400_thresold[i,3]==1){count = count+1}
}
PI_k1Lambda2LASSO400 = c(tmp[1]/150, tmp[2]/150, tmp[3]/150, count/150)
###################################################################
tmp = apply(k1Lambda2LASSO50_thresold, MARGIN = 2, sum)
count=0
for(i in 1:nrow(k1Lambda2LASSO50_thresold)){
  if(k1Lambda2LASSO50_thresold[i,1]==1&k1Lambda2LASSO50_thresold[i,2]==1&k1Lambda2LASSO50_thresold[i,3]==1){count = count+1}
}
PI_k1Lambda2LASSO50 = c(tmp[1]/150, tmp[2]/150, tmp[3]/150, count/150)
###################################################################

PI_lassoQ_ex2 = rbind(PI_k1Lambda2LASSO100, PI_k1Lambda2LASSO200, PI_k1Lambda2LASSO400, PI_k1Lambda2LASSO50)
PI_lassoQ_ex2 = data.frame(PI_lassoQ_ex2)
names(PI_lassoQ_ex2) = c("P1", "P2", "P3", "P-all")



#####################C P 1 , C P 2
CP<-function(x){
  CP1 = c(rep(0,100))
  CP2 = c(rep(0,100))
  for(i in 1:nrow(x)){
    tmp = x[i,]
    CP2[i] = sum(tmp[4:10])/7
    CP1[i] = (3-sum(tmp[1:3]))/3
  }
  return(c(sum(CP1)/100, sum(CP2)/100))
}

CP_k1Lambda2LASSO50 = CP(k1Lambda2LASSO50_thresold)
CP_k1Lambda2LASSO100 = CP(k1Lambda2LASSO100_thresold)
CP_k1Lambda2LASSO200 = CP(k1Lambda2LASSO200_thresold)
CP_k1Lambda2LASSO400 = CP(k1Lambda2LASSO400_thresold)

CP_lassoQ_ex2 = rbind(CP_k1Lambda2LASSO100, CP_k1Lambda2LASSO200,CP_k1Lambda2LASSO400,CP_k1Lambda2LASSO50)

CP_lassoQ_ex2 = data.frame(CP_lassoQ_ex2)
names(CP_lassoQ_ex2) = c("CP1", "CP2")


library(openxlsx)
write.xlsx(CP_lassoQ_ex2, file = "CP_lassoQ_ex2.xlsx", rowNames=TRUE, colNames = TRUE)
write.xlsx(PI_lassoQ_ex2, file = "PI_lassoQ_ex2.xlsx", rowNames=TRUE, colNames = TRUE)



