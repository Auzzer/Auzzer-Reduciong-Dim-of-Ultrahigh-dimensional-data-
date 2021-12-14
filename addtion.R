source("./BFGS_short.R")
gamma1 = 5
h=0.2
p=10

eps = 1e-16
maxit = 66
n1=50;n2=100;n3=200;n4=400
times=100
k=0
###############################################################
tau = 0.5
n=n1
model_select = 2
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(quantreg)
sfExport("eps", "gamma1", "h", "model_select", "n", "p", "tau","maxit","k")
sfExport("deri_cost_eta", "deri_g_lambda", 
         "g_Express", "get_coe", "getRes","gHatCalculate",
         "MEM_BFGS", "NumericalGeneration", 
         "OptimCostFunction", "OptimGradient", "S_Express",
         "thresold","transfrom_eta_to_lambda", "zeroone")

res_n1_ms2_tau5<- sfLapply(1:times, MEM_BFGS)
sfStop()
Lambda_n1_ms2_tau5 = matrix(0, times, p)
Loss_n1_ms2_tau5 = c(rep(0, times))
for (i in 1:times){
  Lambda_n1_ms2_tau5[i, ] = res_n1_ms2_tau5[[i]][-1]
  Loss_n1_ms2_tau5[i] = res_n1_ms2_tau5[[i]][1]
}

###############################################################
tau = 0.5
n=n4
model_select = 2
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(quantreg)
sfExport("eps", "gamma1", "h", "model_select", "n", "p", "tau","maxit","k")
sfExport("deri_cost_eta", "deri_g_lambda", 
         "g_Express", "get_coe", "getRes","gHatCalculate",
         "MEM_BFGS", "NumericalGeneration", 
         "OptimCostFunction", "OptimGradient", "S_Express",
         "thresold","transfrom_eta_to_lambda", "zeroone")

res_n4_ms2_tau5<- sfLapply(1:times, MEM_BFGS)
sfStop()
Lambda_n4_ms2_tau5 = matrix(0, times, p)
Loss_n4_ms2_tau5 = c(rep(0, times))
for (i in 1:times){
  Lambda_n4_ms2_tau5[i, ] = res_n4_ms2_tau5[[i]][-1]
  Loss_n4_ms2_tau5[i] = res_n4_ms2_tau5[[i]][1]
}

###############################################################
tau = 0.25
n=n1
model_select = 2
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(quantreg)
sfExport("eps", "gamma1", "h", "model_select", "n", "p", "tau","maxit","k")
sfExport("deri_cost_eta", "deri_g_lambda", 
         "g_Express", "get_coe", "getRes","gHatCalculate",
         "MEM_BFGS", "NumericalGeneration", 
         "OptimCostFunction", "OptimGradient", "S_Express",
         "thresold","transfrom_eta_to_lambda", "zeroone")

res_n1_ms2_tau25<- sfLapply(1:times, MEM_BFGS)
sfStop()
Lambda_n1_ms2_tau25 = matrix(0, times, p)
Loss_n1_ms2_tau25 = c(rep(0, times))
for (i in 1:times){
  Lambda_n1_ms2_tau25[i, ] = res_n1_ms2_tau25[[i]][-1]
  Loss_n1_ms2_tau25[i] = res_n1_ms2_tau25[[i]][1]
}

###############################################################
tau = 0.75
n=n1
model_select = 2
sfInit(parallel = TRUE, cpus=detectCores()-1)
sfLibrary(quantreg)
sfExport("eps", "gamma1", "h", "model_select", "n", "p", "tau","maxit","k")
sfExport("deri_cost_eta", "deri_g_lambda", 
         "g_Express", "get_coe", "getRes","gHatCalculate",
         "MEM_BFGS", "NumericalGeneration", 
         "OptimCostFunction", "OptimGradient", "S_Express",
         "thresold","transfrom_eta_to_lambda", "zeroone")

res_n1_ms2_tau75<- sfLapply(1:times, MEM_BFGS)
sfStop()
Lambda_n1_ms2_tau75 = matrix(0, times, p)
Loss_n1_ms2_tau75 = c(rep(0, times))
for (i in 1:times){
  Lambda_n1_ms2_tau75[i, ] = res_n1_ms2_tau75[[i]][-1]
  Loss_n1_ms2_tau75[i] = res_n1_ms2_tau75[[i]][1]
}