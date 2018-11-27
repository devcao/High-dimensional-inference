##################################################################################

# # Includes Power simulation of Projection method, etc.

##################################################################################




#######
##### loading required packages
source("~/hdi_path/bin/pathwise_require.R")
source("~/hdi_path/bin/pathwise_simu_setting.R")
library(scalreg)
source(file.path("~/SI/Sim-CI.R"))
source(file.path("~/SI/ST.R"))

######


#ST <- function(X.f, Y.f, sub.size, test.set, M=500, alpha=c(0.2,0.1,0.05,0.01))

################ de-sparsified lasso #########################
ST.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate){
#Return the power of de-sparsified under different settings    
# Args:
# setting: different settings, check 'pathwise_simu_setting.R for details'  
#   rho: related to dependent design setting
# n,p,beta : sample size, features, coefficients
#   iter : # of iterations 
#
# Return:
# A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
#  

  ST.power.nst = matrix(0,iter,4)
  ST.power.st = matrix(0,iter,4)
  
  #pval = matrix(NA, iter, p)
  
  for(s in 1:iter){

    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)



    X_sp = data$X
    Y_sp = data$Y

    sub.size <- n*0.3
    
    results = ST(X.f = X_sp, Y.f = Y_sp, sub.size = sub.size, test.set = which.covariate)
    print(results)
    ST.power.nst[s,] <- results$nst
    ST.power.st[s,] <- results$st
    
    
    if(s %% 10 == 0){  cat("Now computing:", s, "\n")  }
    
  }
 
  simCI.power.nst=apply(ST.power.nst,2,mean)  
  
  simCI.power.st=apply(ST.power.st,2,mean)  
  
  return(list(nst = ST.power.nst, st = ST.power.st))

}  

##############################################################




######
################ de-sparsified lasso #########################
simCI.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate, betaNull){
#Return the power of de-sparsified under different settings    
# Args:
# setting: different settings, check 'pathwise_simu_setting.R for details'  
#   rho: related to dependent design setting
# n,p,beta : sample size, features, coefficients
#   iter : # of iterations 
#
# Return:
# A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
#  

  simCI.power.nst = matrix(0,iter,4)
  simCI.power.st = matrix(0,iter,4)
  
  #pval = matrix(NA, iter, p)
  
  for(s in 1:iter){

    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)



    X_sp = data$X
    Y_sp = data$Y
    
    results = Sim.CI(X = X_sp, Y = Y_sp, set = which.covariate, betaNull = betaNull)

    simCI.power.nst[s,] <- results$nst
    simCI.power.st[s,] <- results$st
    
    
    if(s %% 10 == 0){  cat("Now computing:", s, "\n")  }
    
  }
 
  simCI.power.nst=apply(simCI.power.nst,2,mean)  
  
  simCI.power.st=apply(simCI.power.st,2,mean)  
  
  return(list(nst = simCI.power.nst, st = simCI.power.st))

}  

##############################################################




################ de-sparsified lasso #########################
desparse.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate, betaNull){
#Return the power of de-sparsified under different settings    
# Args:
#	setting: different settings, check 'pathwise_simu_setting.R for details'	
# 	rho: related to dependent design setting
#	n,p,beta : sample size, features, coefficients
# 	iter : # of iterations 
#
# Return:
#	A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
#  
  
  if(len(betaNull) > 1){stop("now only support compute power of 1 coefficients")}

  proj.power = matrix(0,len(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){

    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)



   	X_sp = data$X
    Y_sp = data$Y
    
    fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
      
    pval[s,] = fit.proj$pval
    
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }

  count = 1
  for(j in which.covariate){

    proj.power[count,1] = mean(pval[,j] < 0.2)
    proj.power[count,2] = mean(pval[,j] < 0.1)
    proj.power[count,3] = mean(pval[,j] < 0.05)
    proj.power[count,4] = mean(pval[,j] < 0.01)
  

    count = count + 1
  }
  
  
  return(proj.power)

}  

##############################################################



