##################################################################################

# Includes Power simulation of Path-based test statistics.

##################################################################################


#######
##### loading required packages
source("~/hdi_path/bin/pathwise_require.R")
source("~/hdi_path/bin/pathwise_ts.R")
source("~/hdi_path/bin/pathwise_resample.R")
source("~/hdi_path/bin/pathwise_simu_setting.R")
source("~/hdi_path/bin/pathwise_net_resample.R")
source("~/hdi_path/bin/pathwise_ts.R")

###### 




Path.Resample.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, iter = 500, B = 500, setting = 'dep', which.covariate = 1, betaNull = 1, multiTest = FALSE, ...){
#Return the power fo Path-based test under different settings    
# Args:
#	setting: different settings, check 'pathwise_simu_setting.R for details'	
# 	rho: related to dependent design setting
#	n,p,beta : sample size, features, coefficients
# 	iter : # of iterations 
#
# Return:
#	Simulated power
#

	#TS = matrix(NA, iter, len(which.covariate))
	#b_size_la <-  matrix(0,iter,4)
  
  	path.power = array(0,dim = c(iter,len(which.covariate),4))

	for(s in 1:iter){
    
  		data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)

		X_sp = data$X
		Y_sp = data$Y
		
		results = Path.Resample(X = X_sp, Y = Y_sp, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest, B = B, beta.true = beta, ...)

		######## keep track of MEMORY
		print("After Bootstrap:")
		print(mem_used())  
		########
		
		path.power[s,,] = results$rej

		if(s %% 10 == 0){  cat("Now Computing:",s,"\n") }
        
   
	}    

        
	path.power = apply(path.power,c(2,3),mean)
	
	return(path.power = path.power)
    
   
}    
    







Net.Resample.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, iter = 500, B = 500, setting = 'dep', which.covariate = 1, betaNull = 1, multiTest = FALSE, ...){
#Return the power fo Path-based test under different settings    
# Args:
#	setting: different settings, check 'pathwise_simu_setting.R for details'	
# 	rho: related to dependent design setting
#	n,p,beta : sample size, features, coefficients
# 	iter : # of iterations 
#
# Return:
#	Simulated power
#

	#TS = matrix(NA, iter, len(which.covariate))
	#b_size_la <-  matrix(0,iter,4)
  
  	path.power = array(0,dim = c(iter,len(which.covariate),4))

  	results = list()
  	results$L2.squared = path.power
  	results$L2 = path.power
  	results$L1 = path.power
  	results$L_inf = path.power

	for(s in 1:iter){
    
  		data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)

		X_sp = data$X
		Y_sp = data$Y
		
		res = Net.Resample(X = X_sp, Y = Y_sp, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest, B = B, beta.true = beta, ...)

		######## keep track of MEMORY
		print("After Bootstrap:")
		print(mem_used())  
		########
		
		results$L2.squared[s,,] = res$L2.squared$rej
		results$L2[s,,] = res$L2$rej
		results$L1[s,,] = res$L1$rej
		results$L_inf[s,,] = res$L_inf$rej
		
		#path.power[s,,] = results$rej

		if(s %% 10 == 0){  cat("Now Computing:",s,"\n") }
        
   
	}    

        
	results$L2.squared = apply(results$L2.squared, c(2,3),mean)
	results$L2 = apply(results$L2, c(2,3),mean)
	results$L1 = apply(results$L1, c(2,3),mean)
	results$L_inf = apply(results$L_inf, c(2,3),mean)
	


	return(results)
    
   
}    
    

