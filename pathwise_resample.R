##################################################################################

# Includes all Pathwise test statistics, including end point variant and variate norm 

##################################################################################



#######
##### loading required packages, functions
#source("~/hello/hdi_path/pathwise_require.R")
#source("~/hello/hdi_path/pathwise_ts.R")
###### 
source("~/hdi_path/bin/pathwise_require.R")

source("~/hdi_path/bin/pathwise_ts.R")



Path.Resample = function(X, Y, which.covariate, betaNull, multiTest, B = 500, parallel = FALSE, exact = TRUE, beta.init = 'adaptive', beta.true = beta, ...){
# Bootstrap the null distribution of Path-based statistic, and return reject or not
#
# Args:
#	X, Y, which.covariate : feed in to Path-based TS function
#	B : # of bootstrap replications
#	parallel : run in parallel or nor
#	exact : use exact TS or approx TS
#	beta.true : for simulation only, the true value of beta would be used as initial estimates.
#	
# Return:  
#	Reject or not under alpha = 0.2,0.1,0.05,0.01
# 	p.values of the test

	
  n = nrow(X)
  p = ncol(X)
  rej = matrix(0,len(which.covariate),4)  # 1st : which cov, 2nd: na, 3rd: which alpha, 0.2,0.1,0.05,0.01
  pval = numeric()

  TS = Path.TS(exact = exact, X = X, Y = Y, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest,...)
  
  



  if(p >= n){ # high dimension we can try...

  	if(beta.init == "adaptive"){
  		bhat = adalasso(X = X, y = Y, k = 10, use.Gram = FALSE,both = TRUE, intercept = FALSE)$coefficients.adalasso
  
    	}else if (beta.init == "de-sparse"){
    		
    		bhat = as.vector(lasso.proj(X, Y, standardize = TRUE, parallel = TRUE, ncores = 40)$bhat)

    		}else if (beta.init == "MC+"){
    			
    			bhat = coef(cv.ncvreg(X = X, y = Y, penalty = "MCP",family = "gaussian", nfold= 10))[-1]
          
    			}else if (beta.init == "SCAD"){
    				bhat = coef(cv.ncvreg(X = X, y = Y, penalty = "SCAD",family = "gaussian", nfold= 10))[-1]

    			}else if (beta.init == "Truth"){
    				bhat = beta_true
    			} 

    
    residual = Y - X%*%bhat
    
  	}else{ # low dimenstion just use LSE
    	bhat = ginv(t(X)%*%X)%*%t(X)%*%Y    
    	residual = Y - X%*%bhat
  	
  	}
  
  #TS_null = matrix(NA, nrow = B, ncol = len(which.covariate))

  ################### HERE WE GO ! ! ! ###########################################
  	

  	###################### This part could be parallelized ##################	
  	count = 1

  	for(wc_cov in which.covariate){   

  		b.Null = bhat
  		#b.Null[wc_cov] = 0
 
  		if(multiTest) { 
  			to.which.covariate = list(wc_cov)
  			to.betaNull = list(betaNull[[count]])

  			b.Null[wc_cov] = betaNull[[count]]


  		}else{
  			to.which.covariate = wc_cov
  			to.betaNull = betaNull[count]

  			b.Null[wc_cov] = betaNull[count]

  		} # then run multiple testing


  		TS_null = Path.Resample.Process(X = X, Y = Y, multiTest = multiTest, residual = residual, b.Null = b.Null, betaNull = to.betaNull, 
  												beta.index = to.which.covariate, B = B, exact = exact, parallel = parallel, ...)

  		rej[count,1] = TS[count] > quantile(TS_null,0.8)
    	rej[count,2] = TS[count] > quantile(TS_null,0.9)
    	rej[count,3] = TS[count] > quantile(TS_null,0.95)
    	rej[count,4] = TS[count] > quantile(TS_null,0.99)
    	pval[count] = mean(TS_null > TS[count])

    	count = count + 1
    	
  	}


  	
  
  ##########################################################

  return(list(rej = rej, pval = pval,TS_null = TS_null, TS = TS))
  
} 

#Path.TS.Para = function(mat, list){
# Calculate PATH statistic exactly, could run this in parallel

#	n = nrow(mat)
#	p = ncol(mat) - 1
#	X = mat[,1:p] 
#	Y = mat[,p+1] 

#	return( do.call(ExactPath.TS, c(X = X,Y = Y, list)) )

#}

Path.Resample.Process = function(X, Y, multiTest, residual, b.Null, beta.index, betaNull, B = 500, exact, parallel = FALSE, ...){
# Bootstrap the null distribution of Path-based statistic of coef beta.index
#
# Args:
#	X, Y, which.covariate : feed in to Path-based TS function
#	B : # of bootstrap replications
#	parallel : run in parallel or nor
#	exact : use exact TS or approx TS
#	beta.index : which coef 
#
# Return:  
#	A vector of the bootstrapped null 
	n = nrow(X)
	p = ncol(X)

	TS_null = numeric()
  
  	if(parallel){ # running in parallel
    	mat = list()


    	for(bs in 1:B){
      		ind = sample(1:n,replace = TRUE)
      		boot_residual = residual[ind]
      		#b_null = bhat
      		#b_null[beta.index] = 0
      		Y = X %*% b.Null + boot_residual
      		mat[[bs]] = cbind(X,Y) 
    	}
  		
  		#rgs = list(...)
    	#Args = c(which.covariate = beta.index, betaNull = betaNull, exact = exact, multiTest = multiTest, args)

    	# On a cluster, just use
    	
    	no_cores <- detectCores() 
    	cat("n_cores detected:", no_cores, "\n")
    	# Initiate cluster
    	#cl <- makeCluster(no_cores)
    	cl <- makeCluster(no_cores, type = "FORK")
    	###### if using window ###### 
      # load special packages
    	#clusterEvalQ(cl, .libPaths("~/R"))
    	#clusterEvalQ(cl, library(glmnet))
    	#clusterEvalQ(cl, library(lars))
    	#clusterEvalQ(cl, library(MASS))
    	#clusterEvalQ(cl,library(pryr))
    	#clusterEvalQ(cl,library(plus))
    	#clusterEvalQ(cl,source("~/hdi_path/bin/pathwise_ts.R"))

    	#clusterExport(cl, varlist = c("beta.index", "exact", "betaNull", "multiTest",...), envir = environment())
      	#clusterExport(cl, varlist = 'Args', envir = environment())
      ###### if using window ######

    	re_list = parLapply(cl, mat, Path.TS.Para, exact = exact, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull, ...)
    	#re_list = parLapply(cl, mat, Path.TS.Para, list = Args)
    
    	######## in case run out of MEMORY
    	print("Cluster MEM:")
    	print(mem_used())
    	########
    	stopCluster(cl)  # END parallel bootsrap


    	for(bss in 1:B){

      		TS_null[bss] = re_list[[bss]]
      
    	}

    	return(TS_null)

    }else{ # not parallel, could be slow

  		for(bs in 1:B){   
    		ind = sample(1:n,replace = TRUE)
    		boot_residual = residual[ind]
    		#b_null = bhat 
    		#b_null[beta.index] = 0
    		Y = X %*% b.Null + boot_residual
    		
    		TS_null[bs] = Path.TS(exact = exact, X = X, Y = Y, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull,...)
	     	
  		}

  		return(TS_null)
  	}
	
    

}

    






	

    
    
   

    	
	
