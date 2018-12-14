# Includes all Pathwise test statistics, including end point variant and variate norm 

##################################################################################


#######
##### loading required packages
source("~/hdi_path/bin/pathwise_require.R")
source("~/hdi_path/bin/pathwise_ts.R")
source("~/hdi_path/bin/pathwise_net_ts.R")
source("~/hdi_path/bin/pathwise_resample.R")
source("~/hdi_path/bin/pathwise_net_resample.R")

###### 



Path.Screening <- function(beta.index = NULL, X, Y, which.covariate, betaNull, multiTest, B = 500, parallel = TRUE, exact = TRUE, beta.true = beta, ...){
# Calculate PATH statistic exactly	
#
# Args:
#	X,Y: design matrix and response vector
#	which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do multiple testing.
# 	nrom: indludes "L1", "L2", "L2.squared","L_inf"
# 	path.method: includes "lars" and "plus". "lars" would give full Lasso path and "plus" cound do MC+ and SCAD path
#	normalize: argguments of lars 
# 	ridgePen: ridge penalty term, 0 indicates no ridge penalty
# 	betaNull: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNull. 'Zero' or user sepcify,
#		like betaNull = c(1,1,1), which.covariate = c(1,2,3),betaNull = list(c(1,1,1),c(0,0)), which.covariate = list(c(1,2,3) ,c(5,6))
# Returns:
#	A vector of PATH statistic
#
#
	n = nrow(X)
	p = ncol(X)

	mat = cbind(X, Y)
   		
	if(parallel){
		
		no_cores <- detectCores() 
    	cat("n_cores detected:", no_cores, "\n")
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

    	#re_list = parLapply(cl, as.list(1:p), Path.TS, X = X, Y = Y, exact = exact, multiTest = multiTest, betaNull = 0, ...)
   		re_list = parLapply(cl, as.list(1:p), Path.TS.Para, mat = mat, exact = exact, multiTest = multiTest, betaNull = 0, ...)
   		
   		stopCluster(cl)

   		TS = unlist(re_list)
   		#re_list = parLapply(cl, mat, Path.TS.Para, list = Args)
    

	}else{

		TS = Path.TS(X = X, Y = Y , which.covariate = 1:p, betaNull = rep(0,p), multiTest = FALSE, ...)

	}
	

	
	print(TS)
	screen.out = which(TS == 0)
	screen.in = which(TS != 0)

	if(which.covariate %in% screen.out){

		pval = 1
		rej = matrix(0,1,4)

	}else{

		X_rdc = X[, screen.in]

		print(screen.in)

		results = Path.Resample(X = X_rdc, Y = Y, which.covariate = which.covariate, 
								betaNull = betaNull, multiTest = multiTest, B = B, parallel = parallel, 
								exact = exact, beta.init = 'adaptive', beta.true = beta.true, ...)
		pval = results$pval
		rej = results$rej		

	}

	print(pval)

	if( is.null(beta.index) ){
		return(list(pval = pval, rej = rej)) 	

	}else{

		IR = all( beta.index %in% which(TS > 0) )
		print(IR)

		return(list(pval = pval, rej = rej, IR = IR)) 		


	}
	
}






Net.Screening <- function(beta.index = NULL, X, Y, alpha = 1, which.covariate, betaNull, multiTest, B = 500, parallel = TRUE, exact = TRUE, beta.true = beta, ...){
# Calculate PATH statistic exactly	
#
# Args:
#	X,Y: design matrix and response vector
#	which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do multiple testing.
# 	nrom: indludes "L1", "L2", "L2.squared","L_inf"
# 	path.method: includes "lars" and "plus". "lars" would give full Lasso path and "plus" cound do MC+ and SCAD path
#	normalize: argguments of lars 
# 	ridgePen: ridge penalty term, 0 indicates no ridge penalty
# 	betaNull: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNull. 'Zero' or user sepcify,
#		like betaNull = c(1,1,1), which.covariate = c(1,2,3),betaNull = list(c(1,1,1),c(0,0)), which.covariate = list(c(1,2,3) ,c(5,6))
# Returns:
#	A vector of PATH statistic
#
#
	n = nrow(X)
	p = ncol(X)

	mat = cbind(X, Y)
   		
	if(parallel){
		
		no_cores <- detectCores() 
    	cat("n_cores detected:", no_cores, "\n")
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

    	#re_list = parLapply(cl, as.list(1:p), Path.TS, X = X, Y = Y, exact = exact, multiTest = multiTest, betaNull = 0, ...)
   		re_list = parLapply(cl, as.list(1:p), Net.TS.Para, mat = mat, alpha = alpha, exact = exact, multiTest = multiTest, betaNull = 0, ...)
   		
   		stopCluster(cl)

   		TS = unlist(re_list)
   		#re_list = parLapply(cl, mat, Path.TS.Para, list = Args)
    

	}else{

		TS = Net.TS(X = X, Y = Y, alpha = alpha, which.covariate = 1:p, betaNull = rep(0,p), multiTest = FALSE, ...)

	}
	

	
	print(TS)
	screen.out = which(TS == 0)
	screen.in = which(TS != 0)

	if(which.covariate %in% screen.out){

		pval = 1
		rej = matrix(0,1,4)

	}else{

		X_rdc = X[, screen.in]

		print(screen.in)

		results = Net.Resample(X = X_rdc, Y = Y, alpha = alpha, which.covariate = which.covariate, 
								betaNull = betaNull, multiTest = multiTest, B = B, parallel = parallel, 
								exact = exact, beta.init = 'adaptive', beta.true = beta.true, ...)
		pval = results$pval
		rej = results$rej		

	}

	print(pval)

	if( is.null(beta.index) ){
		return(list(pval = pval, rej = rej)) 	

	}else{

		IR = all( beta.index %in% which(TS > 0) )
		print(IR)

		return(list(pval = pval, rej = rej, IR = IR, TS = TS)) 		


	}
	
}


