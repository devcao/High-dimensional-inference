##################################################################################

# Includes all Pathwise test statistics, including end point variant and variate norm 

##################################################################################



#######
##### loading required packages, functions
#source("~/hello/hdi_path/pathwise_require.R")
#source("~/hello/hdi_path/pathwise_ts.R")
###### 
source("~/hdi_path/bin/pathwise_require.R")

source("~/hdi_path/bin/pathwise_net_ts.R")




Net.Resample = function(X, Y, which.covariate, betaNull, multiTest, alpha, B = 500, parallel = FALSE, exact = TRUE, beta.init = 'adaptive', beta.true = beta, ...){
# Bootstrap the null distribution of Path-based statistic, and return reject or not, 
# now for test for a single test only 
#
# Args:
# X, Y, which.covariate : feed in to Path-based TS function
# B : # of bootstrap replications
# parallel : run in parallel or nor
# exact : use exact TS or approx TS
# beta.true : for simulation only, the true value of beta would be used as initial estimates.
# 
# Return:  
# Reject or not under alpha = 0.2,0.1,0.05,0.01
#   p.values of the test
  
  
  n = nrow(X)
  p = ncol(X)
  rej = matrix(0,len(which.covariate),4)  # 1st : which cov, 2nd: na, 3rd: which alpha, 0.2,0.1,0.05,0.01
  pval = numeric()

  results = list()
  results$L2.squared = list(rej = matrix(0,len(which.covariate),4), pval = numeric(), TS_null = numeric(), TS = numeric())
  results$L2 = list(rej = matrix(0,len(which.covariate),4), pval = numeric(), TS_null = numeric(), TS = numeric())
  results$L1 = list(rej = matrix(0,len(which.covariate),4), pval = numeric(), TS_null = numeric(), TS = numeric())
  results$L_inf = list(rej = matrix(0,len(which.covariate),4), pval = numeric(), TS_null = numeric(), TS = numeric())


  TS = Net.TS(exact = exact, X = X, Y = Y, alpha = alpha, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest,...)[[1]]
  
  



  if(p >= n){ # high dimension we can try...

    if(beta.init == "adaptive"){
      bhat = adalasso(X = X, y = Y, k = 10, use.Gram = FALSE,both = TRUE, intercept = FALSE)$coefficients.adalasso
      
    }else if (beta.init == "MC+"){
      bhat = coef(cv.ncvreg(X = X, y = Y, penalty = "MCP",family = "gaussian", nfold= 10))[-1]
          
    }else if (beta.init == "SCAD"){
      bhat = coef(cv.ncvreg(X = X, y = Y, penalty = "SCAD",family = "gaussian", nfold= 10))[-1]

    }else if (beta.init == "Truth"){
      bhat = beta.true
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

#############
#############
#############
#############
#############
      TS_null = Net.Resample.Process(X = X, Y = Y, alpha = alpha, multiTest = multiTest, residual = residual, b.Null = b.Null, betaNull = to.betaNull, 
                          beta.index = to.which.covariate, B = B, exact = exact, parallel = parallel, ...)




      for (iter_list in 1:4){
        
        results[[iter_list]]$rej[count,1] = TS[[iter_list]] > quantile(TS_null[[iter_list]],0.8)
        results[[iter_list]]$rej[count,2] = TS[[iter_list]] > quantile(TS_null[[iter_list]],0.9)
        results[[iter_list]]$rej[count,3] = TS[[iter_list]] > quantile(TS_null[[iter_list]],0.95)
        results[[iter_list]]$rej[count,4] = TS[[iter_list]] > quantile(TS_null[[iter_list]],0.99)
        results[[iter_list]]$pval[count] = mean(TS_null[[iter_list]] >= TS[[iter_list]])

        results[[iter_list]]$TS = TS[[iter_list]]
        results[[iter_list]]$TS_null = TS_null[[iter_list]]
        


      }
      
      count = count + 1


#############
#############
#############
#############
#############


      
    }


    
  
  ##########################################################

  return(results)
  
} 

#Path.TS.Para = function(mat, list){
# Calculate PATH statistic exactly, could run this in parallel

# n = nrow(mat)
# p = ncol(mat) - 1
# X = mat[,1:p] 
# Y = mat[,p+1] 

# return( do.call(ExactPath.TS, c(X = X,Y = Y, list)) )

#}

Net.Resample.Process = function(X, Y, multiTest, alpha, residual, b.Null, beta.index, betaNull, B = 500, exact, parallel = FALSE, ...){
# Bootstrap the null distribution of Path-based statistic of coef beta.index
#
# Args:
# X, Y, which.covariate : feed in to Path-based TS function
# B : # of bootstrap replications
# parallel : run in parallel or nor
# exact : use exact TS or approx TS
# beta.index : which coef 
#
# Return:  
# A vector of the bootstrapped null 
  n = nrow(X)
  p = ncol(X)

  TS_null = list()
  TS_null$L2.squared = numeric()
  TS_null$L2 = numeric()
  TS_null$L1 = numeric()
  TS_null$L_inf = numeric()
  
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

      re_list = parLapply(cl, mat, Net.TS.Para, exact = exact, alpha = alpha, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull, ...)
      #re_list = parLapply(cl, mat, Path.TS.Para, list = Args)
    
      ######## in case run out of MEMORY
      print("Cluster MEM:")
      print(mem_used())
      ########
      stopCluster(cl)  # END parallel bootsrap


      for(bss in 1:B){

          
          TS_null$L2.squared[bss] = re_list[[bss]][[1]]$L2.squared
          TS_null$L2[bss] = re_list[[bss]][[1]]$L2
          TS_null$L1[bss] = re_list[[bss]][[1]]$L1
          TS_null$L_inf[bss] = re_list[[bss]][[1]]$L_inf
            



          
      }

      return(TS_null)

    }else{ # not parallel, could be slow

      for(bs in 1:B){   
        ind = sample(1:n,replace = TRUE)
        boot_residual = residual[ind]
        #b_null = bhat 
        #b_null[beta.index] = 0
        Y = X %*% b.Null + boot_residual
        
        result = Net.TS(exact = exact, alpha = alpha, X = X, Y = Y, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull,...)[[1]]
        
        TS_null$L2.squared[bs] = result$L2.squared
        TS_null$L2[bs] = result$L2
        TS_null$L1[bs] = result$L1
        TS_null$L_inf[bs] = result$L_inf
          

      }

      return(TS_null)
    }
  
    

}

    





  
    
   

      
  
