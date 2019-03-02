requrie(lars)
requrie(glmnet)

sourceCpp('LOCO_TS.cpp')




ExtractLars.Path <- function(x, y, whichCov, betaNULL = 0, normalize = TRUE, intercept = FALSE){
  # Extract Path info
  #
  # Args:
  #	X,Y: design matrix and response vector
  #	which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do multiple testing.
  # 
  #	normalize: argguments of lars 
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  # Returns:
  #	A list of lambda vector, original path, and the LOCO path
  #
  #
  n <- nrow(x)
  p <- ncol(x)
  
  ################# Here we go ! ! ! #######################################################  

    ##### extract the path #######

    if (length(whichCov) != length(betaNULL)){

      stop("Length of variables being tested must equal the length of their Null hypothesis")

    }else{

      multiTest <- length(whichCov) > 1

    }

    # simultanoues testing 
    if(multiTest){
      
      #### ajust Y by using Y - X_j * beta_j_NULL - X_k * beta_k_NULL - ...
      adjust.X <- rowSums( t( t(x) * betaNULL  ) )    
      adjust.Y <- y - adjust.X
      
      X_scaled <- scale(x)
      
      
      lars.out <- lars(X_scaled, adjust.Y, type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
      lambda.hat <- sort(c(lars.out$lambda,0),decreasing=FALSE)
      
      
      
      lars.j.out <- lars(X_scaled[, -whichCov], adjust.Y, type = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
      lambda.j.hat <- sort(c(lars.j.out$lambda,0),decreasing=FALSE)
      
      
      
      
    # indivdual testing   
    }else if(!multiTest){  #indivdual test
      
      adjust.Y <- y - betaNULL[l] * x[,j]	
      X_scaled <- scale(x)
      
      lars.out <- lars(X_scaled, adjust.Y , type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
      lambda.hat <- sort(c(lars.out$lambda,0),decreasing = FALSE)
      
   
      lars.j.out <- lars(X_scaled[, -whichCov], adjust.Y, type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
      lambda.j.hat <- sort(c(lars.j.out$lambda,0), decreasing = FALSE)
      
      
    }else{
      stop("wrong input of multiTest, must be boolean")
    }	
    

    #### remove the unoverlapping part (test)
    
    if(lambda.hat[1] != lambda.j.hat[1]){
      
      leftmost = c(lambda.hat[1], lambda.j.hat[1])
      whichone = which.max(leftmost)
      
      if(whichone ==1){
        
        lambda.j.hat <- lambda.j.hat[lambda.j.hat >= lambda.hat[1]]
        lambda.j.hat <- c(lambda.hat[1], lambda.j.hat)

      }else{
        
        lambda.hat <- lambda.hat[lambda.hat >= lambda.j.hat[1]]
        lambda.hat <- c(lambda.j.hat[1], lambda.hat)
        
      }
      
      
      
    }
    
    union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
    
    beta.j.hat <- matrix(0, length(union.lambda), p)
    
    beta.hat <- predict(lars.out, s=union.lambda, type="coef", mode="lambda")$coefficients
    
    beta.j.hat[, -j] <- predict(lars.j.out, s=union.lambda, type="coef", mode="lambda")$coefficients
    
    #beta.hat <- coef(lars.out, lam = union.lambda)
    #beta_val <- coef(lars.j.out, lam = union.lambda)



  return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat))	

}





LOCOLars.TS <- function(obj){
  # Calculate PATH statistic exactly  
  #
  # Args:
  # X,Y: design matrix and response vector
  # which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do multiple testing.
  # 
  # normalize: argguments of lars 
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  # Returns:
  # A list of lambda vector, original path, the LOCO path, and Test Statistic
  #
  #  

  M <- length(obj$union.lambda)
  
  Delta <- obj$beta.j.hat - obj$beta.hat

  Delta_1 <- Delta[-M, ]
  
  Delta_2 <- Delta[-1, ]
  
  Lambda <- diff(obj$union.lambda)
  
  Epsilon <- 1/3 * Lambda * (Delta_1 * Delta_1 + Delta_1 * Delta_2 + Delta_2 * Delta_2)
  
  return(sum(Epsilon))
  
}





Single.Lars.TS <- function(...){
  # Wrapper of ExtractLars.Path and LOCOLars.TS 
  #
  # Args: 
  # ...
  # Return:
  # A list of Test statistic, and a list of lambda vector, original path, the LOCO path
  #
  #

  obj = ExtractLars.Path(...)

  TS = LOCOLars.TS(obj)

  return(list( TS = TS, path.info = obj ) )

 }


Multi.Lars.TS <- function(parallel = TRUE, Cov.List, betaNULL.list, ...){
  # Wrapper of Single.Lars.TS, enable parallel computing performance boost
  # Args: 
  # paralle: boolean, parallel computing or not 
  # Cov.List: a list of variables, which we need to calculate their test statistic
  # Return:
  # A list of Test statistic, and a list of lambda vector, original path, the LOCO path
  #
  #
  #TODO: 
  if (parallel == TRUE){

    continue

  }else{

    continue

  }
  

}





## NOW, GLMNET ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#######################################################
#######################################################
################## GLMNET #############################
################## GLMNET #############################
################## GLMNET #############################
#######################################################
#######################################################



ExtractNet.Path <- function(x, y, whichCov, betaNULL = 0, nlambda = 1000,  alpha = 1, family = "gaussian", normalize = TRUE, intercept = FALSE){
  # Extract Path info
  #
  # Args:
  # X,Y: design matrix and response vector
  # which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do multiple testing.
  # 
  # normalize: argguments of lars 
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  # Returns:
  # A list of lambda vector, original path, and the LOCO path
  #
  #
  n <- nrow(x)
  p <- ncol(x)
  
  ################# Here we go ! ! ! #######################################################  

    ##### extract the path #######

    if (length(whichCov) != length(betaNULL)){

      stop("Length of variables being tested must equal the length of their Null hypothesis")

    }else{

      multiTest <- length(whichCov) > 1

    }

    # simultanoues testing 
    if(multiTest){
      
      #### ajust Y by using Y - X_j * beta_j_NULL - X_k * beta_k_NULL - ...
      adjust.X <- rowSums( t( t(x) * betaNULL  ) )    
      adjust.Y <- as.vector( y - adjust.X )
      X_scaled = scale(x)
            
      ## TODO
      max_lam = NA

      net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = intercept, standardize = normalize)
      lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
      
      
      net.j.out <- glmnet(X_scaled[, -whichCov], adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = intercept, standardize = normalize)
      lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
      
      
      
      
    # indivdual testing   
    }else if(!multiTest){  #indivdual test
      
      adjust.Y <- as.vector( y - betaNULL[l] * x[,j] )
      X_scaled <- scale(x)

      ## TODO
      net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = intercept, standardize = normalize)
      lambda.hat <- sort(net.out$lambda, decreasing = FALSE)

      
      net.j.out <- glmnet(X_scaled[, -whichCov], adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = intercept, standardize = normalize)
      lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
      
      

    }else{
      stop("wrong input of multiTest, must be boolean")
    } 
    
    

    
    minLam <- min(lambda.hat, lambda.j.hat)
    
    lambda.hat = lambda.hat[lambda.hat > minLam]

    lambda.j.hat = lambda.hat[lambda.hat > minLam]

    union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
    
    beta.j.hat <- Matrix(0, length(union.lambda), p)  # sparse Matrix
    
    beta.hat <- predict(net.out, s=union.lambda, type="coef")
    beta.hat <- t( beta.hat[-1, ] )  # we don't need intercept
    

    beta.j.tmp <- predict(net.j.out, s=union.lambda, type="coef")
    beta.j.hat[, -j] <- t( beta.j.tmp[-1, ] ) # we don't need intercep
    

  return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 

}





LOCONet.TS <- function(obj){
  

  M <- length(obj$union.lambda)
  
  Delta <- obj$beta.j.hat - obj$beta.hat

  Delta_1 <- Delta[-M, ]
  
  Delta_2 <- Delta[-1, ]
  
  Lambda <- diff(obj$union.lambda)
  
  Epsilon <- 1/3 * Lambda * (Delta_1 * Delta_1 + Delta_1 * Delta_2 + Delta_2 * Delta_2)
  
  return(sum(Epsilon))
  
 
  
}






Single.Net.TS <- function(...){
  # Wrapper of ExtractLars.Path and LOCOLars.TS 
  #
  # Args: 
  # ...
  # Return:
  # A list of Test statistic, and a list of lambda vector, original path, the LOCO path
  #
  #

  obj = ExtractNet.Path(...)

  TS = LOCONet.TS(obj)

  return(list( TS = TS, path.info = obj ) )

 }


Multi.Net.TS <- function(parallel = TRUE, Cov.List, ...){
  # Wrapper of Single.Lars.TS, enable parallel computing performance boost
  # Args: 
  # paralle: boolean, parallel computing or not 
  # Cov.List: a list of variables, which we need to calculate their test statistic
  # Return:
  # A list of Test statistic, and a list of lambda vector, original path, the LOCO path
  #
  #
  #TODO: 
  if (parallel == TRUE){

    continue

  }else{

    TS = lapply(Cov.List, FUN = Single.Net.TS, ...)

  }
  

}






