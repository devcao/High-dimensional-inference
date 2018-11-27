
Net.TS <- function(exact = TRUE, ...){

	if(exact){
		return(ExactNet.TS(...))

	}else{
		return(ApproxNet.TS(...))

	}	

	
}




Net.TS.Para <- function(exact = TRUE,...){

	if(exact){
		return(ExactNet.TS.Para(...))

	}else{
		return(ApproxNet.TS.Para(...))

	}

}


ExactNet.TS <- function(X,Y, alpha = 1, which.covariate, betaNull, multiTest, normalize = TRUE, intercept = FALSE, n_grid = 1000, plot.Diff = FALSE){
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


################# Here we go ! ! ! ##############################################
	


	l <- 1
	
	TS <- list()
		


	for (j in which.covariate){

		#condition, muti-test or not 
		if( multiTest & is.list(which.covariate) & is.list(betaNull) ){
					

			adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
			newY = Y - adjust.X

			X.sc = scale(X)

			a1 = glmnet(X.sc, newY, alpha = alpha, nlambda = 100, family = "gaussian", intercept = intercept, standardize = normalize)
			max_lam = max(a1$lambda)

			#union.lambda = seq(0.001, max_lam + 1, 0.001)
			union.lambda = exp( seq(-5, log(max_lam) + 0.01, length.out = n_grid) )


			full_net = glmnet(X.sc, newY, alpha = alpha, family = "gaussian", intercept = intercept, lambda = union.lambda, standardize = normalize)
			reduce_net = glmnet(X.sc[,-j], newY, alpha = alpha, family = "gaussian", intercept = intercept, lambda = union.lambda, standardize = normalize)

			beta.hat = coef(full_net)[-1,]  # n_beta * n_lambda, not like lars.
			#lambda.hat = sort(full_net$lambda, decreasing = FALSE)


			
			beta.j.hat = Matrix(0, nrow = dim(beta.hat)[1], ncol = dim(beta.hat)[2])
			beta.j.hat[-j, ] = (coef(reduce_net)[-1,])


			if(plot.Diff){

				plotTwo(full_net, reduce_net)
			}	

			
		}else if( (!multiTest) ){  #indivdual test
						
			newY = Y - betaNull[l] * X[,j]
			X.sc = scale(X)


			a1 = glmnet(X.sc, newY, alpha = alpha, nlambda = 100, family = "gaussian", intercept = intercept, standardize = normalize, )
			max_lam = max(a1$lambda)

			#union.lambda = seq(0.01, max_lam + 1, length.out = n_grid)

			
			union.lambda = exp( seq(-5, log(max_lam) + 0.01, length.out = n_grid) )
	
			
			
			

			full_net = glmnet(X.sc, newY, alpha = alpha, family = "gaussian", intercept = intercept, lambda = union.lambda, standardize = normalize)
			reduce_net = glmnet(X.sc[,-j], newY, alpha = alpha, family = "gaussian", intercept = intercept, lambda = union.lambda, standardize = normalize)

			beta.hat = (coef(full_net)[-1,])   # n_beta * n_lambda
			#lambda.hat = sort(full_net$lambda, decreasing = FALSE)


			beta.j.hat = Matrix(0, nrow = dim(beta.hat)[1], ncol = dim(beta.hat)[2])
			beta.j.hat[-j, ] = (coef(reduce_net)[-1,])


			if(plot.Diff){

				plotTwo(full_net, reduce_net)
			}	




		}else{

			stop("wrong input")

		}	

		M <- dim(beta.hat)[2] # cardinality of union of lambda values from both solution paths
 		
		

		TS.k <- list()
		TS.k$L2.squared <- numeric()
		TS.k$L2 <- numeric()
		TS.k$L1 <- numeric()
		TS.k$L_inf <- numeric()
		
		for (k in 1:p){
			
			delta <- rev(beta.hat[k,]) - rev(beta.j.hat[k,])
			
			TS.k$L2[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
			TS.k$L1[k] <- 0.5*sum(diff(union.lambda)*(abs(delta[-M]) + abs(delta[-1])))
			#TS.k$L2[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
			TS.k$L_inf[k] <- max(abs(delta))


		}
 		
		TS[[l]] = list()

		TS[[l]]$L2.squared <- sum(TS.k$L2)
		TS[[l]]$L2 <- sqrt(sum(TS.k$L2))
		TS[[l]]$L1 <- sum(TS.k$L1)
		TS[[l]]$L_inf <- max(TS.k$L_inf)
		

		
		l <- l + 1
 

 
	}
 
	return(TS)
 
}













ApproxNet.TS <- function(X,Y, alpha = 1, which.covariate, betaNull, multiTest, normalize = TRUE, intercept = FALSE, n_grid = 1000){
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


################# Here we go ! ! ! ##############################################
	


	l <- 1
	
	TS <- list()
		


	for (j in which.covariate){

		#condition, muti-test or not 
		if( multiTest & is.list(which.covariate) & is.list(betaNull) ){
					

			adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
			newY = Y - adjust.X

			X.sc = scale(X)

			a1 = glmnet(X.sc, newY, alpha = alpha, nlambda = 100, family = "gaussian", intercept = intercept, standardize = normalize)
			max_lam = max(a1$lambda)

			#union.lambda = seq(0.001, max_lam + 1, 0.001)
			union.lambda = exp( seq(-5, log(max_lam) + 0.01, length.out = n_grid) )


			full_net = glmnet(X.sc, newY, alpha = alpha, family = "gaussian", intercept = intercept, lambda = union.lambda, standardize = normalize)
			reduce_net = glmnet(X.sc[,-j], newY, alpha = alpha, family = "gaussian", intercept = intercept, lambda = union.lambda, standardize = normalize)

			beta.hat = as.matrix(coef(full_net)[-1,])   # n_beta * n_lambda, not like lars.
			#lambda.hat = sort(full_net$lambda, decreasing = FALSE)


			beta.j.hat = matrix(0, nrow = dim(beta.hat)[1], ncol = dim(beta.hat)[2])
			beta.j.hat[-j, ] = as.matrix(coef(reduce_net)[-1,])



			
		}else if( (!multiTest) ){  #indivdual test
						
			newY = Y - betaNull[l] * X[,j]
			X.sc = scale(X)


			a1 = glmnet(X.sc, newY, alpha = alpha, nlambda = 100, family = "gaussian", intercept = intercept, standardize = normalize, )
			max_lam = max(a1$lambda)

			#union.lambda = seq(0.01, max_lam + 1, 0.01)
			union.lambda = exp( seq(-5, log(max_lam) + 0.01, length.out = n_grid) )


			full_net = glmnet(X.sc, newY, alpha = alpha, family = "gaussian", intercept = intercept, lambda = union.lambda, standardize = normalize)
			reduce_net = glmnet(X.sc[,-j], newY, alpha = alpha, family = "gaussian", intercept = intercept, lambda = union.lambda, standardize = normalize)

			beta.hat = (coef(full_net)[-1,])   # n_beta * n_lambda
			#lambda.hat = sort(full_net$lambda, decreasing = FALSE)


			beta.j.hat = Matrix(0, nrow = dim(beta.hat)[1], ncol = dim(beta.hat)[2])
			beta.j.hat[-j, ] = (coef(reduce_net)[-1,])



		}else{

			stop("wrong input")

		}	

		TS.k <- list()
		TS.k$L2.squared <- numeric()
		TS.k$L2 <- numeric()
		TS.k$L1 <- numeric()
		TS.k$L_inf <- numeric()
		
		for (k in 1:p){

			delta <- rev(beta.hat[k,]) - rev(beta.j.hat[k,])
			
			TS.k$L2[k] <- sum(((delta)^2)[-1] * diff(union.lambda))
			TS.k$L1[k] <- sum((abs(delta))[-1] * diff(union.lambda))
			#TS.k$L2[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
			TS.k$L_inf[k] <- max(abs(delta))


		}
 		
		TS[[l]] = list()

		TS[[l]]$L2.squared <- sum(TS.k$L2)
		TS[[l]]$L2 <- sqrt(sum(TS.k$L2))
		TS[[l]]$L1 <- sum(TS.k$L1)
		TS[[l]]$L_inf <- max(TS.k$L_inf)
		
	
		l <- l + 1
 

 
	}
 
	#return(list(TS = TS, TS.k = TS.k))
 	return(TS)
 
}



ExactNet.TS.Para <- function(mat, ...){
# Calculate PATH statistic exactly, could run this in parallel

	n = nrow(mat)
	p = ncol(mat) - 1
	X = mat[,1:p] 
	Y = mat[,p+1] 

	return(ExactNet.TS(X,Y,...))

}

ApproxNet.TS.Para <- function(mat, ...){
# Calculate PATH statistic using interpolating, could run this in parallel

	n = nrow(mat)
	p = ncol(mat) - 1
	X = mat[,1:p] 
	Y = mat[,p+1] 

	return(ApproxNet.TS(X,Y,...))

}





plotTwo = function(full_net, reduce_net){
  par(mfrow=c(1, 2))

  plot(full_net, xvar = "lambda", label = TRUE)
  plot(reduce_net, xvar = "lambda", label = TRUE)


}


  

