##################################################################################

# Includes all Pathwise test statistics, including end point variant and variate norm 

##################################################################################


#######
##### loading required packages
source("~/hdi_path/bin/pathwise_require.R")
###### 

Path.TS <- function(exact = TRUE,...){

	if(exact){
		return(ExactPath.TS(...))

	}else{
		return(ApproxPath.TS(...))

	}

}

Path.TS.Para <- function(exact = TRUE,...){

	if(exact){
		return(ExactPath.TS.Para(...))

	}else{
		return(ApproxPath.TS.Para(...))

	}

}


ridge.Stacking <- function(X, Y, ridgePen){
# Stacking iden matrix to utilize Elastic Net
#
#	
	n = nrow(X)
	p = ncol(X)
	zero_vec = matrix(0, p, 1)
	Iden_mat = diag(rep(1,p))
	X = rbind(X,ridgePen*Iden_mat) 
	Y = rbind(Y,zero_vec)
	
	return(list(X = X, Y = Y))
}




ExactPath.TS <- function(X,Y, which.covariate, betaNull, multiTest, path.method = 'lars', norm = 'L2.squared', normalize = TRUE, intercept = FALSE, ridgePen = 0){
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
	TS <- numeric(length(which.covariate))

	for (j in which.covariate){

 		# first condition : whether test all 0 or nor
 		
		if(path.method == "lars"){
				# 2nd condition, muti-test or not 
			if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)


					lars.out <- lars(X.sc, newY, type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(c(lars.out$lambda,0),decreasing=FALSE)
					beta.hat <- coef(lars.out)[seq(length(lambda.hat),1,-1),] 

					lars.j.out <- lars(X.sc[, -j], newY, type = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(c(lars.j.out$lambda,0),decreasing=FALSE)
					beta_val <- coef(lars.j.out)[seq(length(lambda.j.hat),1,-1),]

			
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)

					lars.out <- lars(X.sc, newY , type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(c(lars.out$lambda,0),decreasing=FALSE)
					beta.hat <- coef(lars.out)[seq(length(lambda.hat),1,-1),] 

					lars.j.out <- lars(X.sc[, -j], newY, type = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(c(lars.j.out$lambda,0),decreasing=FALSE)
					beta_val <- coef(lars.j.out)[seq(length(lambda.j.hat),1,-1),]

				}else{
					stop("wrong input")
				}	

			}else if(path.method == "plus.lasso"){


				if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)

					#### testing
					lars.out <- plus(X.sc, as.vector(newY), method = "lasso",intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path, decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path, decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

			
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)


					lars.out <- plus(X.sc, as.vector(newY) , method = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path, decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path, decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

				}else{
					stop("wrong input")
				}	

				
			}else if (path.method == "plus.mc+"){
				if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)


					lars.out <- plus(X.sc, as.vector(newY), method = "mc+", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)
					
					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "mc+", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)
				
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)

					lars.out <- plus(X.sc, as.vector(newY) , method = "mc+", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "mc+", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

				}else{
					stop("wrong input")
				}




			}else if (path.method == "plus.scad"){

				if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)


					lars.out <- plus(X.sc, as.vector(newY), method = "scad", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "scad", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

					
					
			
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)

					lars.out <- plus(X.sc, as.vector(newY) , method = "scad", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "scad", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)
				
				}else{
					stop("wrong input")
				}




			}

		#### added 0 penalty value (test)
		#lambda.hat = c(0,lambda.hat)
		#lambda.j.hat = c(0,lambda.j.hat)
		#zero.lars.pen = lars(X.sc, newY, type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)			
		#beta.hat[1,] = predict(a,zero.lars.pen=0,mode='lambda',type = "coefficients")$coefficients
		#zero.lars.j.pen = lars(X.sc[,-j], newY, type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)			
		#beta_val[1,] = predict(a,zero.lars.j.pen=0,mode='lambda',type = "coefficients")$coefficients
		#####

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
			
			beta.hat <- coef(lars.out, lam = lambda.hat)
			beta_val <- coef(lars.j.out, lam = lambda.j.hat)
					
			
		}

		####

	

		new_beta <- matrix(0, dim(beta_val)[1], p)
		new_beta[, -j] <- beta_val
		beta.j.hat <- new_beta
 
		union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)),decreasing=FALSE)  # union of lambda values from both solution paths.
		M <- length(union.lambda) # cardinality of union of lambda values from both solution paths
 		
		beta.hat.union.lambda <- beta.j.hat.union.lambda <- matrix(NA,length(union.lambda),p)
 
		TS.k <- numeric()

 ############# issue here ! ! ! ########################
		for (k in 1:p){
			
			# get beta.hat and beta.j.hat at all values of lambda in union lambda:
			beta.hat.union.lambda[,k] <- 
			approx(x = lambda.hat, y = beta.hat[,k], xout = union.lambda,yright=0)$y 
			
			beta.j.hat.union.lambda[,k] <- 
			approx(x = lambda.j.hat, y = beta.j.hat[,k], xout = union.lambda,yright=0)$y
			# get absolute difference between beta.hat and beta.j.hat at all values of lambda in union lambda
			
			delta <- (beta.hat.union.lambda[,k] - beta.j.hat.union.lambda[,k])
				
			if (norm == "L2.squared"){
				TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
			
			}else if (norm == "L1"){
				TS.k[k] <- 0.5*sum(diff(union.lambda)*(abs(delta[-M]) + abs(delta[-1])))

			}else if (norm == "L2"){
				TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))

			}else if (norm == "L_inf"){
				TS.k[k] <- max(abs(delta))

			}
			

		}
 		
 		if (norm == "L_inf"){
 			TS[l] <- max(TS.k)

 		}else if (norm == "L2"){
 			TS[l] <- sqrt(sum(TS.k))

 		}else{
 			TS[l] <- sum(TS.k)

 		}
		
		l <- l + 1
 
	}
 
	return(TS)
 
}









ApproxPath.TS <- function(X,Y, which.covariate, betaNull, multiTest, path.method = "lars", norm = 'L2.squared',normalize = TRUE, intercept = FALSE, ridgePen = 0, numInterpolation = 10000){
# Calculate PATH statistic using interpolating approximation 
#
# Args:
#	X,Y: design matrix and response vector
#	which.covariate: a vector, indicating which covariate we will be computing
# 	normalize: argguments of lars 
# 	nrom: indludes "L1", "L2", "L2.squared","L_inf"
# 	ridgePen: ridge penalty term, 0 indicates no ridge penalty
# 	numInterpolation: number of interpolations, default 10000.
#betaNull: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNull. 'Zero' or user sepcify,
#		like betaNull = c(1,1,1), which.covariate = c(1,2,3),betaNull = list(c(1,1,1),c(0,0)), which.covariate = list(c(1,2,3) ,c(5,6))

# Returns:
#	A vector of PATH statistic
#
#	

	n = nrow(X)
	p = ncol(X)
	
	if (ridgePen > 0){
		zero_vec = rep(0,p)
		Iden_mat = diag(rep(1,p))
		X = rbind(X,ridgePen*Iden_mat) 
		Y = append(Y,zero_vec)
	
	}
################# Here we go ! ! ! ##############################################



	l <- 1
	TS <- numeric(length(which.covariate))


	for (j in which.covariate){
 
		# first condition : whether test all 0 or nor
 		
		if(path.method == "lars"){
				# 2nd condition, muti-test or not 
			if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)


					lars.out <- lars(X.sc, newY, type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(c(lars.out$lambda,0),decreasing=FALSE)
					beta.hat <- coef(lars.out)[seq(length(lambda.hat),1,-1),] 

					lars.j.out <- lars(X.sc[, -j], newY, type = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(c(lars.j.out$lambda,0),decreasing=FALSE)
					beta_val <- coef(lars.j.out)[seq(length(lambda.j.hat),1,-1),]

			
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)

					lars.out <- lars(X.sc, newY , type = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(c(lars.out$lambda,0),decreasing=FALSE)
					beta.hat <- coef(lars.out)[seq(length(lambda.hat),1,-1),] 

					lars.j.out <- lars(X.sc[, -j], newY, type = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(c(lars.j.out$lambda,0),decreasing=FALSE)
					beta_val <- coef(lars.j.out)[seq(length(lambda.j.hat),1,-1),]

				}else{
					stop("wrong input")
				}	

			}else if(path.method == "plus.lasso"){


				if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)

					#### testing
					lars.out <- plus(X.sc, as.vector(newY), method = "lasso",intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path, decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path, decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

			
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)


					lars.out <- plus(X.sc, as.vector(newY) , method = "lasso", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path, decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "lasso", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path, decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

				}else{
					stop("wrong input")
				}	

				
			}else if (path.method == "plus.mc+"){
				if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)


					lars.out <- plus(X.sc, as.vector(newY), method = "mc+", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)
					
					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "mc+", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)
				
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)

					lars.out <- plus(X.sc, as.vector(newY) , method = "mc+", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "mc+", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

				}else{
					stop("wrong input")
				}




			}else if (path.method == "plus.scad"){

				if(multiTest & is.list(which.covariate) & is.list(betaNull)){
					

					adjust.X = rowSums(  t( apply(X[,j], 1, function(x){return(x*betaNull[[l]])}) )  )
					newY = Y - adjust.X

					X.sc = scale(X)


					lars.out <- plus(X.sc, as.vector(newY), method = "scad", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "scad", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)

					
					
			
				}else if((!multiTest)){  #indivdual test
						
					newY = Y - betaNull[l] * X[,j]	
					X.sc =scale(X)

					lars.out <- plus(X.sc, as.vector(newY) , method = "scad", intercept = intercept, use.Gram = FALSE, normalize=normalize)
					lambda.hat <- sort(lars.out$lam.path,decreasing=FALSE)
					beta.hat <- coef(lars.out, lam = lambda.hat)

					lars.j.out <- plus(X.sc[, -j], as.vector(newY), method = "scad", intercept = intercept, use.Gram = FALSE,normalize=normalize)
					lambda.j.hat <- sort(lars.j.out$lam.path,decreasing=FALSE)
					beta_val <- coef(lars.j.out, lam = lambda.j.hat)
				
				}else{
					stop("wrong input")
				}




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
			
			beta.hat <- coef(lars.out, lam = lambda.hat)
			beta_val <- coef(lars.j.out, lam = lambda.j.hat)
					
			
		}

		####


		new_beta <- matrix(0, dim(beta_val)[1], p)
		new_beta[, -j] <- beta_val
		beta.j.hat <- new_beta
 
		union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)),decreasing=FALSE)  # union of lambda values from both solution paths.
		M <- length(union.lambda) # cardinality of union of lambda values from both solution paths
 
		beta.hat.union.lambda <- beta.j.hat.union.lambda <- matrix(NA,length(union.lambda),p)
 
		TS.k <- numeric()
		
 		for (k in 1:p){

 			# using logrithmic scale: 
			#interpolating.lambda = exp(seq(log(min(union.lambda[union.lambda > 0])),log(max(union.lambda)),length.out = numInterpolation+1))   
			# using linear scale: 
			interpolating.lambda = seq(min(union.lambda),max(union.lambda),length.out = numInterpolation+1)   
		
			approx.beta.hat = approx(x = lambda.hat, y = beta.hat[,k], xout = interpolating.lambda,yright = 0)
			approx.beta.j.hat = approx(x = lambda.j.hat, y = beta.j.hat[,k], xout = interpolating.lambda,yright = 0)

			#beta.hat.union.lambda[,k] <- approx.beta.hat$y 
			#beta.j.hat.union.lambda[,k] <- approx.beta.j.hat$y
			# get absolute difference between beta.hat and beta.j.hat at all values of lambda in union lambda
			
			if (norm == "L2.squared"){
				#delta <- (beta.hat.union.lambda[,k] - beta.j.hat.union.lambda[,k])
				#TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
				
				TS.k[k] <- sum(((approx.beta.hat$y - approx.beta.j.hat$y)^2)[-1] * diff(interpolating.lambda))

			}else if (norm == "L1" ){
				TS.k[k] <- sum((abs(approx.beta.hat$y - approx.beta.j.hat$y))[-1] * diff(interpolating.lambda))

			}else if (norm == "L2"){
				TS.k[k] <- sum(((approx.beta.hat$y - approx.beta.j.hat$y)^2)[-1] * diff(interpolating.lambda))

			}else if (norm == "L_inf"){
				TS.k[k] <- max(abs(approx.beta.hat$y - approx.beta.j.hat$y))
				

			}
			

		}

		

		if (norm == "L2.squared" | norm == "L1"){
			TS[l] <- sum(TS.k)	
		}else if (norm == "L2"){
			TS[l] = sqrt(sum(TS.k))
		}else if (norm == "L_inf"){
			TS[l] <-  max(TS.k)
		}

 
		
 
		l <- l + 1





 
	}
 
	return(TS)
 
}




ExactPath.TS.Para <- function(mat, ...){
# Calculate PATH statistic exactly, could run this in parallel

	n = nrow(mat)
	p = ncol(mat) - 1
	X = mat[,1:p] 
	Y = mat[,p+1] 

	return(ExactPath.TS(X,Y,...))

}

ApproxPath.TS.Para <- function(mat, ...){
# Calculate PATH statistic using interpolating, could run this in parallel

	n = nrow(mat)
	p = ncol(mat) - 1
	X = mat[,1:p] 
	Y = mat[,p+1] 

	return(ApproxPath.TS(X,Y,...))

}







