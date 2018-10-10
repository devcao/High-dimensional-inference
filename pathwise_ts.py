##################################################################################

# Includes all Pathwise test statistics, including end point variant and variate norm 
# This is a python version

##################################################################################



import numpy as numpy
from sklearn import linear_model


def ExactPath_TS(X,Y, which_cov, norm):



def ApproxPath_TS():	

# path 
alphas, _, coefs = linear_model.lars_path(X, y, method='lasso', verbose=True)