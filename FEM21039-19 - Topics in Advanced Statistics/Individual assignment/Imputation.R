# -----------------------------------------
# Author: Corné Vriends - 440435
# -----------------------------------------

## Use this code skeleton to implement the procedures for obtaining point
## estimates and valid standard errors via multiple imputation or the 
## bootstrap.  Please submit your implementation together with your report.

## IMPORTANT: Please do not change any function names and make sure that you
##            use the correct input and output as specified for each function.
##            This will simplify grading because each student's code will be
##            similarly structured.  However, feel free to add other arguments
##            to the function definitions as needed.



## Functions for multiple imputation via iterative model-based imputation

#install.packages("VIM")
#install.packages("cellWise")
#install.packages("robustbase")

# setwd("/Users/cavriends/Dropbox/ESE/MSc Econometrics/FEM21039-19 - Topics in Advanced Statistics/Assignment 2/")

library(VIM)
library(cellWise)
library(robustbase)

# library(MASS)

# Generate test data
# 
# set.seed(12345)
# 
# p = 2
# n = 1000
# percentage_NA = 0.20
# 
# xy = mvrnorm(n = n, mu = rep(0,p), Sigma = matrix(c(1,1/2,1/2,1), ncol = p))
# 
# xy_clean = xy
# colnames(xy_clean) = colnames(xy_clean, do.NULL = FALSE)
# colnames(xy_clean)[ncol(xy_clean)] = "y"
# 
# clean_result = lmrob(y ~ ., data=data.frame(xy_clean))
# 
# set_list = list()
# 
# for (i in seq(p)) {
#   set_list = append(set_list, list(sample(1:n, ceiling(n*percentage_NA))))
# }
# 
# set_list_clean = list()
# for (i in seq(p)) {
#   set_list_clean = append(set_list_clean, list(setdiff(set_list[[i]], intersect(set_list[[1]], set_list[[2]]))))
# }
# 
# for (i in seq(p)) {
#   xy[set_list_clean[[i]],i] = NA
# }


# Multiple imputation
# 
# Input:
# xy    data set with missing values
# m     number of imputations
# DDC   a logical indicating whether to run the DetectDeviatingCells algorithm
#       before performing multiple imputation (such that the flagged outlying 
#       cells are also imputed)
# ...   additional arguments to be passed to function irmi() from package VIM
#       (for example, whether to use robust models)
# 
# Output:
# A list with the following components:
# imputed   this should again be a list in which each list element is one 
#           imputed data set
# m         the number of imputations
# any other information you want to keep

multimp <- function(xy, m = 10, DDC = FALSE, ...) {
  
  p = ncol(xy)
  n = nrow(xy)

  if (DDC == TRUE) {
    
    indices = cellWise::DDC(xy[,1:(p-1)], DDCpars = list("silent" = TRUE))$indall
    vectorized_xy = c(xy[,1:(p-1)])
    vectorized_xy[indices] = NA
    xy = cbind(matrix(vectorized_xy, ncol = p-1, nrow = n), xy[,p])
    
  }
  
  imputed_list = list()
  for (i in seq(m)) {
    imputed_list = append(imputed_list, list(irmi(xy)[,1:p]))
  }
  
  return_list = list("imputed" = imputed_list,
                     "m" = m)
  
  return(return_list)
}

# Fit regression models
# 
# Input:
# xyList   list of imputed data sets as returned by function multimp()
# ...      additional arguments to be passed to modeling function (for example,
#          control parameters for the MM-algorithm)
# 
# Output:
# A list with the following components:
# models   this should again be a list in which each list element is a 
#          regression model (fitted to the corresponding imputed data set)
# m        the number of imputations
# any other information you want to keep

fit <- function(xyList, ...) {
  
  m = xyList$m
  
  regression_list = list()
  for (i in seq(m)) {
    
    xy = xyList$imputed[[i]]
    
    # Set last column to y
    colnames(xy) = colnames(xy, do.NULL = FALSE)
    colnames(xy)[ncol(xy)] = "y"
    
    lmcontrol = lmrob.control(k.max = 1000, maxit.scale = 1000)
    
    regression_list = append(regression_list, list(lmrob(y ~ ., data = xy, control = lmcontrol, na.action = NULL)))
  }
  
  return_list = list("models" = regression_list,
                     "m" = m)
  
  return(return_list)
  
}

# Pool point estimates and standard errors
#
# Input:
# fitList  a list as returned by function fit() containting the regression 
#          models fitted to the imputed data sets
# ...      additional arguments you may need to pass down to other functions
#
# Output:
# A matrix that contains the pooled coefficients in the first column, the
# pooled standard errors in the second column, the t-statistic in the third
# column, the estimated degrees of freedom in the fourth column, and the
# p-value in the fifth column (see slide 50 of Lecture 5 for an example)

pool <- function(fitList, ...) {
  
  m = fitList$m
  modelList = fitList$models
  
  p = length(modelList[[1]]$coefficients)
  n = length(modelList[[1]]$residuals)
  
  # Pooled estimate and Pooled variance
  
  coefficient_list = list()
  covariance_list = list()
  
  summary_matrix = matrix(0, nrow = p, ncol = 5)
  
  for (i in seq(m)) {
    
    model = modelList[[i]]
    
    if (!any(is.na(model$cov))) {
      coefficient_list = append(coefficient_list, list(model$coefficients))
      covariance_list = append(covariance_list, list(diag(model$cov)))
    }
  }
  
  if (!is.null(unlist(coefficient_list))) {
    
    pooled_estimate = apply(matrix(unlist(coefficient_list), ncol = p, byrow = TRUE), 2, mean)
    
    within_imputation_variance = apply(matrix(unlist(covariance_list), ncol = p, byrow = TRUE), 2, mean) 
    between_imputation_variance = diag(cov(matrix(unlist(coefficient_list), ncol = p, byrow = TRUE)))
  
    pooled_variance = within_imputation_variance + ((m+1)/m)*between_imputation_variance
    
    # Standard error and t-statistic
    
    standard_error = c()
    t_statistic = c()
    for (i in seq(p)) {
      standard_error = rbind(standard_error, sqrt(pooled_variance)[i])
      t_statistic = rbind(t_statistic, pooled_estimate[i]/standard_error[i])
    }
    
    # Degrees of freedom and significance level
    gamma = ((m+1)/m) * (between_imputation_variance/pooled_variance)
    v_m = (m-1) / (gamma^2)
    v_comp = n - p
    v_obs = ((v_comp+1)/(v_comp+3))*v_comp*(1-gamma)
    v_tilde = (v_m*v_obs)/(v_m+v_obs)
    
    p_values = 2*pt(abs(t_statistic), v_tilde, lower.tail=FALSE)
    
    summary_matrix[,1] = pooled_estimate
    summary_matrix[,2] = standard_error
    summary_matrix[,3] = t_statistic
    summary_matrix[,4] = v_tilde # Shouldn't it be an integer?
    summary_matrix[,5] = p_values
  
  }
  
  return(summary_matrix)
  
}


## Function for the bootstrap with kNN imputation and linear regression


# Input:
#
# x     data set with missing values
# R     number of bootstrap replications
# k     number of neighbors for kNN imputation
# DDC   a logical indicating whether to run the DetectDeviatingCells algorithm
#       before performing imputation (such that the flagged outlying cells are 
#       also imputed)
# ...   additional arguments to be passed to modeling function (for example,
#       control parameters for the MM-algorithm)
#
# Output:
# A list with the following components:
# replicates   a matrix with all coefficient estimates from all replications
# summary      a matrix that contains the point estimates of the coefficients 
#              in the first column, the standard errors in the second column, 
#              the z-statistic in the third column, and the p-value in the 
#              fourth column (see slide 29 of Lecture 5 for an example)

bootstrap <- function(x, R = 100, k = 10, DDC = FALSE, ...) {

  p = ncol(x)
  n = nrow(x)
  
  bootstrapped_estimates = list()
  
  for (i in seq(R)) {
    x_r = x[sample(seq(n), n, replace = TRUE),]
    
    if (DDC == TRUE) {
      
      indices = cellWise::DDC(x_r[,1:(p-1)], DDCpars = list("silent" = TRUE))$indall
      vectorized_x_r = c(x_r[,1:(p-1)])
      vectorized_x_r[indices] = NA
      x_r = cbind(matrix(vectorized_x_r, ncol = p-1, nrow = n), x_r[,p])
      
    } 
    
    x_r_imputed = kNN(data.frame(x_r), k = k)[,1:p]
    
    # Set last column to y
    colnames(x_r_imputed) = colnames(x_r_imputed, do.NULL = FALSE)
    colnames(x_r_imputed)[ncol(x_r_imputed)] = "y"
    
    lmcontrol = lmrob.control(k.max = 1000, maxit.scale = 1000, na.action = NULL)
    
    bootstrapped_estimates = append(bootstrapped_estimates, list(lmrob(y ~ ., data = data.frame(x_r_imputed), control = lmcontrol)$coefficients))
    
    if (i == R) {
      coefficient_matrix = matrix(unlist(bootstrapped_estimates), ncol = p, byrow = TRUE)
      point_estimates = apply(coefficient_matrix, 2, mean)
      standard_errors = sqrt(diag(cov(coefficient_matrix)))
    }
  }
  
  z_statistic = point_estimates/standard_errors
  p_values = 2*pnorm(abs(z_statistic), lower.tail=FALSE)
  
  summary_matrix = matrix(0, nrow = p, ncol = 4)
  
  summary_matrix[,1] = point_estimates
  summary_matrix[,2] = standard_errors
  summary_matrix[,3] = z_statistic
  summary_matrix[,4] = p_values
  
  return_list = list("replicates" = coefficient_matrix,
                     "summary" = summary_matrix)
  
  return(return_list)
  
}