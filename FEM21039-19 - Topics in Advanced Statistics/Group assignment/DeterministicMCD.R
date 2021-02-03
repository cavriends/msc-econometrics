# --------------------------------------------------------------------
# Author: Corné Vriends
# Group: 25
# Willemijn Brus - 484505
# Martijn Hofstra - 387276
# Brian Leenen - 476908
# CornÃ© Vriends - 440435
# --------------------------------------------------------------------

## Use this code skeleton to implement the deterministic MCD and the plug-in
## robust regression estimator.  Please submit your implementation together
## with your report.

## IMPORTANT: Please do not change any function names and make sure that you
##            use the correct input and output as specified for each function.
##            This will simplify grading because each group's code will be
##            similarly structured.  However, feel free to add other arguments
##            to the function definitions as needed.

#install.packages("robustbase")
library(robustbase) # covOGK()

# User specific code, used for testing
#setwd("/Users/cavriends/Dropbox/ESE/MSc Econometrics/FEM21039-19 - Topics in Advanced Statistics/Assignment - Topics in AS/")
#load("Eredivisie28.RData")
#x = data.matrix(Eredivisie28)

## Functions for initial estimators

# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators

# correlation matrix based on the hyperbolic tangent
corHT <- function(z) {
  s_1 = cor(tanh(z))
  return(s_1)
}

# spearman correlation matrix
corSpearman <- function(z) {
  s_2 = cor(z, method="spearman")
  return(s_2)
}

# correlation matrix based on normal scores of the ranks
corNSR <- function(z) {
  
  p = ncol(z)
  t = c()
  
  for (i in 1:p) {
    t = cbind(t, qnorm(((rank(z[,i]))-1/3)/(nrow(z)+1/3)))
  }
  
  s_3 = cor(t)
  return(s_3)
}

# modified spatial sign covariance matrix
covMSS <- function(z) {
  k = c()
  for (i in 1:nrow(z)) {
    if (sqrt(sum(z[i,]^2)) != 0) {
      k = rbind(k, z[i,]/sqrt(sum(z[i,]^2)))
    } else { # Stond niks over beschreven in de paper wat te doen met 0-vectors, voor nu zo aangepakt.
      k = rbind(k, rep(0, ncol(z)))
      }
  }
  
  s_4 = matrix(data = 0, ncol=ncol(z), nrow=ncol(z))
  
  for (i in 1:nrow(z)) {
    s_4 = s_4 + k[i,] %*% t(k[i,]) 
    if (i == nrow(z)) {
      s_4 = 1/nrow(z)*s_4    
    }
  }
  return(s_4)
}

# covariance matrix based on first step of BACON
covBACON1 <- function(z) {
  
  n_2_ceil = ceiling((nrow(z)/2)) 
  norm_z = apply(z, 1, norm, type="2")
  sorted_z = sort(norm_z, index.return=TRUE)$ix[1:n_2_ceil]
  s_5 = cov(as.matrix(z[sorted_z,]))
  return(s_5)
}

# raw OGK estimator of the covariance matrix with median and Qn
rawCovOGK <- function(z) {
  s_6 = covOGK(z, sigmamu = s_Qn)$cov
  return(s_6)
}

# Own implementation of Qn
Qn_own <- function(x) {
  
  correction_factor = 2.21914
  n = length(x)
  h = floor(n/2+1)
  
  distances = dist(x, method="manhattan", diag=TRUE)
  robust_scatter = correction_factor*sort(distances)[choose(h,2)]
  
  return(robust_scatter)
  
}

## Main function for deterministic MCD algorithm

# Input:
# x ....... data matrix
# alpha ... proportion of observations to be used for subset size
# anything else you need

# Output
# A list with the following components:
# center ....... mean of the reweighted estimator
# cov .......... covariance matrix of the reweighted estimator
# weights ...... binary weights of the observations used for the reweighted
#                estimator (0 if outlier, 1 if used for estimation)
# raw.center ... mean of the raw estimator
# raw.cov ...... covariance matrix of the raw estimator
# best ......... indices of the observations in the best subset found by the
#                raw estimator
# any other output you want to return

covDetMCD <- function(x, alpha, ...) {
  
  # Function is generalized to more dimensions than 2
  n = nrow(x)
  p = ncol(x)

  # Standardize the data
  z = c()
  for (i in 1:p) {
    z = cbind(z, (x[,i] - median(x[,i]))/Qn_own(x[,i]))
  }
  h = h.alpha.n(alpha, nrow(z), ncol(z))

  # Determine initial estimates
  scatter_list = list(corHT(z), corSpearman(z), corNSR(z), covMSS(z), covBACON1(z), rawCovOGK(z))
  number_of_estimates = length(scatter_list)
  
  # General procudure  
  eigen_list = list()
  for (scatter in scatter_list) {
    eigen_list = append(eigen_list, list(eigen(scatter)$vectors))
  }

  cov_list = list()
  for (eigen_matrix in eigen_list) {
    v = z %*% eigen_matrix
    
    diagonal_qn = c()
    for (i in 1:p) {
      diagonal_qn = cbind(diagonal_qn, Qn_own(v[,i])^2)
    }
    
    L = diag(c(diagonal_qn))
    cov_list = append(cov_list, list(eigen_matrix %*% L %*% t(eigen_matrix)))
  }
  
  mu_list = list()
  for(cov in cov_list) {
    mu_list = append(mu_list, list(chol(cov) %*% apply(z %*% chol(solve(cov)), 2, median)))
  }
  
  subset_h_list = list()
  subset_index_list = list() 
  for (i in seq(number_of_estimates)) {
    subset_z = sort(mahalanobis(z, mu_list[[i]], cov_list[[i]]), index.return=TRUE)$ix[1:h]
    subset_index_list = append(subset_index_list, list(subset_z))
    subset_h_list = append(subset_h_list, list(x[subset_z,]))
  }
  
  cov_list_c = lapply(subset_h_list, cov)
  mu_list_c = lapply(subset_h_list, apply, 2, mean)
  
  # C-steps
  convergence = FALSE
  c_steps = 0
  
  while (!convergence) {
    c_steps = c_steps + 1
    det_list = lapply(cov_list_c, det)
    subset_h_list = list()
    subset_index_list = list()
    for (i in seq(number_of_estimates)) {
      sorted_h = sort(mahalanobis(x, mu_list_c[[i]], cov_list_c[[i]]), index.return=TRUE)$ix[1:h]
      subset_index_list = append(subset_index_list, list(sorted_h))
      subset_h_list = append(subset_h_list, list(x[sorted_h,]))
    }
    
    cov_list_c = lapply(subset_h_list, cov)
    mu_list_c = lapply(subset_h_list, apply, 2, mean)
    
    # Check if convergence criteria is met (no change in determinants of covariance matrices for one iteration)
    convergence = all(unlist(det_list) == unlist(lapply(cov_list_c, det)))
    
  }
  
  min_index = min(unlist(lapply(cov_list_c, det)), index.return=TRUE)
  best_index = subset_index_list[[min_index]]
  
  fisher_correction <- function(alpha, p) {
    chi2 <- qchisq(alpha, p)
    correction <- alpha / pgamma(chi2/2, p/2 + 1)
    return(correction)
  }
  
  true_alpha = h/n
  fisher_raw_correction = fisher_correction(true_alpha, p)
  
  rawMCD = fisher_raw_correction*cov_list_c[[min_index]]
  rawMu = mu_list_c[[min_index]]
  
  # Reweighting procedure
  weights = as.numeric(as.matrix(mahalanobis(x, rawMu, rawMCD)) < qchisq(0.975, p))
  
  rwgtMu = apply(x[weights == 1,], 2, mean)
  rwgtMCD = cov(x[weights == 1,])
  
  one_minus_delta = 0.975 # Default value, it is different in covMcd, sum(weights)/nrow(x) to be precise (robustbase)
  fisher_rwgt_correction = fisher_correction(one_minus_delta, p)
  rwgtMCD = fisher_rwgt_correction*rwgtMCD
  
  return_list = list("center" = rwgtMu,
                     "cov" = rwgtMCD,
                     "weights" = weights,
                     "raw.center" = rawMu,
                     "raw.cov" = rawMCD,
                     "best" = best_index)
  
  return(return_list)

}

## Function for regression based on the deterministic MCD

# Input:
# x ........ matrix of explanatory variables
# y ........ response variable
# alpha .... proportion of observations to be used for the subset size in the 
#            MCD estimator
# anything else you need

# Output
# A list with the following components:
# coefficients .... estimated regression coefficients based on the reweighted
#                   deterministic MCD covariance matrix
# fitted.values ... fitted values for all observations in the data
# residuals ....... residuals for all observations in the data
# MCD ............. R object for the deterministic MCD (entire output from
#                   function covDetMCD())
# any other output you want to return

lmDetMCD <- function(x, y, alpha, ...) {
  
  # Generalized to more dimensions
  x = as.matrix(x)
  y = as.matrix(y)
  resultMCD = covDetMCD(cbind(y,x), alpha=alpha)
  p = ncol(x)
  
  estimated_cov = resultMCD$cov
  estimated_center = resultMCD$center
  
  robust_beta = solve(estimated_cov[2:(1+p),2:(1+p)])%*%estimated_cov[2:(1+p),1]
  robust_alpha = estimated_center[1] - t(estimated_center[2:(1+p)])%*%robust_beta
  
  coefficients = c("alpha" = robust_alpha, "beta" = robust_beta)
  fitted_values = t(robust_alpha%*%t(as.matrix(rep(1, nrow(x)), nrow=nrow(x)))) + x%*%robust_beta
  residuals = y - fitted_values
  
  return_list = list("coefficients" = coefficients,
                     "fitted.values" = fitted_values,
                     "residuals" = residuals,
                     "MCD" = resultMCD)
  
  return(return_list)
  
}
