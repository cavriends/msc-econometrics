# -----------------------------------------
# Author:
# Corné Vriends - 440435
# -----------------------------------------

## Use this empty file to implement your simulation study.  Please submit your 
## code together with your report.

library(MASS)
library(cellWise)

# Load R script with your implementation of the methods

#setwd("/Users/cavriends/Dropbox/ESE/MSc Econometrics/FEM21039-19 - Topics in Advanced Statistics/Assignment 2/")

source("Imputation.R")

# Put the code for your simulations below here

# Create missing value mechanisms

# MCAR
MCAR <- function(x, percentage_missing) {
  
  n = nrow(x)
  
  # Select the missing indices at random
  x[sample(seq(n), floor(n*percentage_missing)),1] = NA
  
  return(x)
  
}

# MAR
MAR <- function(x, percentage_missing) {
  
  # Select the missing indices based on a threshold value
  # quantile() makes sure that the number of missing values is equal to the group from which it is missing (thus equal probability)
  threshold = quantile(x[,2], 1-percentage_missing)
  x[x[,2] > threshold,1] = NA
  
  return(x)
  
}

# MNAR
MNAR <- function(x, percentage_missing) {
  
  threshold = quantile(x[,1], 1-percentage_missing)
  x[x[,1] > threshold,1] = NA
  
  return(x)
  
}

run_simulation_missing_values <- function(repetitions, percentage_missing, silence = TRUE) {
  
  set.seed(12345)
  
  p = 2
  n = 1000
  
  Sigma <- matrix(runif(p^2)*2-1, ncol=p)
  Sigma <- t(Sigma) %*% Sigma + diag(p)
  
  mcar_list_multiple = list()
  mar_list_multiple = list()
  mnar_list_multiple = list()
  
  mcar_list_bootstrap = list()
  mar_list_bootstrap = list()
  mnar_list_bootstrap = list()

  for (i in seq(repetitions)) {
    
    Betas = c(2,3)
    
    x = mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
    y = x %*% Betas + rnorm(n)
    
    mcar_x = MCAR(x, percentage_missing)
    mar_x  = MAR(x, percentage_missing)
    mnar_x = MNAR(x, percentage_missing)
    
    mcar_xy = cbind(mcar_x, y)
    mar_xy = cbind(mar_x, y)
    mnar_xy = cbind(mnar_x, y)
    
    mcar_list_multiple = append(mcar_list_multiple, list(pool(fit(multimp(mcar_xy, m = 20)))))
    mar_list_multiple = append(mar_list_multiple, list(pool(fit(multimp(mar_xy, m = 20)))))
    mnar_list_multiple = append(mnar_list_multiple, list(pool(fit(multimp(mnar_xy, m = 20)))))
    
    mcar_list_bootstrap = append(mcar_list_bootstrap, list(bootstrap(mcar_xy)$summary))
    mar_list_bootstrap = append(mar_list_bootstrap, list(bootstrap(mar_xy)$summary))
    mnar_list_bootstrap = append(mnar_list_bootstrap, list(bootstrap(mnar_xy)$summary))
    
    
    if (silence == FALSE) {
    
      if ((i %% (repetitions/10)) == 0) {
        
        print(paste("Progress: ", (i/repetitions)*100, '%', sep=""))
        
      } else if (i == repetitions) {
        
        print("Progress: 100%")
        print("Done!")
        
      }
    
    }
    
    
  }
  
  return_list = list("MCAR_multiple" = mcar_list_multiple,
                     "MCAR_bootstrap" = mcar_list_bootstrap,
                     "MAR_multiple" = mar_list_multiple,
                     "MAR_bootstrap" = mar_list_bootstrap,
                     "MNAR_multiple" = mnar_list_multiple,
                     "MNAR_bootstrap" = mnar_list_bootstrap)
  
  return(return_list)
  
}


run_simulation_cellwise_outliers <- function(repetitions, percentage_outlier, silence = TRUE) {
  
  set.seed(12345)
  
  p = 2
  n = 1000
  
  Sigma <- matrix(runif(p^2)*2-1, ncol=p)
  Sigma <- t(Sigma) %*% Sigma + diag(p)
  
  multimp_list = list()
  bootstrap_list= list()
  clean_list = list()
  
  for (i in seq(repetitions)) {
    
    Betas = c(2,3)
    
    x = mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
    y = x %*% Betas + rnorm(n)
    
    xy = cbind(x,y)
    
    colnames(xy ) = colnames(xy , do.NULL = FALSE)
    colnames(xy )[ncol(xy)] = "y"
    lmcontrol = lmrob.control(k.max = 1000, maxit.scale = 1000, na.action = NULL)
    
    clean_list = list()
    
    clean_matrix= matrix(0, nrow = p+1, ncol = 2)
    
    clean_result = lmrob(y ~ ., data = data.frame(xy), control = lmcontrol)
    
    if(!any(is.na(diag(clean_result$cov)))) {
      
      clean_matrix[,1] = clean_result$coefficients
      clean_matrix[,2] = sqrt(diag(clean_result$cov))
      
      clean_list = append(clean_list, list(clean_matrix))
      
    }
    
    clean_list = append(clean_list, list(clean_matrix))
    
    vectorized_x = c(x)
    number_of_cells = floor((n*p)*percentage_outlier)
    cellwise_outlier_indices = sample(seq(n*p), number_of_cells)
    vectorized_x[cellwise_outlier_indices] = rnorm(number_of_cells, mean = 10, sd = 1)
    
    cellwise_x = matrix(vectorized_x, nrow = n, ncol = p)
    
    cellwise_xy = cbind(cellwise_x, y)
    
    multimp_list = append(multimp_list, list(pool(fit(multimp(cellwise_xy, DDC = TRUE)))))
    bootstrap_list = append(bootstrap_list, list(bootstrap(cellwise_xy, DDC = TRUE)$summary))
    
    
    colnames(cellwise_xy) = colnames(cellwise_xy, do.NULL = FALSE)
    colnames(cellwise_xy)[ncol(cellwise_xy)] = "y"
    
    dirty_list = list()
    
    dirty_matrix= matrix(0, nrow = p+1, ncol = 2)
    
    dirty_result = lmrob(y ~ ., data = data.frame(cellwise_xy), control = lmcontrol)
    
    if(!any(is.na(diag(dirty_result$cov)))) {
      
      dirty_matrix[,1] = dirty_result$coefficients
      dirty_matrix[,2] = sqrt(diag(dirty_result$cov))
      
      dirty_list = append(dirty_list, list(dirty_matrix))
      
    }
    
    dirty_list = append(dirty_list, list(dirty_matrix))
      
    
    if (silence == FALSE) {
      
      if ((i %% (repetitions/10)) == 0) {
        
        print(paste("Progress: ", (i/repetitions)*100, '%', sep=""))
        
      } else if (i == repetitions) {
        
        print("Progress: 100%")
        print("Done!")
        
      }
      
    }
    
  }
  
  return_list = list("multimp" = multimp_list,
                     "bootstrap" = bootstrap_list,
                     "dirty" = dirty_list,
                     "clean" = clean_list)
  
  return(return_list)
  
}


missing_values_result = run_simulation_missing_values(repetitions = 200, percentage_missing = 0.20, silence = FALSE)

#save(missing_values_result, file = "missing_value_2902.RData")

cellwise_result = list()

cellwise_outlier_list = c(0.40)

for (outlier in cellwise_outlier_list) {
  
  cellwise_result = append(cellwise_result, list(run_simulation_cellwise_outliers(repetitions = 200, outlier, silence = FALSE)))
  
}


run_simulation_cellwise_increase <- function(repetitions, outlier_extremity, silence = TRUE) {
  
  set.seed(12345)
  
  percentage_outlier = 0.30
  
  p = 2
  n = 1000
  
  Sigma <- matrix(runif(p^2)*2-1, ncol=p)
  Sigma <- t(Sigma) %*% Sigma + diag(p)
  
  multimp_list = list()
  bootstrap_list= list()
  dirty_list = list()
  accuracy_list = list()
  
  for (i in seq(repetitions)) {
    
    Betas = c(2,3)
    
    x = mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
    y = x %*% Betas + rnorm(n)
    
    xy = cbind(x,y)
    
    vectorized_x = c(x)
    number_of_cells = floor((n*p)*percentage_outlier)
    cellwise_outlier_indices = sample(seq(n*p), number_of_cells)
    vectorized_x[cellwise_outlier_indices] = rnorm(number_of_cells, mean = 1, sd = 1)*outlier_extremity
    
    cellwise_x = matrix(vectorized_x, nrow = n, ncol = p)
    
    cellwise_xy = cbind(cellwise_x, y)
    
    multimp_list = append(multimp_list, list(pool(fit(multimp(cellwise_xy, DDC = TRUE)))))
    #bootstrap_list = append(bootstrap_list, list(bootstrap(cellwise_xy, DDC = TRUE)$summary))
    
    dirty_matrix= matrix(0, nrow = p+1, ncol = 2)
    
    prediction_cw_outlier = as.numeric(seq(n*p) %in% cellWise::DDC(cellwise_xy, DDCpars = list("silent" = TRUE))$indall)
    true_cw_outlier = as.numeric(seq(n*p) %in% cellwise_outlier_indices)
    
    true_negative = sum(as.numeric((true_cw_outlier == 0) & (prediction_cw_outlier == 0)))
    true_positive = sum(as.numeric((true_cw_outlier == 1) & (prediction_cw_outlier == 1)))
    false_positive = sum(as.numeric((true_cw_outlier == 0) & (prediction_cw_outlier == 1)))
    false_negative = sum(as.numeric((true_cw_outlier == 1) & (prediction_cw_outlier == 0)))
    
    specificity = true_positive/(false_negative+true_positive)
    sensitivity = true_negative/(false_positive+true_negative)
    accuracy = (true_negative+true_positive)/(n*p)
    
    accuracy_list = append(accuracy_list , list(cbind(specificity,sensitivity,accuracy)))
    
    colnames(cellwise_xy) = colnames(cellwise_xy, do.NULL = FALSE)
    colnames(cellwise_xy)[ncol(cellwise_xy)] = "y"
    
    lmcontrol = lmrob.control(k.max = 1000, maxit.scale = 1000, na.action = NULL)
    dirty_result = lmrob(y ~ ., data = data.frame(cellwise_xy), control = lmcontrol)
    
    if(!any(is.na(dirty_result$cov))) {
      
      dirty_matrix[,1] = dirty_result$coefficients
      dirty_matrix[,2] = sqrt(diag(dirty_result$cov))
      
      dirty_list = append(dirty_list, list(dirty_matrix))
      
    }
    
    dirty_list = append(dirty_list, list(dirty_matrix))
    
    if (silence == FALSE) {
      
      if ((i %% (repetitions/10)) == 0) {
        
        print(paste("Progress: ", (i/repetitions)*100, '%', sep=""))
        
      } else if (i == repetitions) {
        
        print("Progress: 100%")
        print("Done!")
        
      }
      
    }
    
  }
  
  return_list = list("multimp" = multimp_list,
                     "bootstrap" = bootstrap_list,
                     "dirty_list" = dirty_list,
                     "accuracy_list" = accuracy_list)
  
  return(return_list)
  
}


cellwise_result_extremity = list()

cellwise_outlier_list = c(1, 2.5, 5)

for (outlier in cellwise_outlier_list) {
  
  cellwise_result_extremity = append(cellwise_result_extremity, list(run_simulation_cellwise_increase(repetitions = 50, outlier, silence = FALSE)))
  
}

#save(cellwise_result_extremity, file = "cellwise_outlier_multimp.RData")
