### Generates a prediction interval based on assumed-normal variation about
### the grand outcome estimate of the mean, where the standard deviation
### is computed via our method-of-moments based estimator.

### `constant_helper()` is used when computing the moment based 
### estimator of gamma squared.


constant_helper <- function(s, K, var_cov_matrix, 
                            w){
  
  non_s_covs <- var_cov_matrix[s, -s]
  
  term_1 <- -2*(1-w[s])*sum(w[-s]*non_s_covs)
  
  term_2 <- 0
  
  for (i in setdiff(1:K, s)){
    
    summand <- w[i]*sum(w[-c(i, s)]*var_cov_matrix[i, -c(i, s)])
    
    term_2 <- term_2+summand
    
  } 
  
  C_s <- term_1 + term_2
  
  return(C_s)
}

get_interval_gamma <- function(target_data, 
                               trial_data, 
                               number_studies){
  
  ### First, we need to compute our estimate of gamma squared, which requires, 
  ### among other things, a bootstrapped distribution of the transported 
  ### potential outcome estimates.
  
  n_boot <- 1000
  
  trt_trial_data <- lapply(trial_data, 
                           function(x){x[assignment == 1]})
  
  transported_means <- lapply(trt_trial_data, 
                              function(x){make_prediction(target_data = target_data, 
                                                          trial_data = list(x))}) %>%
    unlist() %>%
    as.numeric()
  
  bootstrap_data <- matrix(0, nrow = n_boot, ncol = number_studies)
  
  for (b in 1:n_boot){
    
    ### Taking a bootstrap sample from the target population
    
    target_boot <- target_data[sample(1:nrow(target_data), 
                                      size = nrow(target_data), replace = TRUE)]
    
    ### Taking bootstrap samples of each of the observed studies
    
    study_boots <- lapply(trt_trial_data, 
                          function(x){x[sample(1:nrow(x), 
                                               size = nrow(x), 
                                               replace = TRUE)]})
    
    ### Computing transported outcome estimates based on those bootstrapped
    ### datasets.
    
    boot_transported_means <- lapply(study_boots, 
                                     function(x){make_prediction(target_data = target_boot, 
                                                                 trial_data = list(x))}) %>%
      unlist() %>%
      as.numeric()
    
    for(s in 1:number_studies){bootstrap_data[b, s] <- boot_transported_means[s]}
    
  }
  
  ### Great! Now we have our bootstrap samples. 
  
  
  # Computing our gamma squared ------------------------------------------
  
  ### Generating an estimated variance-covariance matrix
  
  var_cov_matrix <- cov(bootstrap_data)
  
  ### Renaming some things so the naming conventions in this script
  ### match that of the manuscript.
  
  y <- transported_means
  
  w <- rep(1, number_studies)
  
  w <- w/sum(w)
  
  y_w_sum <- sum(w * y)
  
  Q <- sum((y - y_w_sum) ^ 2)
  
  K <- number_studies
  
  study_specific_sds <- sqrt(diag(var_cov_matrix))
  
  sum_top   <- sum( ((w^2) * K + (1 - 2 * w)) * study_specific_sds ^ 2 )
  sum_denom <- sum( ((w^2) * K + (1 - 2 * w)) )
  
  ### Computing the somewhat complicated term arising from the correlation 
  ### using the helper function defined above
  
  C_s <- sapply(1:K, constant_helper, K = K, 
                var_cov_matrix = var_cov_matrix, 
                w = w, 
                simplify = TRUE)
  
  gamma_square_hat <- max(0,(Q - sum_top - sum(C_s)) / sum_denom)
  
  ### Getting the bootstrap estimate for the variance of the sample mean
  
  grand_mean_var <- bootstrap_data %>% 
    as.data.table() %>%
    .[, .(grand_mean = rowMeans(.SD))] %>%
    .[, var(grand_mean)]
  
  ### Ok! Now computing a prediction interval based on the 
  ### t-distribution with number_studies-2 degrees of freedom, per the 
  ### Higgins et al. paper on random effects meta-analysis
  
  t_quantile <- qt(p = .975, df = number_studies - 2, 
                   lower.tail = TRUE)
  
  pred_interval <- data.frame(lwr = mean(transported_means)-t_quantile*sqrt(grand_mean_var+gamma_square_hat), 
                              upper = mean(transported_means)+t_quantile*sqrt(grand_mean_var+gamma_square_hat))
  
  return(data.table(`Gamma Squared Interval` = pred_interval))
  
}