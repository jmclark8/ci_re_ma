### Generates a prediction interval based on our simpler bootstrap approach. 
### Note that this generates 1000 bootstrap samples, whenever applicable.

upper_lower <- function(x, alpha  = 0.05){
  data.frame(lwr = quantile(x, probs = alpha/2, names = FALSE), 
             upper = quantile(x, probs = 1-alpha/2, names = FALSE))
}

get_interval_simple <- function(target_data, trial_data, number_studies
                                ### `trial data` is a LIST of data from the
                                ### various studies
){
  ### For now, hard-coding this at 1000; in the future, we might change this
  ### and make it a parameter for `get_interval_simple()`
  
  n_boot <- 1000
  
  ### First, only retaining the trial data from participants assigned to 
  ### treatment
  
  trt_trial_data <- lapply(trial_data, 
                           function(x){x[assignment == 1]})
  
  est_transported_outcomes <- lapply(trial_data, 
                                     function(x){make_prediction(
                                       target_data = target_data, 
                                       trial_data = list(x[assignment == 1])
                                     )})
  
  bootstrap_data <- matrix(0, nrow = n_boot, ncol = 2)
  
  for (b in 1:n_boot){
    
    ### Step 1: Choosing a random study to be "unobserved"
    
    unobs_study <- sample(1:number_studies, size = 1)
    
    unobs_study_data <- trt_trial_data[[unobs_study]]
    
    obs_study_data <- trt_trial_data[setdiff(1:n_studies, 
                                             unobs_study)]
    
    ### Step 2: Transporting the effect from that study to our target pop
    ###         This is the "true"value of the prediction we're making
    
    unobs_study_boot <- unobs_study_data[sample(1:nrow(unobs_study_data), 
                                                size = nrow(unobs_study_data), replace = TRUE)]
    
    target_boot_1 <- target_data[sample(1:nrow(target_data), 
                                        size = nrow(target_data), replace = TRUE)]
    
    true_value <- make_prediction(target_data = target_boot_1, 
                                  trial_data = list(unobs_study_boot))
    
    ### Step 3: Drawing bootstrap samples for all of the observed study data,
    ### drawing another bootstrap sample from the target population, and 
    ### then applying our method to those data in order to estimate `true_value`
    
    target_boot_2 <- target_data[sample(1:nrow(target_data), 
                                        size = nrow(target_data), replace = TRUE)]
    
    obs_study_boot <- lapply(obs_study_data, 
                             function(x){x[sample(1:nrow(x), 
                                                  size = nrow(x), replace = TRUE)]})
    
    true_value_pred <- make_prediction(target_data = target_boot_2, 
                                       trial_data = obs_study_boot)
    
    ### Step 4: Compute and save the bootstrap residual between the `true_value`
    ### and `true_value_pred`
    
    delta_boot <- true_value_pred-true_value
    
    ### Step 5: Computing our estimator using a bootstrapped sample from 
    ### ALL the data
    
    all_study_boot <- lapply(trt_trial_data, 
                             function(x){x[sample(1:nrow(x), 
                                                  size = nrow(x), replace = TRUE)]})
    
    target_boot_3 <- target_pop_data[sample(1:nrow(target_data), 
                                            size = nrow(target_data), replace = TRUE)]
    
    all_data_pred <- make_prediction(target_data = target_boot_3, 
                                     trial_data = all_study_boot)
    
    pred_boot <- all_data_pred - delta_boot
    
    bootstrap_data[b, 1] <- delta_boot
    bootstrap_data[b, 2] <- all_data_pred-delta_boot
    
  }
  
  ### Constructing a prediction interval using the bootstrapped values. 
  ### As detailed in the manuscript, this includes both a normal apporximation
  ### using the SE of the bootstrap estimates and the quantiles from the bootstrap
  ### distribution itself.
  
  pred_interval_quantile <- upper_lower(x = bootstrap_data[, 2], 
                                        alpha = 0.05)
  
  est_se <- sd(bootstrap_data[, 2])
  
  normal_bootstrap_interval <- data.frame(`lwr` = mean(unlist(est_transported_outcomes)) - 1.96*est_se, 
                                          `upper` = mean(unlist(est_transported_outcomes)) + 1.96*est_se)
  
  
  return(data.table(`Simple Bootstrap Interval Quantiles` = pred_interval_quantile, 
                    `Simple Bootstrap Interval Normal Approx` = normal_bootstrap_interval))
  
}