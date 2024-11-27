### Generates a prediction based on the wild bootstrap approach.

### For now, we're only computing the influence function for the outcome-model
### based estimate, i.e., not the DR estimate. 

upper_lower <- function(x, alpha  = 0.05){
  data.frame(lwr = quantile(x, probs = alpha/2, names = FALSE), 
             upper = quantile(x, probs = 1-alpha/2, names = FALSE))
}

sample_rademacher <- function(n){
  
  samples <- sample(x = c(-1, 1), 
                    size = n, 
                    prob = c(0.5, 0.5), 
                    replace = TRUE)
  
  return(samples)
  
  
}

get_interval_wild <- function(target_data, trial_data, number_studies
                              ### `trial data` is a LIST of data from the
                              ### various studies
){
  
  ### First things first: we need to estimate various quantities in 
  ### the influence function for expected potential outcome under treatment
  ### a in study s, transported to the target population
  
  ### Models relating outcome under treatment to covariates specific to each
  ### study
  
  study_outcome_models <- lapply(trial_data, 
                                 function(x){lm(y_a ~ X_1+X_2+X_3, 
                                                data = x[assignment == 1])})
  
  ### Getting predicted outcomes transported from each study to the target
  ### population. In practice, this involves simply making a study-specific
  ### prediction using each of the outcome models fit above on data from
  ### the target population
  
  est_transported_outcomes <- lapply(trial_data, 
                                     function(x){make_prediction(
                                       target_data = target_data, 
                                       trial_data = list(x[assignment == 1])
                                     )})
  
  ### Now we compute everyone's influence function.
  ### We start by combining the data.
  
  # Computing Influence Functions -------------------------------------------
  
  
  ### Note that we're copying data here so we don't make changes to the
  ### underlying data tables. We rbind the target data and the data from 
  ### each study.
  
  infl_func_data <- rbindlist(c(list(copy(target_data)[, `:=` (setting = 0, 
                                                               y_a = NA)]), 
                                mapply(function(data, index){
                                  copy(data)[, setting := index] %>%
                                    .[assignment == 1, 
                                      .(X_1, X_2, X_3, setting, y_a)]
                                }, 
                                data = trial_data, 
                                index = 1:length(trial_data), 
                                SIMPLIFY = FALSE)), 
                              use.names = TRUE)
  
  ### We compute the influence function separately for each setting-specific 
  ### outcome model. Note that `setting_index` should be a number in 1:n_studies.
  
  ### I'm defining this function within `get_interval_wild` so we can use
  ### the objects defined within `get_interval_wild` w/o passing them as
  ### arguments.
  
  compute_if_s <- function(setting_index){
    
    betas <- study_outcome_models[[setting_index]]$coefficients %>%
      as.matrix() %>%
      unname()
    
    if_s_data <- infl_func_data[setting == setting_index | setting == 0]
    
    ### First, getting some non person-specific things we'll need in the
    ### computation
    
    total_n <- Reduce(`+`, x = lapply(trial_data, nrow))+nrow(target_data)
    
    p_0 <- nrow(target_data)/(total_n)
    
    p_s_a <- trial_data[[setting_index]][assignment == 1, .N]/total_n
    
    avg_target_x <- if_s_data[setting == 0, lapply(.SD, mean), 
                              .SDcols = c("X_1", "X_2", "X_3")] %>%
      unlist() %>%
      matrix() %>%
      ### Adding the intercept
      rbind(1, .) %>%
      t()
    
    x_s <- if_s_data[setting == setting_index] %>%
      .[, .(X_1, X_2, X_3)] %>%
      as.matrix() %>%
      cbind(1, .)
    
    x_x_t_avg <- t(x_s) %*% x_s %>%
      {./if_s_data[setting == setting_index, .N]} %>%
      unname()
    
    ### Now defining a function that'll be used on each individual row
    ### of if_s_data, i.e., for each participant 
    
    compute_row_wise_if <- function(row_wise_data){
      
      setting <- row_wise_data[, setting]
      
      if(setting == 0){
        outcome_model_pred <- predict(study_outcome_models[[setting_index]], 
                                      newdata = row_wise_data[, .(X_1, X_2, X_3)]) %>%
          as.numeric()
        
        if_piece <- (1/p_0)*(outcome_model_pred - mean(unlist(est_transported_outcomes)))
        
        return(if_piece)
      }
      
      if(setting == setting_index){
        
        x_i <- row_wise_data[, .(X_1, X_2, X_3)] %>%
          as.matrix() %>%
          unname() %>%
          cbind(1, .) %>%
          t()
        
        y_i <- row_wise_data[, y_a]
        
        x_i_y_i_resid <- x_i %*% (y_i-t(x_i)%*%betas)
        
        if_piece <- (1/p_s_a) * avg_target_x %*% solve(x_x_t_avg) %*% x_i_y_i_resid
        
        return(if_piece)
      }
    }
    
    ### Computing the piece of the influence function for each row:
    
    if_s_data <- if_s_data[, paste0("if_", setting_index) := compute_row_wise_if(.SD), 
                           by = seq_len(nrow(if_s_data))]
    
    return(if_s_data)
    
  }
  
  ### Before computing all of the various influence functions, I'm adding
  ### a `study_id` for each participant that will be helpful when merging things together again
  
  infl_func_data[, study_id := .I]
  
  sep_infl_func <- lapply(1:number_studies, 
                          compute_if_s) 
  
  wild_bootstrap_data <- lapply(sep_infl_func, 
                                function(x){x[, !c("X_1", "X_2", "X_3", "y_a", 
                                                   "setting")]}) %>%
    Reduce(function(x, y){merge(x, y, by = "study_id", 
                                all = TRUE)}, .)
  
  ### Replacing NAs with zero for those rows where both of the indicators
  ### in a particular influence function weren't satisfied
  
  wild_bootstrap_data[is.na(wild_bootstrap_data)] <- 0
  
  
  # Performing the bootstrap ------------------------------------------------
  
  bootstrap_data <- matrix(0, nrow = 1000, ncol = 1)
  
  for (b in 1:1000){
    
    unobs_study <- sample(1:number_studies, size = 1)
    
    unobs_study_col <- paste0("if_", unobs_study)
    
    obs_study_cols <- paste0("if_", setdiff(1:number_studies, unobs_study))
    
    ### Computing the transported mean outcome for the study we don't observe, 
    ### and factoring in sampling variability of that quanity via the wild bootstrap
    
    unobs_study_resid <- wild_bootstrap_data[, mean(sample_rademacher(nrow(infl_func_data))*get(unobs_study_col))]
    
    boot_study_outcome <- est_transported_outcomes[[unobs_study]]-unobs_study_resid
    
    ### Computing the grand mean for the studies we do observe and 
    ### factoring in sampling variability of that quantity via the wild bootstrap
    
    obs_study_grand_mean <- est_transported_outcomes[-c(unobs_study)] %>%
      as.numeric() %>%
      mean()
    
    ### This involves taking the mean of n_studies-1 residuals like
    ### `unobs_study_resid`
    
    obs_study_mean_resid <- copy(wild_bootstrap_data) %>%
      .[, mean_infl_func := rowMeans(.SD), 
        .SDcols = obs_study_cols] %>%
      .[, mean(sample_rademacher(nrow(infl_func_data))*mean_infl_func)]
    
    boot_grand_mean <- obs_study_grand_mean-obs_study_mean_resid
    
    delta <- boot_grand_mean-boot_study_outcome
    
    prediction <- mean(unlist(est_transported_outcomes))-delta
    
    bootstrap_data[b, 1] <- prediction
    
    
  }
  
  ### Computing an interval with the quantiles of the bootstrap data.
  
  pred_interval_quantile <- upper_lower(x = bootstrap_data[, 1], 
                                        alpha = 0.05)
  
  ### Computing an interval using a normal approximation incorporating
  ### the standard error of the bootstrap distribution.
  
  est_se <- sd(bootstrap_data[, 1])
  
  normal_bootstrap_interval <- data.frame(`lwr` = mean(unlist(est_transported_outcomes)) - 1.96*est_se, 
                                          `upper` = mean(unlist(est_transported_outcomes)) + 1.96*est_se)
  
  
  return(data.table(`Wild Bootstrap Interval Quantiles` = pred_interval_quantile, 
                    `Wild Bootstrap Interval Normal Approx` = normal_bootstrap_interval))
  
}