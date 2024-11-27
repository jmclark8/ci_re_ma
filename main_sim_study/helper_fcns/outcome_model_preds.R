### Defines a simple function that generates predictions for an outcome
### model. Note that this does so by taking the mean of the average predictions
### across `trial_data`. Things still work if we want to get the transported
### effect from a single study: just make sure `trial_data` is one element
### in a list, i.e., `trial_data = list(one_study_data_frame)`

make_prediction <- function(target_data, 
                            trial_data){
  
  preds <- c()
  
  for (i in 1:length(trial_data)){
    
    one_study_data <- trial_data[[i]]
    
    one_study_model <- lm(y_a~X_1+X_2+X_3, data = one_study_data)
    
    target_preds <- predict(one_study_model, 
                            newdata = target_data[, .(X_1, X_2, X_3)])
    
    ### Adding the study-specific prediction to the vector of predictions
    
    preds <- c(preds, mean(target_preds))
    
  }
  
  return(mean(preds))
  
}