### The manuscript references a small change to our gamma-squared estimator
### necessary to carry out the data analysis. This script includes a code
### snippet that illustrates that change. 

### The quantity `gamma_square_hat_boot` computed in this snipped would
### take the place of `gamma_square_hat` defined in 
### `main_sim_study/helper_fcns/gamma_squared_interval.R`. It would follow
### the definition of `C_s` in that script. Finally, `boot_cols` refers to
### a character vector of the column names in `bootstrap_data()`.

bootstrap_q_data <- data.table(bootstrap_data) %>%
  ### Computing Q for each of the bootstrap samples
  .[, Q_boot := rowSums((.SD-rowMeans(.SD))^2), 
    .SDcols = boot_cols] %>%
  .[, q := (Q_boot - sum_top - sum(C_s)) / sum_denom] %>%
  .[, gamma_square_boot := ifelse(q <= 0, 0, q)]

gamma_square_hat_boot <- bootstrap_q_data[, mean(gamma_square_boot)]