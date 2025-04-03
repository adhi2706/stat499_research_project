
install.packages("caret")
install.packages("dplyr")
install.packages("MCMCpack")
library(caret)
library(dplyr)
library(MCMCpack)

set.seed(123)

# Generating a random sample where each y is a one-hot encoded vector
generate_one_hot <- function(k, weights=NULL) {
  if (is.null(weights)) {
    weights <- rep(1/k, k) }
  
  vec <- rep(0, k)
  class_index <- sample(k, size = 1, prob = weights)
  vec[class_index] <- 1
  
  return(vec)
}

# Generate k random weights sum to 1
generate_weights <- function(k){
  x <- runif(k)
  return(x / sum(x))
}


# Function to calculate Bayes Conformity scores
# Returns vector of n+1 scores
bayes_cs <- function(y_candidate, Y_calib, K, alpha) {
  
  combined_data <- rbind(Y_calib, y_candidate)
  n <- nrow(Y_calib)
  
  scores <- numeric(n + 1)
  
  for (i in 1:(n + 1)) {
    leave_one_out <- combined_data[-i, ]
    observed_counts <- colSums(leave_one_out)
    category_i <- which(combined_data[i, ] == 1)
    
    scores[i] <- (alpha[category_i] + observed_counts[category_i]) / 
      (nrow(leave_one_out) + sum(alpha))
  }
  
  return(scores)
}

# Function to implement full cp using Bayes Conformity Score
bayes_full_cp <- function(Y_calib, K, alpha, error_prob) {
  n <- nrow(Y_calib)
  
  prediction_set <- matrix(0, nrow = 0, ncol = K)
  
  p_values <- numeric(K)
  
  for (i in 1:K) {
    
    y_candidate <- rep(0, K)
    y_candidate[i] <- 1
    
    scores <- bayes_cs(y_candidate, Y_calib, K, alpha)
    
    candidate_score <- scores[n + 1]  # c_{n+1}(y_{n+1})
    p_i <- sum(scores <= candidate_score) / (n + 1)
    p_values[i] <- p_i
    
    # If p_i > alpha, add to prediction set
    if (p_i > error_prob) {
      prediction_set <- rbind(prediction_set, y_candidate)
    }
  }
  
  return(list(calibration_set = Y_calib,
              prediction_set = prediction_set,
              p_values = p_values))
}


# Function to calculate Distance to Average (DTA) Conformity score
dta_cs <- function(y_candidate, Y_calib, K) {
  
  combined_data <- rbind(Y_calib, y_candidate)
  n <- nrow(Y_calib)
  
  scores <- numeric(n + 1)
  
  for (i in 1:(n + 1)) {
    
    leave_one_out <- combined_data[-i, ]
    y_bar <- colSums(leave_one_out) / n
    
    category_i <- which(combined_data[i, ] == 1)
    response_vec <- rep(0, K)
    response_vec[category_i] <- 1
    scores[i] <- norm(response_vec - y_bar, type = "2")
  }
  
  return(scores)
}


# Function to implement full cp using DTA Conformity Score |y_i - y_bar|
dta_full_cp <- function(Y_calib, K, error_prob) {
  
  n <- nrow(Y_calib)
  
  prediction_set <- matrix(0, nrow = 0, ncol = K)
  p_values <- numeric(K)
  
  for (i in 1:K) {
    
    y_candidate <- rep(0, K)
    y_candidate[i] <- 1
    
    scores <- dta_cs(y_candidate, Y_calib, K)
    
    candidate_score <- scores[n + 1]  # c_{n+1}(y_{n+1})
    p_i <- sum(scores >= candidate_score) / (n + 1)
    p_values[i] <- p_i
    
    # If p_i > alpha, add to prediction set
    if (p_i > error_prob) {
      prediction_set <- rbind(prediction_set, y_candidate)
    }
  }
  
  return(list(calibration_set = Y_calib,
              prediction_set = prediction_set,
              p_values = p_values))
}


validity_simulation <- function(category_lim, error_prob, n) {
  coverage_rates <- numeric()
  
  # Iterate through different numbers of categories
  for (K in 3:category_lim) {
    coverages <- numeric()
    
    # Repeat 2000 MC simulations for fixed K
    for (sim_rep in 1:2000) {
      
      # Generate alpha for Dirichlet prior
      alpha <- runif(K, min=0.1, max=10)
      
      rand_weights <- generate_weights(K)
      Y_cal <- t(sapply(1:n, function(i) generate_one_hot(K, rand_weights)))
      
      # Run Bayes full conformal prediction
      result <- bayes_full_cp(Y_cal, K, alpha, error_prob)
      #result <- dta_full_cp(Y_cal, K, error_prob)
      pred_set <- result$prediction_set
      
      # Calculate theoretical coverage based on weights
      pred_coverage <- 0
      for (i in 1:nrow(pred_set)) {
        category_num <- which(pred_set[i,] == 1)
        pred_coverage <- pred_coverage + rand_weights[category_num]
      }
      
      coverages[sim_rep] <- pred_coverage
    }
    
    coverage_rates[K-2] <- mean(coverages)
  }
  
  return(coverage_rates)
}


# Plot coverage rates for different sizes of calibration set
category_lim <- 20
epsilon <- 0.10
n_values <- c(20, 50, 100)
colors <- c("blue", "green", "black")

plot(NULL, xlim = c(3, category_lim), ylim = c(0.85, 1), 
     xlab = 'Number of Categories (K)', 
     ylab = 'Coverage Rate', 
     main = 'Coverage Rate of Full conformal prediction using Bayes Conformity Score')

abline(h = 1 - epsilon, col = 'red', lty = 2)

for (i in seq_along(n_values)) {
  n <- n_values[i]
  results <- validity_simulation(category_lim, error_prob = epsilon, n)
  lines(3:category_lim, results, type = 'b', col = colors[i], pch = 16)
}

legend("bottomright", legend = paste("n =", n_values), col = colors, lty = 1, pch = 16)


# Checking coverage for large K

n <- 20 # Number of observations in calibration set
epsilon <- 0.10 # User-specified error probability
category_lim <- 50

# Run the simulation
results <- validity_simulation(category_lim, error_prob = epsilon, n)

# Plot the results
plot(3:category_lim, results, 
     type = 'b', 
     col = "blue",
     xlab = 'Number of Categories (K)', 
     ylab = 'Coverage Rate', 
     main = 'Coverage Rate of Full conformal prediction using Bayes Conformity Score',
     xlim = c(20, category_lim),
     ylim = c(0.85, 1),
     pch = 16)
abline(h = 1 - epsilon, col = 'red', lty = 2)



##### Comparison of size of prediction set from different conformity scores #####

comparision_by_category <- function(category_lim, error_prob, n) {
  
  bayes_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 2000)
  dta_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 2000)
  
  for (K in 3:category_lim) {
    
    # Repeat 2000 MC simulations for fixed K
    for (sim_rep in 1:2000) {
      
      #alpha <- runif(K, min = 0.1, max = 1)
      alpha <- (rep(1,K))
      
      class_prob <- rdirichlet(1, alpha)[1,]
      Y <- t(rmultinom(n, 1, class_prob))
      
      result_bayes <- bayes_full_cp(Y, K, alpha, error_prob)
      result_dta <- dta_full_cp(Y, K, error_prob)
      
      bayes_prediction_sizes[K-2, sim_rep] <- nrow(result_bayes$prediction_set)
      dta_prediction_sizes[K-2, sim_rep] <- nrow(result_dta$prediction_set)
      
    }
  }
  
  return(list(bayes_sizes = apply(bayes_prediction_sizes, 1, mean),
              dta_sizes = apply(dta_prediction_sizes, 1, mean) ))
}



# Plot to visualize results of comparison by category count
plot_prediction_set_sizes <- function(results, category_lim) {
  
  plot(3:category_lim, results$bayes_sizes, 
       type = 'b', 
       col = 'blue', 
       xlab = 'Number of Categories (K)', 
       ylab = 'Average Prediction Set Size',
       main = 'Comparison of Prediction Set Sizes by Category Count',
       ylim = c(0, max(c(results$bayes_sizes, results$dta_sizes)) * 1.2))
  
  lines(3:category_lim, results$dta_sizes, 
        type = 'b', 
        col = 'red')
  
  legend('bottomright', 
         legend = c('Bayes Conformity Score', 'DTA Conformity Score'), 
         col = c('blue', 'red'), 
         lty = 1, 
         pch = 1)
  
  grid()
}

n <- 50 # Number of observations in calibration data
epsilon <- 0.05 # User-specified error probability
category_lim <- 20

# Run the simulation
compare_results <- comparision_by_category(category_lim, error_prob = epsilon, n)
plot_prediction_set_sizes(compare_results, category_lim)


# Function to produce a skewed Dirichlet concentration parameter vector
create_skewed_alpha <- function(K, skew_ratio) {
  # Initialize an empty vector for alpha
  alpha <- numeric(K)
  
  denominator <- (0.3 * skew_ratio + 0.7) * K
  # Set 30% of entries to higher value
  high_count <- round(0.3 * K)
  high_indices <- sample(1:K, high_count)
  
  # Assign values
  alpha[high_indices] <- skew_ratio/denominator
  alpha[-high_indices] <- 1/denominator
  #alpha[high_indices] <- skew_ratio
  #alpha[-high_indices] <- 1
  
  return(alpha)
}


# Comparing prediction sets when skewed prior assumption is true
compare_skewed_prior_true <- function(category_lim, error_prob, n, skew_ratio) {
  
  bayes_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 3000)
  dta_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 3000)
  
  for (K in 3:category_lim) {
    
    # Repeat 3000 times for fixed K
    for (sim_rep in 1:3000) {
      
      true_alpha <- create_skewed_alpha(K, skew_ratio)
      class_prob <- rdirichlet(1, true_alpha)[1,]
      Y <- t(rmultinom(n, 1, class_prob))
      
      assumed_alpha <- true_alpha
      result_bayes <- bayes_full_cp(Y, K, assumed_alpha, error_prob)
      result_dta <- dta_full_cp(Y, K, error_prob)
      
      bayes_prediction_sizes[K-2, sim_rep] <- nrow(result_bayes$prediction_set)
      dta_prediction_sizes[K-2, sim_rep] <- nrow(result_dta$prediction_set)
      
    }
  }
  
  return(list(bayes_sizes = apply(bayes_prediction_sizes, 1, mean),
              dta_sizes = apply(dta_prediction_sizes, 1, mean) ))
}


n <- 15 # Number of observations in calibration data
epsilon <- 0.1 # User-specified error probability
category_lim <- 20

# Run the simulation
compare_results <- compare_skewed_prior_true(category_lim, error_prob = epsilon, n, 4)
plot_prediction_set_sizes(compare_results, category_lim)


n <- 50 # Number of observations in calibration data
epsilon <- 0.1 # User-specified error probability
category_lim <- 20

# Run the simulation
compare_results <- compare_skewed_prior_true(category_lim, error_prob = epsilon, n, 4)
plot_prediction_set_sizes(compare_results, category_lim)



# Function to assigns higher alpha_j to categories different from true alpha used to generate data
create_incorrect_skew <- function(alpha) {
  unique_values <- unique(alpha)
  
  if(length(unique_values) <= 1) {
    return(alpha)
  }
  
  sorted_values <- sort(unique_values)
  high_value <- sorted_values[length(sorted_values)]
  low_value <- sorted_values[1]
  
  high_indices <- which(alpha == high_value)
  low_indices <- which(alpha == low_value)
  n_high <- length(high_indices)
  n_low <- length(low_indices)
  
  # Create a new alpha vector with shuffled categories
  new_alpha <- alpha
  new_high_indices <- sample(low_indices, n_high)
  new_alpha[] <- low_value
  new_alpha[new_high_indices] <- high_value
  
  # Check if the new assignment is the same as the original
  while(all(new_alpha == alpha)) {
    new_high_indices <- sample(low_indices, n_high)
    new_alpha[] <- low_value
    new_alpha[new_high_indices] <- high_value
  }
  
  return(new_alpha)
}


# Comparing prediction sets when skewed prior assumption is false
compare_skewed_prior_false <- function(category_lim, error_prob, n, skew_ratio) {
  
  bayes_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 3000)
  dta_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 3000)
  
  for (K in 3:category_lim) {
    
    # Repeat 3000 times for fixed K
    for (sim_rep in 1:3000) {
      
      true_alpha <- create_skewed_alpha(K, skew_ratio)
      #true_alpha <- (rep(0.1,K))
      class_prob <- rdirichlet(1, true_alpha)[1,]
      Y <- t(rmultinom(n, 1, class_prob))
      
      assumed_alpha <- create_incorrect_skew(true_alpha)
      #assumed_alpha <- create_skewed_alpha(K, skew_ratio)
      result_bayes <- bayes_full_cp(Y, K, assumed_alpha, error_prob)
      result_dta <- dta_full_cp(Y, K, error_prob)
      
      bayes_prediction_sizes[K-2, sim_rep] <- nrow(result_bayes$prediction_set)
      dta_prediction_sizes[K-2, sim_rep] <- nrow(result_dta$prediction_set)
      
    }
  }
  
  return(list(bayes_sizes = apply(bayes_prediction_sizes, 1, mean),
              dta_sizes = apply(dta_prediction_sizes, 1, mean) ))
}

n <- 15 # Number of observations in calibration data
epsilon <- 0.1 # User-specified error probability
category_lim <- 20

# Run the simulation
compare_results <- compare_skewed_prior_false(category_lim, error_prob = epsilon, n, 8)
plot_prediction_set_sizes(compare_results, category_lim)


n <- 50 # Number of observations in calibration data
epsilon <- 0.1 # User-specified error probability
category_lim <- 20

# Run the simulation
compare_results <- compare_skewed_prior_false(category_lim, error_prob = epsilon, n, 4)
plot_prediction_set_sizes(compare_results, category_lim)