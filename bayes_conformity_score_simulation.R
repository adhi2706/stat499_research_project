
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


# Function to calculate Bayes Conformity score
bayes_cs <- function(calibration_set_size, observed_counts, K, alpha) {
  
  scores <- numeric(K)
  
  # Calculate score for each possible category
  for (category_num in 1:K) {
    scores[category_num] <- (alpha[category_num] + observed_counts[category_num]) / 
      (calibration_set_size + sum(alpha))
  }
  return(scores)
}


# Function to implement split cp using Bayes Conformity Score
bayes_split_cp <- function(Y, K, alpha, error_prob){
  
  # Splitting into 70% training and 30% calibration sets
  #train_index <- sample(1:nrow(Y), size = ceiling(0.7 * nrow(Y)))
  #training_set <- Y[train_index,]
  calibration_set <- Y
  
  n <- nrow(calibration_set)
  observed_counts <- colSums(calibration_set)
  category_scores <- bayes_cs(calibration_set_size = n, observed_counts, K, alpha)
  
  # Calculate conformity scores for each observation in the calibration set
  conformity_scores <- numeric(n)
  for (i in 1:n) {
    category_num <- which(calibration_set[i,] == 1)
    conformity_scores[i] <- category_scores[category_num]
  }
  
  # Calculate q_hat
  q_hat_index <- max(1, floor((n + 1) * error_prob))
  q_hat <- sort(conformity_scores)[q_hat_index]
  
  
  # Generate prediction set
  prediction_set <- matrix(0, nrow = 0, ncol = K)
  
  # Iterate through each category
  for (i in 1:K) {
    # Generate response vector representing ith category
    response_vec <- rep(0, K)
    response_vec[i] <- 1
    
    # Calculate Bayes conformity score
    score <- category_scores[i]
    
    # If score is greater than or equal to q_hat, add response_vec to prediction set
    if (score >= q_hat) {
      prediction_set <- rbind(prediction_set, response_vec)
    }
  }
  
  return(list(calibration_set = calibration_set,
              conformity_scores = conformity_scores,
              prediction_set = prediction_set,
              q_hat = q_hat))
}


# Function to implement split cp using DTA Conformity Score |y_i - y_bar|
dta_split_cp <- function(Y, K, error_prob){
  
  # Splitting into 70% training and 30% calibration sets
  #train_index <- sample(1:nrow(Y), size = ceiling(0.7 * nrow(Y)))
  #training_set <- Y[train_index,]
  calibration_set <- Y
  
  n <- nrow(calibration_set)
  y_bar <- colSums(calibration_set) / n
  
  # Calculate conformity scores for each observation in the calibration set
  conformity_scores <- numeric(n)
  for (i in 1:n) {
    response_vec <- calibration_set[i,]
    conformity_scores[i] <- norm(response_vec - y_bar , type = "2")
  }
  
  # Calculate q_hat
  q_hat_index <- min(n, ceiling((n + 1) * (1 - error_prob)))
  q_hat <- sort(conformity_scores)[q_hat_index]
  #q_hat_quantile <- ceiling((n+1) * (1-error_prob)) / n
  #q_hat <- quantile(conformity_scores, probs = q_hat_quantile)
  
  # Generate prediction set
  prediction_set <- matrix(0, nrow = 0, ncol = K)
  
  # Iterate through each category
  for (i in 1:K) {
    # Generate response vector representing ith category
    response_vec <- rep(0, K)
    response_vec[i] <- 1
    
    # Calculate conformity score as ||y_i - y_bar||
    score <- norm(response_vec-y_bar , type = "2")
    
    # If score is less than or equal to q_hat, add response_vec to prediction set
    if (score <= q_hat) {
      prediction_set <- rbind(prediction_set, response_vec)
    }
  }
  
  return(list(calibration_set = calibration_set,
              conformity_scores = conformity_scores,
              prediction_set = prediction_set,
              q_hat = q_hat))
}


##### Validity Test #####

# Simulation to test Validity of Bayes-Optimal Conformal Prediction
validity_simulation <- function(category_lim, error_prob, n) {
  
  coverage_rates <- numeric()
  
  # Iterate through different numbers of categories
  for (K in 3:category_lim) {
    
    coverages <- numeric()
    
    # Repeat 200 times for fixed K
    for (sim_rep in 1:500) {
      
      alpha <- runif(K, min = 0.1, max = 10)
      #alpha <- (rep(1,K))
      rand_weights <- generate_weights(K)
      Y <- t(sapply(1:n, function(i) generate_one_hot(K, rand_weights)))
      
      # Run Bayes split conformal prediction
      result <- bayes_split_cp(Y, K, alpha, error_prob)
      #result <- dta_split_cp(Y, K, error_prob)
      
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

n <- 2000 # Number of observations in calibration set
epsilon <- 0.10 # User-specified error probability
category_lim <- 20
num_test_obs <- 500
# k <- 15 # Number of categories
# Y <- t(sapply(1:n, function(i) generate_one_hot(k)))

# Run the simulation
results <- validity_simulation(category_lim, error_prob = epsilon, n)

# Plot the results
plot(3:category_lim, results, 
     type = 'b', 
     xlab = 'Number of Categories (K)', 
     ylab = 'Coverage Rate', 
     main = 'Coverage rate of split conformal prediction using Bayes Conformity Score',
     ylim = c(0.8, 1))
abline(h = 1 - epsilon, col = 'red', lty = 2)


# Plot coverage rates for different sizes of calibration set
category_lim <- 20
epsilon <- 0.10
n_values <- c(500, 1000, 2000)  # Different values of n
colors <- c("blue", "green", "black")  # Colors for the lines

plot(NULL, xlim = c(3, category_lim), ylim = c(0.85, 1), 
     xlab = 'Number of Categories (K)', 
     ylab = 'Coverage Rate', 
     main = 'Coverage rate of split conformal prediction using Bayes Conformity Score')

abline(h = 1 - epsilon, col = 'red', lty = 2)  # Target coverage rate

for (i in seq_along(n_values)) {
  n <- n_values[i]
  results <- validity_simulation(category_lim, error_prob = epsilon, n)
  lines(3:category_lim, results, type = 'b', col = colors[i], pch = 16)
}

legend("bottomright", legend = paste("n =", n_values), col = colors, lty = 1, pch = 16)


##### Comparison of size of prediction set from different conformity scores #####

comparision_by_category <- function(category_lim, error_prob, n) {
  
  bayes_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 500)
  dta_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 500)
  
  for (K in 3:category_lim) {
    
    # Repeat 500 times for fixed K
    for (sim_rep in 1:500) {
      
      alpha <- runif(K, min = 0.1, max = 1)
      #alpha <- (rep(0.01,K))
      
      class_prob <- rdirichlet(1, alpha)[1,]
      Y <- t(rmultinom(n, 1, class_prob))
      
      result_bayes <- bayes_split_cp(Y, K, alpha, error_prob)
      result_dta <- dta_split_cp(Y, K, error_prob)
      
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
       xlab = 'Number of Categories', 
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

n <- 60 # Number of observations in calibration data
epsilon <- 0.05 # User-specified error probability
category_lim <- 20

# Run the simulation
compare_results <- comparision_by_category(category_lim, error_prob = epsilon, n)
plot_prediction_set_sizes(compare_results, category_lim)
