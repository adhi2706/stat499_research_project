
install.packages("caret")
install.packages("dplyr")
install.packages("MCMCpack")
library(caret)
library(dplyr)
library(MCMCpack)

set.seed(NULL)

# Generating a random sample where each y is a one-hot encoded vector
generate_one_hot <- function(k) {
  y_cat <- sample(1:k, 1)
  y_vec <- rep(0, k)
  y_vec[y_cat] <- 1
  return(y_vec) }


# Function to calculate Bayes Conformity score
bayes_cs <- function(calibration_set_size, observed_counts, response_vec, alpha) {
  # Find the category number (position of 1 in the one-hot coded response vector)
  category_num <- which(response_vec == 1)
  
  # Calculate the Bayesian conformity score
  score <- (alpha[category_num] + observed_counts[category_num]) / (calibration_set_size + sum(alpha))
  
  return(score)
}


# Function to implement split cp using Bayes Conformity Score
bayes_split_cp <- function(Y, category_count, alpha, error_prob){
  
  # Splitting into 70% training and 30% calibration sets
  train_index <- sample(1:nrow(Y), size = ceiling(0.7 * nrow(Y)))
  training_set <- Y[train_index,]
  calibration_set <- Y[-train_index,]
  
  n <- nrow(calibration_set)
  observed_counts <- colSums(calibration_set)
  
  # Calculate conformity scores for each observation in the calibration set
  conformity_scores <- numeric(n)
  for (i in 1:n) {
    response_vec <- calibration_set[i,]
    conformity_scores[i] <- bayes_cs(calibration_set_size = n, observed_counts, 
                                     response_vec, alpha)
  }
  
  # Calculate q_hat
  q_hat_index <- max(1, floor((n + 1) * error_prob))
  q_hat <- sort(conformity_scores)[q_hat_index]
  
  
  # Generate prediction set
  prediction_set <- matrix(0, nrow = 0, ncol = category_count)
  
  # Iterate through each category
  for (i in 1:category_count) {
    # Generate response vector representing ith category
    response_vec <- rep(0, category_count)
    response_vec[i] <- 1
    
    # Calculate Bayes conformity score
    score <- bayes_cs(calibration_set_size = n, observed_counts, 
                      response_vec, alpha)
    
    # If score is greater than or equal to q_hat, add response_vec to prediction set
    if (score >= q_hat) {
      prediction_set <- rbind(prediction_set, response_vec)
    }
  }

    return(list(training_set = training_set,
                calibration_set = calibration_set,
                conformity_scores = conformity_scores,
                prediction_set = prediction_set,
                q_hat = q_hat))
}


# Function to implement split cp using trivial Conformity Score |y_i - y_bar|
trivial_split_cp <- function(Y, category_count, error_prob){
  
  # Splitting into 70% training and 30% calibration sets
  train_index <- sample(1:nrow(Y), size = ceiling(0.7 * nrow(Y)))
  training_set <- Y[train_index,]
  calibration_set <- Y[-train_index,]
  
  n <- nrow(calibration_set)
  y_bar <- colSums(calibration_set) / n
  
  # Calculate conformity scores for each observation in the calibration set
  conformity_scores <- numeric(n)
  for (i in 1:n) {
    response_vec <- calibration_set[i,]
    conformity_scores[i] <- norm(response_vec-y_bar , type = "2")
  }
  
  # Calculate q_hat
  q_hat_index <- min(n, ceiling((n + 1) * (1 - error_prob)))
  q_hat <- sort(conformity_scores)[q_hat_index]
  
  # Generate prediction set
  prediction_set <- matrix(0, nrow = 0, ncol = category_count)
  
  # Iterate through each category
  for (i in 1:category_count) {
    # Generate response vector representing ith category
    response_vec <- rep(0, category_count)
    response_vec[i] <- 1
    
    # Calculate conformity score as ||y_i - y_bar||
    score <- norm(response_vec-y_bar , type = "2")
    
    # If score is less than or equal to q_hat, add response_vec to prediction set
    if (score <= q_hat) {
      prediction_set <- rbind(prediction_set, response_vec)
    }
  }
  
  return(list(
    training_set = training_set,
    calibration_set = calibration_set,
    conformity_scores = conformity_scores,
    prediction_set = prediction_set,
    q_hat = q_hat
  ))
}

##### Validity Test #####

# Simulation to test Validity of Bayes-Optimal Conformal Prediction
validity_simulation <- function(category_lim, num_test_obs, error_prob, n) {
  
  coverage_rates <- numeric()
  
  # Iterate through different numbers of categories
  for (K in 3:category_lim) {

    coverages <- numeric()
    
    # Repeat multiple times to get stable estimate
    for (sim_rep in 1:200) {
      
      # (a) Choose random alpha vector with strictly positive components
      alpha <- runif(K, min = 0.1, max = 10)
      #alpha <- (rep(1,K))
      
      # (b) Generate random sample
      Y <- t(sapply(1:n, function(i) generate_one_hot(K)))
      
      # Run Bayes split conformal prediction
      result <- bayes_split_cp(Y, K, alpha, error_prob)
      #result <- trivial_split_cp(Y, K, alpha, error_prob)
      
      # (c) Generate test data
      test_data <- t(sapply(1:num_test_obs, function(i) generate_one_hot(K)))
      
      # Check coverage for test data
      coverage_count <- 0
      for (j in 1:nrow(test_data)) {
        # Check if the test data observation is in the prediction set
        coverage <- any(apply(result$prediction_set, 1, function(pred_vec) 
          all(pred_vec == test_data[j,])))
        coverage_count <- coverage_count + as.numeric(coverage)
      }
      
      # (d) Record proportion of prediction sets containing true value
      coverages[sim_rep] <- coverage_count / num_test_obs
    }
    
    # Store average coverage rate for this K
    coverage_rates[K-2] <- mean(coverages)
  }
  
  return(coverage_rates)
}

n <- 1000 # Number of observations in entire data set
epsilon <- 0.05 # User-specified error probability
category_lim <- 20
num_test_obs <- 500
# k <- 15 # Number of categories
# Y <- t(sapply(1:n, function(i) generate_one_hot(k)))

# Run the simulation
results <- validity_simulation(category_lim, num_test_obs, error_prob = epsilon, n)

# Plot the results
plot(3:category_lim, results, 
     type = 'b', 
     xlab = 'Number of Categories (K)', 
     ylab = 'Coverage Rate', 
     main = 'Coverage rate of split conformal prediction using Bayes Conformity Score',
     ylim = c(0.6, 1))
abline(h = 1 - epsilon, col = 'red', lty = 2)

##### Comparison of size of prediction set from different conformity scores #####

comparision_simulation <- function(category_lim, error_prob, n) {
  
  bayes_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 500)
  trivial_prediction_sizes <- matrix(0, nrow = category_lim - 2, ncol = 500)
  
  # Iterate through different numbers of categories
  for (K in 3:category_lim) {
    
    # Repeat multiple times to get stable estimate
    for (sim_rep in 1:500) {
      
      # (a) Choose random alpha vector with strictly positive components
      alpha <- runif(K, min = 0.1, max = 1)
      #alpha <- (rep(0.01,K))
      
      # (b) Generate random sample from multinomial dirichlet
      class_prob <- rdirichlet(1, alpha)[1,]
      Y <- t(rmultinom(n, 1, class_prob))
      
      # (c) Run both types of split conformal prediction
      result_bayes <- bayes_split_cp(Y, K, alpha, error_prob)
      result_trivial <- trivial_split_cp(Y, K, error_prob)
      
      # Store prediction set sizes
      bayes_prediction_sizes[K-2, sim_rep] <- nrow(result_bayes$prediction_set)
      trivial_prediction_sizes[K-2, sim_rep] <- nrow(result_trivial$prediction_set)
      
    }
  }
  
  return(list(bayes_sizes = apply(bayes_prediction_sizes, 1, mean),
              trivial_sizes = apply(trivial_prediction_sizes, 1, mean) ))
}



# Plot to visualize results of comparison by category count
plot_prediction_set_sizes <- function(results, category_lim) {
  
  plot(3:category_lim, results$bayes_sizes, 
       type = 'b', 
       col = 'blue', 
       xlab = 'Number of Categories', 
       ylab = 'Average Prediction Set Size',
       main = 'Comparison of Prediction Set Sizes by Category Count',
       ylim = c(0, max(c(results$bayes_sizes, results$trivial_sizes)) * 1.2))
  
  # Add trivial method results
  lines(3:category_lim, results$trivial_sizes, 
        type = 'b', 
        col = 'red')
  
  # Add legend
  legend('bottomright', 
         legend = c('Bayes Conformity Score', 'DTA Conformity Score'), 
         col = c('blue', 'red'), 
         lty = 1, 
         pch = 1)
  
  # Add grid
  grid()
}

n <- 200 # Number of observations in entire data set
epsilon <- 0.05 # User-specified error probability
category_lim <- 20

# Run the simulation
compare_results <- comparision_simulation(category_lim, error_prob = epsilon, n)
plot_prediction_set_sizes(compare_results, category_lim)

