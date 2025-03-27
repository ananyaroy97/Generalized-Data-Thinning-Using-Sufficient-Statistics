# Optimal value of epsilon for Data thinned samples

## Required Packages
library(datathin)
## Main code for Analysis and Visualization
e = seq(0, 1, by = 0.01) # A sequence of epsilon values
m = NULL
for (i in e) {
  m = c(m, prediction_err(data, epsilon = i))
}
plot(e, m, main = "Prediction Error of Thinned Data for Different values of e", xlab = "e", ylab = "Mean Squared Error (MSE)")

# ===================================================================================================================================================================================
# Comparison of Regular and Data-Thinned Cross-Validation

## Required Packages
library(datathin)
library(caret)
library(ggplot2)
## Main code for Analysis and Visualization
# Step 2: Perform Multiple Iterations for Thinning and Regular Cross-Validation
n_iterations <- 1000
thinning_gap <- 2  # Number of subsets to thin into

# Containers for storing results
thinned_mse <- numeric(n_iterations)
regular_mse <- numeric(n_iterations)

set.seed(123)  # For reproducibility
for (i in 1:n_iterations) {
  # Thinning-based Cross-Validation
  thinning_result <- datathin(data$y, K = thinning_gap, epsilon = c(1/2, 1/2), family = "normal", arg = 1)  # Thin data
  thinned_datasets <- list(thinning_result[ , , 1], thinning_result[ , , 2])
  
  mse_scores <- numeric()  # For storing MSE scores from thinned subsets
  for (subset_idx in seq_along(thinned_datasets)) {
    subset_data <- data.frame(y = thinned_datasets[[subset_idx]])
    subset_data$lag1 <- lag(subset_data$y, 1)
    subset_data$lag2 <- lag(subset_data$y, 2)
    subset_data <- na.omit(subset_data)
    
    # Split into training and testing sets
    train_index <- 1:(0.8 * nrow(subset_data))
    train <- subset_data[train_index, ]
    test <- subset_data[-train_index, ]
    
    # Fit a regression model
    model <- lm(y ~ lag1 + lag2, data = train)
    predictions <- predict(model, newdata = test)
    
    # Calculate MSE
    mse_scores <- c(mse_scores, mean((test$y - predictions)^2))
  }
  thinned_mse[i] <- mean(mse_scores)
  
  # Regular Cross-Validation
  folds <- createFolds(data$y, k = 5)
  regular_mse_scores <- sapply(folds, function(indices) {
    train <- data[-indices, ]
    test <- data[indices, ]
    
    # Fit the model
    model <- lm(y ~ lag1 + lag2, data = train)
    predictions <- predict(model, newdata = test)
    mean((test$y - predictions)^2)
  })
  regular_mse[i] <- mean(regular_mse_scores)
}

# Step 3: Visualize the Results
# Combine results into a data frame
results <- data.frame(
  Method = c(rep("Thinned Cross-Validation", n_iterations), 
             rep("Regular Cross-Validation", n_iterations)),
  MSE = c(thinned_mse, regular_mse)
)

# Plot the results
ggplot(results, aes(x = Method, y = MSE)) +
  geom_boxplot(aes(fill = Method)) +
  labs(title = "Comparison of Cross-Validation Methods",
       x = "Method",
       y = "Mean Squared Error (MSE)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# Step 4: Print Summary of Results
cat("Thinned Cross-Validation Average MSE:", mean(thinned_mse), "\n")
cat("Regular Cross-Validation Average MSE:", mean(regular_mse), "\n")

# Violin plot for visualizing the distribution of MSE
ggplot(results, aes(x = Method, y = MSE, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "Comparison of MSE values under Regular CV and Thinned CV",
       x = "Method",
       y = "Mean Squared Error (MSE)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(legend.position = c(0.8, 0.8),
        plot.title = element_text(face = "bold", hjust = 0.5)) # Adjust position within the plot

# ===================================================================================================================================================================================
# Changepoint Detection for Japan GDP Data

## Required Packages
library(ggplot2)
library(patchwork)
library(dplyr)
library(extraDistr)
library(datathin)
## Main code for Analysis and Simulation
### Application

run_analysis(sx, "application")

### Type 1 error simulation

set.seed(1)
n = 57
rep = 1000
methods = c("naive", "thinned", "split")
result = list()

for (method in methods) {
  
  p_values = c()
  
  for (r in 1:rep) {
    
    x = rnorm(n)
    sx = x^2
    
    x_thinned = datathin(sx, family="gamma", arg=1/2)[,1,]
    x1 = x_thinned[, 1]
    x2 = x_thinned[, 2]
    
    if (method == "thinned") {
      cpt_data = x1
      test_data = x2
      minseglength = 10
    } else if (method == "naive") {
      cpt_data = sx
      test_data = sx
      minseglength = 10
    } else if (method == "split") {
      cpt_data = sx[seq(1, 57, 2)]
      test_data = sx[seq(2, 57, 2)]
      minseglength = 5
    }
    
    changepoints = estimate_changepoints(cpt_data, minseglen = minseglength)
    
    if (method == "split") {
      changepoints = changepoints - 1
      changepoints[1] = 0
      changepoints[length(changepoints)] = length(test_data)
    }
    
    if (length(changepoints) > 2) { 
      for (i in sample(1:(length(changepoints) - 2), 1)) {
        p_values = c(p_values, gamma_glm(test_data, changepoints[i], changepoints[i + 1], changepoints[i + 2]))
      }
    }
    
  }
  
  result[[method]] = p_values
  
}

pdf("simulation_T1E.pdf", height = 3, width = 8)
par(mfrow = c(1, 3))
plot(sort(result$naive), qunif(ppoints(length(result$naive))), 
     main = bquote("Naive"~""), xlab = "Significance Level", ylab = "Rejection Probability")
abline(0, 1, lty = "dashed")
plot(sort(result$split), qunif(ppoints(length(result$split))), 
     main = bquote("Order-preserved sample splitting"~""), xlab = "Significance Level", ylab = "Rejection Probability")
abline(0, 1, lty = "dashed")
plot(sort(result$thinned), qunif(ppoints(length(result$thinned))), 
     main = bquote("Generalized data thinning"~""), xlab = "Significance Level", ylab = "Rejection Probability")
abline(0, 1, lty = "dashed")
dev.off()

### Visualization of simulation under null and alternative

for (type in c("null", "alt")) {
  
  if (type == "null") {
    x = rnorm(n)
  } else {
    x = c(rnorm(20, sd = 2), rnorm(17, sd = 5), rnorm(20, sd = 1))
  }
  
  sx = x^2
  
  run_analysis(sx, paste0("simulation_", type), stability = FALSE)
  
}

# ===================================================================================================================================================================================
# Changepoint Detection for LGA Airport International Passenger Data

## Required Packages
library(ggplot2)
library(patchwork)
library(dplyr)
library(extraDistr)
library(datathin)
## Main code for Analysis and Simulation
run_analysis(sx, "application")

### Type 1 error simulation

set.seed(1)
n = 57
rep = 1000
methods = c("naive", "thinned", "split")
result = list()

for (method in methods) {
  
  p_values = c()
  
  for (r in 1:rep) {
    
    x = rnorm(n)
    sx = x^2
    
    x_thinned = datathin(sx, family="gamma", arg=1/2)[,1,]
    x1 = x_thinned[, 1]
    x2 = x_thinned[, 2]
    
    if (method == "thinned") {
      cpt_data = x1
      test_data = x2
      minseglength = 10
    } else if (method == "naive") {
      cpt_data = sx
      test_data = sx
      minseglength = 10
    } else if (method == "split") {
      cpt_data = sx[seq(1, 57, 2)]
      test_data = sx[seq(2, 57, 2)]
      minseglength = 5
    }
    
    changepoints = estimate_changepoints(cpt_data, minseglen = minseglength)
    
    if (method == "split") {
      changepoints = changepoints - 1
      changepoints[1] = 0
      changepoints[length(changepoints)] = length(test_data)
    }
    
    if (length(changepoints) > 2) { 
      for (i in sample(1:(length(changepoints) - 2), 1)) {
        p_values = c(p_values, gamma_glm(test_data, changepoints[i], changepoints[i + 1], changepoints[i + 2]))
      }
    }
    
  }
  
  result[[method]] = p_values
  
}

pdf("simulation_T1E.pdf", height = 3, width = 8)
par(mfrow = c(1, 3))
plot(sort(result$naive), qunif(ppoints(length(result$naive))), 
     main = bquote("Naive"~""), xlab = "Significance Level", ylab = "Rejection Probability")
abline(0, 1, lty = "dashed")
plot(sort(result$split), qunif(ppoints(length(result$split))), 
     main = bquote("Order-preserved sample splitting"~""), xlab = "Significance Level", ylab = "Rejection Probability")
abline(0, 1, lty = "dashed")
plot(sort(result$thinned), qunif(ppoints(length(result$thinned))), 
     main = bquote("Generalized data thinning"~""), xlab = "Significance Level", ylab = "Rejection Probability")
abline(0, 1, lty = "dashed")
dev.off()

### Visualization of simulation under null and alternative

for (type in c("null", "alt")) {
  
  if (type == "null") {
    x = rnorm(n)
  } else {
    x = c(rnorm(20, sd = 2), rnorm(17, sd = 5), rnorm(20, sd = 1))
  }
  
  sx = x^2
  
  run_analysis(sx, paste0("simulation_", type), stability = FALSE)
  
}