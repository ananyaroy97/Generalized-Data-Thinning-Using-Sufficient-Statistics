# Optimal value of epsilon for Data thinned samples

prediction_err <- function(data = data, epsilon = 0.5){
  y <- data$response
  # Step 2: Thin the data into two independent parts
  # epsilon <- 0.5  # Proportion for thinning
  y_scaled <- y - mean(y)
  thinned_data <- datathin(y_scaled, family = "normal", K = 2, epsilon = c(epsilon, 1 - epsilon), arg = 1)
  
  # Extract thinned parts
  y_train <- thinned_data[,1,][,1] + mean(y)
  y_test <- thinned_data[,1,][,2] + mean(y)
  X <- as.matrix(data[-5])
  # Split covariates into the same training and test sets
  X_train <- X
  X_test <- X
  
  # Step 4: Fit a linear regression model
  model <- lm(y_train ~ X_train - 1)  # No intercept since the model is centered
  
  # Step 5: Predict on the test set
  y_pred <- predict(model, newdata = data.frame(X_train = X_test))
  
  # Step 6: Calculate prediction error (Mean Squared Error)
  mse <- mean((y_test - y_pred)^2)
  
  # Output results
  cat("Model Coefficients:\n")
  print(coef(model))
  cat("\nMean Squared Error on Test Data:", mse, "\n")
  return(mse)
}

# ===================================================================================================================================================================================
# Changepoint Detection for Japan GDP Data

estimate_changepoints = function(x, pen.value = 10, minseglen = 10) {
  
  c(0, changepoint.np::cpt.np(x, minseglen = minseglen, penalty = "BIC")@cpts)
  
}

# testing for difference in means between (left+1):middle and (middle+1):right indices
gamma_glm = function(x, left, middle, right) {
  
  z = numeric(length(x))
  z[(middle + 1):right] = 1
  
  x = x[(left + 1):right]
  z = z[(left + 1):right]
  
  summary(glm(x ~ z, family = Gamma()))$coef[2, 4]
  
}

run_analysis = function(sx, filename, penalty = 10, minseglength = 10, stability = TRUE, zoom = NULL) {
  
  naive_p_values = c()
  naive_changepoints = estimate_changepoints(sx, pen.value = penalty, minseglen = minseglength)
  for (i in 1:(length(naive_changepoints) - 2)) {
    naive_p_values = c(naive_p_values, gamma_glm(sx, naive_changepoints[i], naive_changepoints[i + 1], naive_changepoints[i + 2]))
  }
  
  split_p_values = c()
  split_changepoints = estimate_changepoints(sx[seq(1, 57, 2)], pen.value = penalty, minseglen = minseglength / 2)
  split_changepoints_even = split_changepoints - 1
  split_changepoints_even[1] = 0
  split_changepoints_even[length(split_changepoints_even)] = length(sx) / 2
  for (i in 1:(length(split_changepoints) - 2)) {
    split_p_values = c(split_p_values, gamma_glm(sx[seq(2, 57, 2)], split_changepoints_even[i], split_changepoints_even[i + 1], split_changepoints_even[i + 2]))
  }
  split_changepoints = split_changepoints * 2 - 1
  split_changepoints[1] = 0
  split_changepoints[length(split_changepoints)] = length(sx)
  
  x_thinned = datathin(sx, family="gamma", arg=1/2)[,1,]
  x1 = x_thinned[, 1]
  x2 = x_thinned[, 2]
  
  thinned_p_values = c()
  thinned_changepoints = estimate_changepoints(x2, pen.value = penalty, minseglen = minseglength)
  for (i in 1:(length(thinned_changepoints) - 2)) {
    thinned_p_values = c(thinned_p_values, gamma_glm(x1, thinned_changepoints[i], thinned_changepoints[i + 1], thinned_changepoints[i + 2]))
  }
  
  naive_p_values <- p.adjust(naive_p_values, method="bonferroni")
  split_p_values <- p.adjust(split_p_values, method="bonferroni")
  thinned_p_values <- p.adjust(thinned_p_values, method="bonferroni")
  
  for (xlim in unique(list(NULL, zoom))) {
    
    if (!is.null(xlim)) {
      stability = FALSE
      filename = paste0(filename, "_zoomed")
    } else {
      xlim = c(0, length(x))
    }
    
    nrow = ifelse(stability, 5, 3)
    
    pdf(paste0(filename, ".pdf"), height = 1.5 * nrow, width = 8)
    par(mfcol = c(nrow, 1), mar = rep(2.5, 4))
    plot(x, type = "l", xlab = "", ylab = "", bty = "n", xlim = xlim, 
         main = bquote("Naive (Algorithm 2)":.(length(naive_changepoints) - 2)~"changepoints,"~.(sum(naive_p_values < 0.05))~"with p-value < 0.05 /"~.(length(naive_changepoints) - 2)))
    for (i in naive_changepoints[2:(length(naive_changepoints) - 1)]) {
      lines(x = c(i, i), y = c(min(x), max(x)), col = "red", lwd = 1)
    }
    for (i in naive_changepoints[2:(length(naive_changepoints) - 1)][which(naive_p_values < 0.05)]) {
      points(x = c(i, i), y = c(max(x), max(x)), col = "red", lwd = 3, pch = "*", cex = 2)
    }
    plot(x, type = "l", xlab = "", ylab = "", bty = "n", xlim = xlim, 
         main = bquote("Order-preserved sample splitting (Algorithm 3)":.(length(split_changepoints) - 2)~"changepoints,"~.(sum(split_p_values < 0.05))~"with p-value < 0.05 /"~.(length(split_changepoints) - 2)))
    for (i in split_changepoints[2:(length(split_changepoints) - 1)]) {
      lines(x = c(i, i), y = c(min(x), max(x)), col = "purple", lwd = 1)
    }
    for (i in split_changepoints[2:(length(split_changepoints) - 1)][which(split_p_values < 0.05)]) {
      points(x = c(i, i), y = c(max(x), max(x)), col = "purple", lwd = 3, pch = "*", cex = 2)
    }
    plot(x, type = "l", xlab = "", ylab = "", bty = "n", xlim = xlim, 
         main = bquote("Generalized data thinning (Algorithm 4)":.(length(thinned_changepoints) - 2)~"changepoints,"~.(sum(thinned_p_values < 0.05))~"with p-value < 0.05 /"~.(length(thinned_changepoints) - 2)))
    for (i in thinned_changepoints[2:(length(thinned_changepoints) - 1)]) {
      lines(x = c(i, i), y = c(min(x), max(x)), col = "blue", lwd = 1)
    }
    for (i in thinned_changepoints[2:(length(thinned_changepoints) - 1)][which(thinned_p_values < 0.05)]) {
      points(x = c(i, i), y = c(max(x), max(x)), col = "blue", lwd = 3, pch = "*", cex = 2)
    }
    
    repeated_thinned_changepoints = list()
    repeated_rejected_thinned_changepoints = list()
    for (rep in 1:100) {
      x_thinned = datathin(sx, family="gamma", arg=1/2)[,1,]
      x1 = x_thinned[, 1]
      x2 = x_thinned[, 2]
      cpts = estimate_changepoints(x2, pen.value = penalty, minseglen = minseglength)
      if (length(cpts) <= 2) {
        repeated_rejected_thinned_changepoints = repeated_rejected_thinned_changepoints
      } else {
        repeated_thinned_changepoints = c(repeated_thinned_changepoints, list(cpts[2:(length(cpts) - 1)]))
        pvals = c()
        for (i in 1:(length(cpts) - 2)) {
          pvals = c(pvals, gamma_glm(x1, cpts[i], cpts[i + 1], cpts[i + 2]))
        }
        pvals <- p.adjust(pvals, method="bonferroni")
        repeated_rejected_thinned_changepoints = c(repeated_rejected_thinned_changepoints, list(cpts[2:(length(cpts) - 1)][which(pvals < 0.05)]))
        
      }
    }
    if (stability) {
      hist(unlist(repeated_thinned_changepoints), breaks = length(x) / 10, 
           main = bquote("Percentage of data thinning replicates with changepoint"~""), xlab = NULL, freq = T, xlim = xlim)
      for (i in 1:5) {
        abline(v = i * 365, lty = "dashed")
      }
      hist(unlist(repeated_rejected_thinned_changepoints), breaks = length(x) / 10, 
           main = bquote("Percentage of data thinning replicates with changepoint p-value < 0.05 / (number of detected changepoints)"~""), xlab = NULL, freq = T, xlim = xlim)
      for (i in 1:5) {
        abline(v = i * 365, lty = "dashed")
      }
    }
    
    dev.off()
    
  }
  
}

# ===================================================================================================================================================================================
# Changepoint Detection for LGA Airport International Passenger Data

estimate_changepoints = function(x, pen.value = 10, minseglen = 10) {
  
  c(0, changepoint.np::cpt.np(x, minseglen = minseglen, penalty = "BIC")@cpts)
  
}

# testing for difference in means between (left+1):middle and (middle+1):right indices
gamma_glm = function(x, left, middle, right) {
  
  z = numeric(length(x))
  z[(middle + 1):right] = 1
  
  x = x[(left + 1):right]
  z = z[(left + 1):right]
  
  summary(glm(x ~ z, family = Gamma()))$coef[2, 4]
  
}

run_analysis = function(sx, filename, penalty = 10, minseglength = 10, stability = TRUE, zoom = NULL) {
  
  naive_p_values = c()
  naive_changepoints = estimate_changepoints(sx, pen.value = penalty, minseglen = minseglength)
  for (i in 1:(length(naive_changepoints) - 2)) {
    naive_p_values = c(naive_p_values, gamma_glm(sx, naive_changepoints[i], naive_changepoints[i + 1], naive_changepoints[i + 2]))
  }
  
  split_p_values = c()
  split_changepoints = estimate_changepoints(sx[seq(1, 1500, 2)], pen.value = penalty, minseglen = minseglength / 2)
  split_changepoints_even = split_changepoints - 1
  split_changepoints_even[1] = 0
  split_changepoints_even[length(split_changepoints_even)] = length(sx) / 2
  for (i in 1:(length(split_changepoints) - 2)) {
    split_p_values = c(split_p_values, gamma_glm(sx[seq(2, 1500, 2)], split_changepoints_even[i], split_changepoints_even[i + 1], split_changepoints_even[i + 2]))
  }
  split_changepoints = split_changepoints * 2 - 1
  split_changepoints[1] = 0
  split_changepoints[length(split_changepoints)] = length(sx)
  
  x_thinned = datathin(sx, family="gamma", arg=1/2)[,1,]
  x1 = x_thinned[, 1]
  x2 = x_thinned[, 2]
  
  thinned_p_values = c()
  thinned_changepoints = estimate_changepoints(x2, pen.value = penalty, minseglen = minseglength)
  for (i in 1:(length(thinned_changepoints) - 2)) {
    thinned_p_values = c(thinned_p_values, gamma_glm(x1, thinned_changepoints[i], thinned_changepoints[i + 1], thinned_changepoints[i + 2]))
  }
  
  naive_p_values <- p.adjust(naive_p_values, method="bonferroni")
  split_p_values <- p.adjust(split_p_values, method="bonferroni")
  thinned_p_values <- p.adjust(thinned_p_values, method="bonferroni")
  
  for (xlim in unique(list(NULL, zoom))) {
    
    if (!is.null(xlim)) {
      stability = FALSE
      filename = paste0(filename, "_zoomed")
    } else {
      xlim = c(0, length(x))
    }
    
    nrow = ifelse(stability, 5, 3)
    
    pdf(paste0(filename, ".pdf"), height = 1.5 * nrow, width = 8)
    par(mfcol = c(nrow, 1), mar = rep(2.5, 4))
    plot(x, type = "l", xlab = "", ylab = "", bty = "n", xlim = xlim, 
         main = bquote("Naive (Algorithm 2)":.(length(naive_changepoints) - 2)~"changepoints,"~.(sum(naive_p_values < 0.05))~"with p-value < 0.05 /"~.(length(naive_changepoints) - 2)))
    for (i in naive_changepoints[2:(length(naive_changepoints) - 1)]) {
      lines(x = c(i, i), y = c(min(x), max(x)), col = "red", lwd = 1)
    }
    for (i in naive_changepoints[2:(length(naive_changepoints) - 1)][which(naive_p_values < 0.05)]) {
      points(x = c(i, i), y = c(max(x), max(x)), col = "red", lwd = 3, pch = "*", cex = 2)
    }
    plot(x, type = "l", xlab = "", ylab = "", bty = "n", xlim = xlim, 
         main = bquote("Order-preserved sample splitting (Algorithm 3)":.(length(split_changepoints) - 2)~"changepoints,"~.(sum(split_p_values < 0.05))~"with p-value < 0.05 /"~.(length(split_changepoints) - 2)))
    for (i in split_changepoints[2:(length(split_changepoints) - 1)]) {
      lines(x = c(i, i), y = c(min(x), max(x)), col = "purple", lwd = 1)
    }
    for (i in split_changepoints[2:(length(split_changepoints) - 1)][which(split_p_values < 0.05)]) {
      points(x = c(i, i), y = c(max(x), max(x)), col = "purple", lwd = 3, pch = "*", cex = 2)
    }
    plot(x, type = "l", xlab = "", ylab = "", bty = "n", xlim = xlim, 
         main = bquote("Generalized data thinning (Algorithm 4)":.(length(thinned_changepoints) - 2)~"changepoints,"~.(sum(thinned_p_values < 0.05))~"with p-value < 0.05 /"~.(length(thinned_changepoints) - 2)))
    for (i in thinned_changepoints[2:(length(thinned_changepoints) - 1)]) {
      lines(x = c(i, i), y = c(min(x), max(x)), col = "blue", lwd = 1)
    }
    for (i in thinned_changepoints[2:(length(thinned_changepoints) - 1)][which(thinned_p_values < 0.05)]) {
      points(x = c(i, i), y = c(max(x), max(x)), col = "blue", lwd = 3, pch = "*", cex = 2)
    }
    
    repeated_thinned_changepoints = list()
    repeated_rejected_thinned_changepoints = list()
    for (rep in 1:100) {
      x_thinned = datathin(sx, family="gamma", arg=1/2)[,1,]
      x1 = x_thinned[, 1]
      x2 = x_thinned[, 2]
      cpts = estimate_changepoints(x2, pen.value = penalty, minseglen = minseglength)
      if (length(cpts) <= 2) {
        repeated_rejected_thinned_changepoints = repeated_rejected_thinned_changepoints
      } else {
        repeated_thinned_changepoints = c(repeated_thinned_changepoints, list(cpts[2:(length(cpts) - 1)]))
        pvals = c()
        for (i in 1:(length(cpts) - 2)) {
          pvals = c(pvals, gamma_glm(x1, cpts[i], cpts[i + 1], cpts[i + 2]))
        }
        pvals <- p.adjust(pvals, method="bonferroni")
        repeated_rejected_thinned_changepoints = c(repeated_rejected_thinned_changepoints, list(cpts[2:(length(cpts) - 1)][which(pvals < 0.05)]))
        
      }
    }
    if (stability) {
      hist(unlist(repeated_thinned_changepoints), breaks = length(x) / 10, 
           main = bquote("Percentage of data thinning replicates with changepoint"~""), xlab = NULL, freq = T, xlim = xlim)
      for (i in 1:5) {
        abline(v = i * 365, lty = "dashed")
      }
      hist(unlist(repeated_rejected_thinned_changepoints), breaks = length(x) / 10, 
           main = bquote("Percentage of data thinning replicates with changepoint p-value < 0.05 / (number of detected changepoints)"~""), xlab = NULL, freq = T, xlim = xlim)
      for (i in 1:5) {
        abline(v = i * 365, lty = "dashed")
      }
    }
    
    dev.off()
    
  }
  
}