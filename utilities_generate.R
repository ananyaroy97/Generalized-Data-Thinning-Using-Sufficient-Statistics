# Optimal value of epsilon for Data thinned samples

set.seed(123)
data("iris")
data <- iris[-5]
n <- nrow(data)
# Step 2: Simulate response variable
beta <- c(2, 1, 0.5, 1.5)  # True coefficients for covariates
error_variance <- 1  # Variance of the noise term
epsilon <- rnorm(n, mean = 0, sd = sqrt(error_variance))  # Error term
y <- as.matrix(data) %*% beta + epsilon  # Linear response variable with noise
data$response = y

# ===================================================================================================================================================================================
# Comparison of Regular and Data-Thinned Cross-Validation

# Step 1: Generate Synthetic Time-Series Data
set.seed(42)
n <- 100
x <- seq(0, 10, length.out = n)
y <- 5 * sin(x) + rnorm(n)

# Create a lagged dataset
data <- data.frame(y = y)
data$lag1 <- lag(data$y, 1)
data$lag2 <- lag(data$y, 2)
data <- na.omit(data)  # Remove rows with NAs

# Plot the synthetic time series
ggplot(data, aes(x = seq_along(y), y = y)) +
  geom_line() +
  labs(title = "Synthetic Time-Series Data", x = "Index", y = "y") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# ===================================================================================================================================================================================
# Changepoint Detection for Japan GDP Data

set.seed(2023)
data = read.csv("/Users/admin/Downloads/NCSU PhD/Semester I/ST793/Project/Japan GDP ChangePoint Detection/gdp.csv")
x = diff(as.vector(unlist(data[data$Country.Name == "Japan",][1, 5:62])))
noise = rnorm(length(x), 0, 1)
x = x + noise
sx = x^2

# ===================================================================================================================================================================================
# Changepoint Detection for LGA Airport International Passenger Data

set.seed(2023)
data = read.csv("air-passenger-traffic-per-month-port-authority-of-ny-nj-beginning-1977.csv")
x = diff(data$International.Passengers)
noise = rnorm(length(x), 0, 1)
x = x + noise
sx = x^2
