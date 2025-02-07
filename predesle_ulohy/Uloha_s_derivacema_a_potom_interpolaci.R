# F(x) = cos(5 * arccos(x))
# Compute first & second derivatives
# Solve f(x) = 0
# Perform interpolation

# Function definition
f <- function(x) {
  return(cos(5 * acos(x)))
}

# Numerical first derivative using central difference method
numerical_derivative <- function(f, x, h = 1e-5) {
  return((f(x + h) - f(x - h)) / (2 * h))
}

# Numerical second derivative using central difference method
numerical_second_derivative <- function(f, x, h = 1e-5) {
  return((f(x + h) - 2 * f(x) + f(x - h)) / (h^2))
}

# Generate 10 points in [-1, 1]
x_values <- seq(-1, 1, length.out = 10)

# Compute derivatives
first_derivative <- sapply(x_values, function(x) numerical_derivative(f, x))
second_derivative <- sapply(x_values, function(x) numerical_second_derivative(f, x))

# Print computed derivatives
print(data.frame(x = x_values, first_derivative = first_derivative, second_derivative = second_derivative))

# Bisection Method for finding roots of f(x) = 0
bisection_method <- function(f, a, b, tol = 1e-5, max_iter = 100) {
  if (f(a) * f(b) > 0) {
    stop("No root in the given interval (function does not change sign)")
  }
  
  for (i in 1:max_iter) {
    c <- (a + b) / 2  # Midpoint
    fc <- f(c)
    
    if (abs(fc) < tol) {
      return(c)  # Root found
    }
    
    if (f(a) * fc < 0) {
      b <- c  # Root is in left half
    } else {
      a <- c  # Root is in right half
    }
  }
  
  return((a + b) / 2)  # Return best approximation
}

# Midpoint Rule for Numerical Integration
midpointrule <- function(f, dolni_hranice, vrchni_hranice, step_size = 1000) {
  h <- (vrchni_hranice - dolni_hranice) / step_size
  midpoints <- dolni_hranice + h * (1:step_size) - h / 2
  integral_approx <- h * sum(f(midpoints))
  return(integral_approx)
}

# Find zero-crossing intervals
x_fine <- seq(-1, 1, length.out = 1000)  # Fine grid for root detection
y_fine <- sapply(x_fine, f)

zero_crossing_intervals <- list()
for (i in 1:(length(x_fine) - 1)) {
  if (y_fine[i] * y_fine[i + 1] < 0) {  # Sign change detected
    zero_crossing_intervals <- append(zero_crossing_intervals, list(c(x_fine[i], x_fine[i + 1])))
  }
}

# Compute roots using Bisection Method
roots <- c()
for (interval in zero_crossing_intervals) {
  root <- bisection_method(f, interval[1], interval[2])  # Apply bisection
  roots <- c(roots, root)
}

# Print found roots
print("Roots of F(x) = 0:")
print(roots)

# Lagrange Interpolation Function
lagrange_interpolation <- function(x, y, xi) {
  n <- length(x)  # Number of known data points
  result <- 0     # Initialize the result variable
  
  for (i in 1:n) {
    li <- 1  # Initialize Lagrange basis polynomial L_i(x)
    
    for (j in 1:n) {
      if (i != j) {  # Avoid division by itself
        li <- li * (xi - x[j]) / (x[i] - x[j])
      }
    }
    
    result <- result + li * y[i]
  }
  
  return(result)  # Return interpolated value at xi
}

# Compute interpolation for finer points
x_interp <- seq(-1, 1, length.out = 50)
y_interp <- sapply(x_interp, function(x) lagrange_interpolation(x_values, f(x_values), x))

# Plot Function, Roots, Derivatives, and Interpolation
plot(x_fine, y_fine, type="l", col="blue", lwd=2, xlab="x", ylab="y",
     main="Function, Derivatives, Roots & Interpolation")

# Add interpolated function in red
lines(x_interp, y_interp, col="red", lwd=2, lty=2)

# Mark found roots
points(roots, rep(0, length(roots)), col="black", pch=19)

# Mark derivatives
points(x_values, first_derivative, col="green", pch=4)
points(x_values, second_derivative, col="purple", pch=3)

legend("topright", legend=c("F(x)", "Interpolation", "Roots", "First Derivative", "Second Derivative"),
       col=c("blue", "red", "black", "green", "purple"), lty=c(1, 2, NA, NA, NA), pch=c(NA, NA, 19, 4, 3), lwd=2)

#Takhle by to asi mělo fungovat?