euler_method <- function(f, y0, t0, tn, h) {
  t <- seq(t0, tn, by = h)
  y <- numeric(length(t))
  y[1] <- y0  # Initial condition
  
  for (i in 1:(length(t) - 1)) {
    y[i + 1] <- y[i] + h * f(t[i], y[i])  # Euler update rule
  }
  
  return(data.frame(t = t, y = y))
}

# Midpoint rule for numerical integration
midpointrule <- function(f, lower, upper, step_size = 1e6) {
  h <- (upper - lower) / step_size  # Step size
  midpoints <- lower + h * (1:step_size) - h / 2
  integral_approx <- h * sum(f(midpoints))  # Midpoint rule formula
  
  return(integral_approx)
}

# Bisection method for root finding
bisection_method <- function(f, a, b, tol = 1e-6, max_iter = 1000) {
  if (f(a) * f(b) >= 0) {
    stop("Function values at a and b must have opposite signs")
  }
  
  iter <- 0
  c <- (a + b) / 2  # Midpoint
  
  while ((b - a) / 2 > tol && iter < max_iter) {
    c <- (a + b) / 2
    
    if (f(c) == 0) {
      break  # Exact root found
    } else if (f(a) * f(c) < 0) {
      b <- c  # Root is in the left half
    } else {
      a <- c  # Root is in the right half
    }
    
    iter <- iter + 1
  }
  
  return(list(root = c, iterations = iter, error = (b - a) / 2))
}

# Define the differential equation
dydx <- function(t, y) {
  return(y - t)
}

# Function to compute integral using the provided Midpoint Rule
compute_integral <- function(y_values, x_values) {
  f <- approxfun(x_values, y_values, method = "linear")  # Create function from data
  return(midpointrule(f, min(x_values), max(x_values), length(x_values)))
}

# Function to find correct y(0) using Bisection Method
find_correct_y0 <- function(target_integral, x0, xn, h) {
  y0_function <- function(y0) {
    solution <- euler_method(dydx, y0, x0, xn, h)
    integral_value <- compute_integral(solution$y, solution$t)
    return(integral_value - target_integral)
  }
  
  result <- bisection_method(y0_function, 0, 5)
  return(result$root)
}

# Solve for the correct y(0)
y0_correct <- find_correct_y0(2, 0, 1, 0.01)
solution <- euler_method(dydx, y0_correct, 0, 1, 0.01)

# Print the result
print(paste("Correct initial condition y(0):", y0_correct))
print(paste("Computed integral:", compute_integral(solution$y, solution$t)))

# Plot the solution
plot(solution$t, solution$y, type = "l", col = "blue", lwd = 2, 
     xlab = "x", ylab = "y(x)", main = "Euler's Method Solution")
