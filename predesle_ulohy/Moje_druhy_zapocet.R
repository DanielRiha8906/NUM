midpointrule <- function(f, lower, upper, step_size = 1e6) {
  h <- (upper - lower) / step_size
  midpoints <- lower + h * (1:step_size) - h / 2
  integral_approx <- h * sum(f(midpoints))
  return(integral_approx)
}

# Gaussian elimination for solving Ax = b
gauss_elimination <- function(A, b) {
  n <- nrow(A)
  for (k in 1:(n - 1)) {
    max_row <- which.max(abs(A[k:n, k])) + (k - 1)
    if (max_row != k) {
      A[c(k, max_row), ] <- A[c(max_row, k), ]
      b[c(k, max_row)] <- b[c(max_row, k)]
    }
    for (i in (k + 1):n) {
      factor <- A[i, k] / A[k, k]
      A[i, k:n] <- A[i, k:n] - factor * A[k, k:n]
      b[i] <- b[i] - factor * b[k]
    }
  }
  x <- numeric(n)
  for (i in n:1) {
    if (i < n) {
      x[i] <- (b[i] - sum(A[i, (i + 1):n] * x[(i + 1):n])) / A[i, i]
    } else {
      x[i] <- b[i] / A[i, i]
    }
  }
  return(x)
}

# Function to compute definite integral of t^(i+j-2) from 0 to 1
compute_Aij <- function(i, j) {
  f <- function(t) t^(i + j - 2)
  return(midpointrule(f, 0, 1, 1e6))  # Using midpoint rule for integration
}

# Function to construct matrix A
construct_A <- function(n) {
  A <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      A[i, j] <- compute_Aij(i, j)
    }
  }
  return(A)
}

# Function to compute vector b
compute_b <- function(A) {
  return(rowSums(A))
}

# Iterate over n = 2 to 10 and solve A * x = b
for (n in 2:10) {
  A <- construct_A(n)
  b <- compute_b(A)
  x <- gauss_elimination(A, b)
  
  print(paste("Solution for n =", n))
  print(x)
}
