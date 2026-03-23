# ============================================================
# Optimal split boundary using IG and D_KL
# Job 1: Labeled data (X, Y) — supervised optimal cut
# Job 2: Unlabeled data X   — unsupervised optimal cut via GMM
# ============================================================

# ── Shared helper functions ──

H <- function(y) {
  p <- table(y) / length(y)
  -sum(p * log2(p))
}

DKL <- function(p, q) {
  e <- 1e-10
  p <- p / sum(p) + e
  q <- q / sum(q) + e
  sum(p * log(p / q))
}

best_split <- function(x, y, method = c("IG", "KL")) {
  method <- match.arg(method)
  idx <- order(x); xs <- x[idx]; ys <- y[idx]; n <- length(x)
  cls <- sort(unique(y))
  p0 <- table(factor(y, levels = cls)) / n
  cuts <- which(diff(xs) != 0)
  if (length(cuts) == 0) return(list(boundary = NA, score = 0))
  scores <- sapply(cuts, function(i) {
    L <- ys[1:i]; R <- ys[(i + 1):n]
    if (method == "IG") {
      H(y) - (i / n) * H(L) - ((n - i) / n) * H(R)
    } else {
      pL <- table(factor(L, levels = cls)) / length(L)
      pR <- table(factor(R, levels = cls)) / length(R)
      (i / n) * DKL(pL, p0) + ((n - i) / n) * DKL(pR, p0)
    }
  })
  best <- which.max(scores)
  list(boundary = (xs[cuts[best]] + xs[cuts[best] + 1]) / 2, score = scores[best])
}

# ── EM for 2-component Gaussian Mixture Model ──

gmm_em <- function(x, max_iter = 200, tol = 1e-6) {
  n <- length(x)
  # Initialize: split at median
  mu <- c(mean(x[x <= median(x)]), mean(x[x > median(x)]))
  sigma <- c(sd(x), sd(x))
  pi_k <- c(0.5, 0.5)

  for (iter in 1:max_iter) {
    # E-step: compute responsibilities
    d1 <- pi_k[1] * dnorm(x, mu[1], sigma[1])
    d2 <- pi_k[2] * dnorm(x, mu[2], sigma[2])
    total <- d1 + d2
    gamma <- d1 / total  # responsibility for component 1

    # M-step: update parameters
    n1 <- sum(gamma); n2 <- n - n1
    mu_new <- c(sum(gamma * x) / n1, sum((1 - gamma) * x) / n2)
    sigma_new <- c(
      sqrt(sum(gamma * (x - mu_new[1])^2) / n1),
      sqrt(sum((1 - gamma) * (x - mu_new[2])^2) / n2)
    )
    pi_new <- c(n1 / n, n2 / n)

    # Check convergence
    if (max(abs(mu_new - mu)) < tol) break
    mu <- mu_new; sigma <- sigma_new; pi_k <- pi_new
  }

  # Ensure component 1 has the smaller mean
  if (mu[1] > mu[2]) {
    mu <- rev(mu); sigma <- rev(sigma); pi_k <- rev(pi_k)
  }

  list(mu = mu, sigma = sigma, pi = pi_k, iterations = iter)
}

gmm_boundary <- function(x) {
  fit <- gmm_em(x)
  # Find crossing point: where pi1*f1(x) = pi2*f2(x)
  # Search in a fine grid between the two means
  grid <- seq(min(x), max(x), length.out = 10000)
  w1 <- fit$pi[1] * dnorm(grid, fit$mu[1], fit$sigma[1])
  w2 <- fit$pi[2] * dnorm(grid, fit$mu[2], fit$sigma[2])
  diff_w <- w1 - w2
  # Find where the sign changes (crossing point between means)
  crossings <- which(diff(sign(diff_w)) != 0)
  # Pick the crossing closest to midpoint of the two means
  mid <- mean(fit$mu)
  if (length(crossings) == 0) {
    boundary <- mid
  } else {
    boundary <- grid[crossings[which.min(abs(grid[crossings] - mid))]]
  }
  list(boundary = boundary, fit = fit)
}

# ============================================================
# JOB 1: Labeled data (X, Y) — Drug activities
# ============================================================
cat("========================================\n")
cat("JOB 1: Labeled data (X, Y)\n")
cat("========================================\n\n")

x1 <- c(2, 1, 2, 5)
y1 <- c("A", "B", "A", "B")

cat("X =", x1, "\nY =", y1, "\n\n")

ig1 <- best_split(x1, y1, "IG")
kl1 <- best_split(x1, y1, "KL")

cat(sprintf("IG:   boundary = %.1f,  score = %.4f\n", ig1$boundary, ig1$score))
cat(sprintf("D_KL: boundary = %.1f,  score = %.4f\n", kl1$boundary, kl1$score))

# ============================================================
# JOB 2: Unlabeled data X — divide into two distributions
# ============================================================
cat("\n========================================\n")
cat("JOB 2: Unlabeled data X\n")
cat("========================================\n\n")

x2 <- c(3.2, 1.5, 4.8, 2.1, 6.3, 1.8, 5.5, 7.1, 2.9, 6.8)

cat("X =", x2, "\n\n")

# Step 1: Fit GMM to discover two groups
gmm_result <- gmm_boundary(x2)
cat(sprintf("GMM component 1: mu=%.2f, sigma=%.2f, pi=%.2f\n",
    gmm_result$fit$mu[1], gmm_result$fit$sigma[1], gmm_result$fit$pi[1]))
cat(sprintf("GMM component 2: mu=%.2f, sigma=%.2f, pi=%.2f\n",
    gmm_result$fit$mu[2], gmm_result$fit$sigma[2], gmm_result$fit$pi[2]))
cat(sprintf("GMM crossing-point boundary: %.4f\n\n", gmm_result$boundary))

# Step 2: Assign labels based on the GMM boundary
labels <- ifelse(x2 <= gmm_result$boundary, "G1", "G2")
cat("Assigned labels:", labels, "\n\n")

# Step 3: Score the boundary using IG and D_KL
ig2 <- best_split(x2, labels, "IG")
kl2 <- best_split(x2, labels, "KL")

cat(sprintf("IG:   boundary = %.4f,  score = %.4f\n", ig2$boundary, ig2$score))
cat(sprintf("D_KL: boundary = %.4f,  score = %.4f\n", kl2$boundary, kl2$score))

cat("\n========================================\n")
cat("NOTE: For Job 2, replace x2 with your actual data.\n")
cat("========================================\n")
