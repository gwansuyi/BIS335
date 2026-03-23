# BIS335 - Optimal Boundary Finder using IG and D_KL
# For a single continuous variable X with class labels Y

# ── Helper functions ─────────────────────────────────────────────

shannon_entropy <- function(labels) {
  p <- table(labels) / length(labels)
  -sum(p * log2(p))
}

kl_divergence <- function(p, q) {
  # Add small epsilon to avoid log(0)
  eps <- 1e-10
  p <- p + eps
  q <- q + eps
  p <- p / sum(p)
  q <- q / sum(q)
  sum(p * log(p / q))
}

# ── IG-based optimal boundary ───────────────────────────────────

find_best_split_IG <- function(x, y) {
  sorted_idx <- order(x)
  x_sorted <- x[sorted_idx]
  y_sorted <- y[sorted_idx]
  n <- length(x)
  h_parent <- shannon_entropy(y)

  candidates <- which(diff(x_sorted) != 0)
  best_ig <- -Inf
  best_threshold <- NA

  for (i in candidates) {
    threshold <- (x_sorted[i] + x_sorted[i + 1]) / 2
    left_y <- y_sorted[1:i]
    right_y <- y_sorted[(i + 1):n]

    h_left <- shannon_entropy(left_y)
    h_right <- shannon_entropy(right_y)
    ig <- h_parent - (i / n) * h_left - ((n - i) / n) * h_right

    if (ig > best_ig) {
      best_ig <- ig
      best_threshold <- threshold
    }
  }

  list(threshold = best_threshold, score = best_ig)
}

# ── D_KL-based optimal boundary ─────────────────────────────────

find_best_split_KL <- function(x, y) {
  sorted_idx <- order(x)
  x_sorted <- x[sorted_idx]
  y_sorted <- y[sorted_idx]
  n <- length(x)

  # Get all unique class labels
  classes <- sort(unique(y))

  # Parent class distribution
  p_parent <- table(factor(y, levels = classes)) / n

  candidates <- which(diff(x_sorted) != 0)
  best_kl <- -Inf
  best_threshold <- NA

  for (i in candidates) {
    threshold <- (x_sorted[i] + x_sorted[i + 1]) / 2
    left_y <- y_sorted[1:i]
    right_y <- y_sorted[(i + 1):n]

    p_left <- table(factor(left_y, levels = classes)) / length(left_y)
    p_right <- table(factor(right_y, levels = classes)) / length(right_y)

    # Weighted symmetric KL: D_KL(left || parent) + D_KL(right || parent)
    dkl <- (i / n) * kl_divergence(p_left, p_parent) +
           ((n - i) / n) * kl_divergence(p_right, p_parent)

    if (dkl > best_kl) {
      best_kl <- dkl
      best_threshold <- threshold
    }
  }

  list(threshold = best_threshold, score = best_kl)
}

# ── Example: X = (2, 1, 2, 5) ───────────────────────────────────

x <- c(2, 1, 2, 5)
y <- c("A", "B", "A", "B")  # example class labels

cat("============================================\n")
cat("Optimal Boundary: X = (2, 1, 2, 5)\n")
cat("Class labels Y =", paste(y, collapse = ", "), "\n")
cat("============================================\n\n")

result_ig <- find_best_split_IG(x, y)
cat("IG-based split:\n")
cat("  Best boundary:", result_ig$threshold, "\n")
cat("  IG score:     ", round(result_ig$score, 4), "\n\n")

result_kl <- find_best_split_KL(x, y)
cat("D_KL-based split:\n")
cat("  Best boundary:", result_kl$threshold, "\n")
cat("  D_KL score:   ", round(result_kl$score, 4), "\n\n")

# ── Show all candidates for comparison ───────────────────────────

cat("--------------------------------------------\n")
cat("All candidate splits:\n")
cat("--------------------------------------------\n")

x_sorted <- sort(x)
candidates <- which(diff(x_sorted) != 0)

for (i in candidates) {
  threshold <- (x_sorted[i] + x_sorted[i + 1]) / 2

  sorted_idx <- order(x)
  y_sorted <- y[sorted_idx]
  n <- length(x)

  left_y <- y_sorted[1:i]
  right_y <- y_sorted[(i + 1):n]

  # IG
  h_parent <- shannon_entropy(y)
  h_left <- shannon_entropy(left_y)
  h_right <- shannon_entropy(right_y)
  ig <- h_parent - (i / n) * h_left - ((n - i) / n) * h_right

  # D_KL
  classes <- sort(unique(y))
  p_parent <- table(factor(y, levels = classes)) / n
  p_left <- table(factor(left_y, levels = classes)) / length(left_y)
  p_right <- table(factor(right_y, levels = classes)) / length(right_y)
  dkl <- (i / n) * kl_divergence(p_left, p_parent) +
         ((n - i) / n) * kl_divergence(p_right, p_parent)

  cat(sprintf("\n  Threshold %.1f: Left={%s} Right={%s}\n",
              threshold,
              paste(left_y, collapse = ","),
              paste(right_y, collapse = ",")))
  cat(sprintf("    IG = %.4f | D_KL = %.4f\n", ig, dkl))
}
