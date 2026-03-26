# =============================================================================
# Binomial-Beta Relationship and Clopper-Pearson Confidence Interval
# =============================================================================

# --- Parameters ---
k <- 3      # number of successes
n <- 4      # number of trials
p <- 0.5    # probability of success
a <- 0.05   # significance level (alpha)

# =============================================================================
# Part 1: The Binomial-Beta Identity
# =============================================================================
# Key identity: pbinom(k, n, p) = pbeta(1 - p, n - k, k + 1)
#
# This connects the CDF of the Binomial to the CDF of the Beta distribution.
# Intuition: the incomplete beta function is the "tail sum" of binomial
# probabilities.

cat("=== Part 1: Binomial-Beta Identity ===\n")
cat(sprintf("pbinom(k=%d, n=%d, p=%.1f)     = %.6f\n", k, n, p, pbinom(k, n, p)))
cat(sprintf("pbeta(1-p, n-k, k+1)           = %.6f\n", pbeta(1 - p, n - k, k + 1)))
cat(sprintf("These are equal: %s\n\n",
            all.equal(pbinom(k, n, p), pbeta(1 - p, n - k, k + 1))))

# =============================================================================
# Part 2: Clopper-Pearson Exact Confidence Interval
# =============================================================================
# For k successes in n trials, the 100(1-a)% Clopper-Pearson CI for p is:
#
#   Lower = qbeta(a/2,   k,   n - k + 1)
#   Upper = qbeta(1-a/2, k+1, n - k)
#
# With a = 0.05:
#   Confidence Level = 1 - a = 0.95  (i.e., 95% CI)
#   Significance Level = a = 0.05

cat("=== Part 2: Clopper-Pearson Interval ===\n")

lower <- qbeta(a / 2, k, n - k + 1)
upper <- qbeta(1 - a / 2, k + 1, n - k)

cat(sprintf("k = %d, n = %d, alpha = %.2f\n", k, n, a))
cat(sprintf("Lower = qbeta(%.4f, %d, %d) = %.4f\n", a / 2, k, n - k + 1, lower))
cat(sprintf("Upper = qbeta(%.4f, %d, %d) = %.4f\n", 1 - a / 2, k + 1, n - k, upper))
cat(sprintf("95%% Clopper-Pearson CI for p: [%.4f, %.4f]\n\n", lower, upper))

cat(sprintf("Confidence Level  = 1 - alpha = 1 - %.2f = %.2f  (NOT 0.8)\n", a, 1 - a))
cat(sprintf("Significance Level = alpha     = %.2f             (NOT 0.1)\n\n", a))

# =============================================================================
# Part 3: WHY These Formulas Work (Derivation from Inverting the Binomial Test)
# =============================================================================
# The Clopper-Pearson interval inverts two one-sided binomial tests.
#
# --- Lower bound p_L ---
# Find p_L such that P(X >= k | p_L) = a/2:
#   1 - pbinom(k - 1, n, p_L) = a/2
#   pbinom(k - 1, n, p_L) = 1 - a/2
#
# Apply the binomial-beta identity: pbinom(k-1, n, p_L) = pbeta(1 - p_L, n-k+1, k)
#   pbeta(1 - p_L, n - k + 1, k) = 1 - a/2
#   1 - p_L = qbeta(1 - a/2, n - k + 1, k)
#   p_L     = 1 - qbeta(1 - a/2, n - k + 1, k)
#
# By beta symmetry: qbeta(p, a, b) = 1 - qbeta(1-p, b, a), so:
#   p_L = qbeta(a/2, k, n - k + 1)                              ... (*)
#
# --- Upper bound p_U ---
# Find p_U such that P(X <= k | p_U) = a/2:
#   pbinom(k, n, p_U) = a/2
#
# Apply the identity: pbinom(k, n, p_U) = pbeta(1 - p_U, n-k, k+1)
#   pbeta(1 - p_U, n - k, k + 1) = a/2
#   1 - p_U = qbeta(a/2, n - k, k + 1)
#   p_U     = 1 - qbeta(a/2, n - k, k + 1)
#
# By beta symmetry:
#   p_U = qbeta(1 - a/2, k + 1, n - k)                          ... (**)

cat("=== Part 3: Verification of the Derivation ===\n")

# Verify lower bound: P(X >= k | p_L) should equal a/2
prob_lower <- 1 - pbinom(k - 1, n, lower)
cat(sprintf("P(X >= %d | p_L = %.4f) = %.6f  (should be a/2 = %.4f)\n",
            k, lower, prob_lower, a / 2))

# Verify upper bound: P(X <= k | p_U) should equal a/2
prob_upper <- pbinom(k, n, upper)
cat(sprintf("P(X <= %d | p_U = %.4f) = %.6f  (should be a/2 = %.4f)\n",
            k, upper, prob_upper, a / 2))

# Verify beta symmetry property
cat(sprintf("\nBeta symmetry check:\n"))
cat(sprintf("  qbeta(a/2, k, n-k+1)           = %.6f\n", qbeta(a / 2, k, n - k + 1)))
cat(sprintf("  1 - qbeta(1-a/2, n-k+1, k)     = %.6f\n",
            1 - qbeta(1 - a / 2, n - k + 1, k)))
cat(sprintf("  Equal: %s\n\n",
            all.equal(qbeta(a / 2, k, n - k + 1),
                      1 - qbeta(1 - a / 2, n - k + 1, k))))

# =============================================================================
# Part 4: Summary of Key Relationships
# =============================================================================
cat("=== Part 4: Summary ===\n")
cat("
+------------------------------------------------------------------+
| Binomial-Beta Identity:                                          |
|   pbinom(k, n, p) = pbeta(1 - p, n - k, k + 1)                 |
+------------------------------------------------------------------+
| Clopper-Pearson 100(1-a)% CI for binomial proportion p:         |
|                                                                  |
|   Lower = qbeta(a/2,   k,   n - k + 1)   ... from inverting    |
|   Upper = qbeta(1-a/2, k+1, n - k)            binomial test +   |
|                                                beta symmetry     |
+------------------------------------------------------------------+
| Beta Symmetry:                                                   |
|   qbeta(p, a, b) = 1 - qbeta(1-p, b, a)                        |
+------------------------------------------------------------------+
| For a = 0.05:                                                    |
|   Confidence Level   = 1 - a = 0.95  (95%)                      |
|   Significance Level = a     = 0.05  (5%)                        |
+------------------------------------------------------------------+
")

# =============================================================================
# Part 5: Cross-check with binom.test()
# =============================================================================
cat("=== Part 5: Cross-check with binom.test() ===\n")
result <- binom.test(k, n, conf.level = 1 - a)
cat(sprintf("binom.test(%d, %d) 95%% CI: [%.4f, %.4f]\n",
            k, n, result$conf.int[1], result$conf.int[2]))
cat(sprintf("Our Clopper-Pearson CI:    [%.4f, %.4f]\n", lower, upper))
cat(sprintf("Match: %s\n",
            all.equal(as.numeric(result$conf.int), c(lower, upper))))
