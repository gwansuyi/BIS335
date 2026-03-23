# Optimal split boundary using IG and D_KL for X=(2,1,2,5), Y=(A,B,A,B)

H <- function(y) { p <- table(y)/length(y); -sum(p * log2(p)) }

DKL <- function(p, q) { e <- 1e-10; p <- p/sum(p)+e; q <- q/sum(q)+e; sum(p * log(p/q)) }

best_split <- function(x, y, method = c("IG", "KL")) {
  method <- match.arg(method)
  idx <- order(x); xs <- x[idx]; ys <- y[idx]; n <- length(x)
  cls <- sort(unique(y)); p0 <- table(factor(y, levels=cls)) / n
  cuts <- which(diff(xs) != 0)
  scores <- sapply(cuts, function(i) {
    L <- ys[1:i]; R <- ys[(i+1):n]
    if (method == "IG") {
      H(y) - (i/n)*H(L) - ((n-i)/n)*H(R)
    } else {
      pL <- table(factor(L, levels=cls))/length(L)
      pR <- table(factor(R, levels=cls))/length(R)
      (i/n)*DKL(pL, p0) + ((n-i)/n)*DKL(pR, p0)
    }
  })
  best <- which.max(scores)
  list(boundary = (xs[cuts[best]] + xs[cuts[best]+1]) / 2, score = scores[best])
}

# ── Run ──
x <- c(2, 1, 2, 5)
y <- c("A", "B", "A", "B")

cat("X =", x, " Y =", y, "\n\n")

ig  <- best_split(x, y, "IG")
kl  <- best_split(x, y, "KL")

cat(sprintf("IG:   boundary = %.1f, score = %.4f\n", ig$boundary, ig$score))
cat(sprintf("D_KL: boundary = %.1f, score = %.4f\n", kl$boundary, kl$score))
