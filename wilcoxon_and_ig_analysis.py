"""
BIS335 - Wilcoxon Test & Information Gain Analysis
"""

import math
from collections import Counter


# ── Information Gain ────────────────────────────────────────────────

def entropy(values):
    """Compute Shannon entropy H(X) = -Σ p(x) log2 p(x)."""
    n = len(values)
    counts = Counter(values)
    return -sum((c / n) * math.log2(c / n) for c in counts.values())


def information_gain(parent, left, right):
    """IG = H(parent) - weighted avg of children entropies."""
    n = len(parent)
    h_parent = entropy(parent)
    h_children = (len(left) / n) * entropy(left) + (len(right) / n) * entropy(right)
    return h_parent - h_children


# Given: X = (2, 1, 2, 5)
X = [2, 1, 2, 5]

print("=" * 50)
print("Information Gain Analysis for X = (2, 1, 2, 5)")
print("=" * 50)

# Parent entropy
h = entropy(X)
print(f"\nValue distribution: {dict(Counter(X))}")
print(f"Parent Entropy H(X) = {h:.4f} bits")

# Show step-by-step
print("\nStep-by-step:")
n = len(X)
for val, count in sorted(Counter(X).items()):
    p = count / n
    contrib = -p * math.log2(p)
    print(f"  Value {val}: p={count}/{n}={p:.2f}, "
          f"-p*log2(p) = -{p:.2f}*{math.log2(p):.2f} = {contrib:.4f}")

print(f"\nH(X) = {h:.4f} bits")

# Possible binary splits (sorted: 1, 2, 2, 5)
print("\n" + "-" * 50)
print("Possible splits (sorted X = [1, 2, 2, 5]):")
print("-" * 50)

sorted_x = sorted(X)
split_points = []
for i in range(len(sorted_x) - 1):
    if sorted_x[i] != sorted_x[i + 1]:
        midpoint = (sorted_x[i] + sorted_x[i + 1]) / 2
        split_points.append((midpoint, i + 1))

for threshold, idx in split_points:
    left = sorted_x[:idx]
    right = sorted_x[idx:]
    ig = information_gain(sorted_x, left, right)
    h_left = entropy(left) if len(set(left)) > 1 else 0.0
    h_right = entropy(right) if len(set(right)) > 1 else 0.0
    print(f"\n  Split at {threshold}:")
    print(f"    Left  = {left}, H = {h_left:.4f}")
    print(f"    Right = {right}, H = {h_right:.4f}")
    print(f"    IG = {h:.4f} - ({len(left)}/{n})*{h_left:.4f} - ({len(right)}/{n})*{h_right:.4f} = {ig:.4f}")


# ── Wilcoxon Signed-Rank Test ──────────────────────────────────────

print("\n" + "=" * 50)
print("Wilcoxon Signed-Rank Test (example)")
print("=" * 50)


def wilcoxon_signed_rank(x, y):
    """Manual Wilcoxon signed-rank test for paired samples."""
    diffs = [xi - yi for xi, yi in zip(x, y)]
    abs_diffs = [(abs(d), 1 if d > 0 else -1, i)
                 for i, d in enumerate(diffs) if d != 0]
    abs_diffs.sort(key=lambda t: t[0])

    # Assign ranks
    ranks = list(range(1, len(abs_diffs) + 1))

    # Handle ties by averaging ranks
    i = 0
    while i < len(abs_diffs):
        j = i
        while j < len(abs_diffs) and abs_diffs[j][0] == abs_diffs[i][0]:
            j += 1
        avg_rank = sum(ranks[i:j]) / (j - i)
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j

    w_plus = sum(r for r, (_, sign, _) in zip(ranks, abs_diffs) if sign > 0)
    w_minus = sum(r for r, (_, sign, _) in zip(ranks, abs_diffs) if sign < 0)

    print(f"\n  Differences: {diffs}")
    print(f"  |Differences| with signs and ranks:")
    for rank, (abs_d, sign, idx) in zip(ranks, abs_diffs):
        print(f"    Pair {idx+1}: diff={diffs[idx]:+d}, |diff|={abs_d}, rank={rank}")
    print(f"\n  W+ (positive ranks) = {w_plus}")
    print(f"  W- (negative ranks) = {w_minus}")
    print(f"  Test statistic W = min(W+, W-) = {min(w_plus, w_minus)}")

    return min(w_plus, w_minus)


# Example paired data
sample_x = [2, 1, 2, 5]
sample_y = [1, 3, 1, 4]
W = wilcoxon_signed_rank(sample_x, sample_y)
