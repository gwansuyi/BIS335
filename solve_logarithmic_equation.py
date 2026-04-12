"""
Solve the equation: 2 * ln(m) * ln(m) = m   where m is a positive integer.
i.e., 2 * [ln(m)]^2 = m

This is a transcendental equation that cannot be solved in closed form
using elementary functions. We solve it numerically and then check
which positive integers satisfy (or most closely satisfy) the equation.

Analysis:
  Define f(m) = 2*[ln(m)]^2 - m for m > 0.

  As m -> 0+:  ln(m) -> -inf, so [ln(m)]^2 -> +inf, while m -> 0.
               Therefore f(m) -> +inf.

  At m = 1:    f(1) = 2*(0)^2 - 1 = -1 < 0.

  f has a local maximum at m ~ 8.6 (where f'(m) = 4*ln(m)/m - 1 = 0).

  As m -> +inf: m dominates [ln(m)]^2, so f(m) -> -inf.

  By the Intermediate Value Theorem, f has THREE positive real roots:
    Root 1: in (0, 1)    -- where f goes from +inf to -1
    Root 2: in (4, 5)    -- where f crosses from negative to positive
    Root 3: in (13, 14)  -- where f crosses from positive to negative

  Since m must be an integer, no exact solution exists (transcendental equation).
  We find the integer(s) that most closely satisfy the equation.
"""

import math


def f(m):
    """f(m) = 2*[ln(m)]^2 - m"""
    return 2 * math.log(m) ** 2 - m


def f_prime(m):
    """f'(m) = 4*ln(m)/m - 1"""
    return 4 * math.log(m) / m - 1


def bisection(a, b, tol=1e-15, max_iter=200):
    """Find root of f in [a, b] by bisection."""
    for _ in range(max_iter):
        mid = (a + b) / 2
        if f(mid) == 0 or (b - a) / 2 < tol:
            return mid
        if f(a) * f(mid) < 0:
            b = mid
        else:
            a = mid
    return (a + b) / 2


def newton(m0, tol=1e-15, max_iter=100):
    """Find root of f using Newton's method starting from m0."""
    m = m0
    for _ in range(max_iter):
        fm = f(m)
        fpm = f_prime(m)
        if abs(fpm) < 1e-30:
            break
        m_new = m - fm / fpm
        if m_new <= 0:
            m_new = m / 2  # safeguard: stay positive
        if abs(m_new - m) < tol:
            return m_new
        m = m_new
    return m


if __name__ == "__main__":
    print("=" * 60)
    print("Solving:  2 * [ln(m)]^2 = m   (m is a positive integer)")
    print("=" * 60)

    # --- Find all three continuous roots ---
    print("\n--- Continuous (real-valued) roots ---\n")

    root1 = newton(bisection(0.01, 0.99))
    root2 = newton(bisection(4.0, 5.0))
    root3 = newton(bisection(13.0, 14.0))

    for i, root in enumerate([root1, root2, root3], 1):
        lhs = 2 * math.log(root) ** 2
        print(f"  Root {i}:  m = {root:.10f}   "
              f"(verification: |LHS - RHS| = {abs(lhs - root):.2e})")

    # --- Check all positive integers from 1 to 30 ---
    print("\n--- Integer search: 2*[ln(m)]^2 vs m ---\n")
    print(f"  {'m':>4s}   {'2*[ln(m)]^2':>14s}   {'m':>6s}   {'difference':>12s}")
    print(f"  {'----':>4s}   {'--------------':>14s}   {'------':>6s}   {'----------':>12s}")

    best_integers = []
    for m in range(1, 31):
        lhs = 2 * math.log(m) ** 2
        diff = lhs - m
        marker = ""
        if abs(diff) < 0.5:
            marker = "  <-- close"
            best_integers.append((m, lhs, diff))
        print(f"  {m:4d}   {lhs:14.6f}   {m:6d}   {diff:+12.6f}{marker}")

    # --- Summary ---
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"\nThe three real-valued roots are:")
    print(f"  m1 = {root1:.10f}  (not an integer)")
    print(f"  m2 = {root2:.10f}  (nearest integers: {math.floor(root2)} and {math.ceil(root2)})")
    print(f"  m3 = {root3:.10f}  (nearest integers: {math.floor(root3)} and {math.ceil(root3)})")

    print(f"\nNo positive integer exactly satisfies 2*[ln(m)]^2 = m.")
    print(f"The closest integer solutions are:\n")
    for m, lhs, diff in best_integers:
        print(f"  m = {m:2d}:  2*[ln({m})]^2 = {lhs:.6f},  "
              f"difference = {diff:+.6f}")

    print("\n" + "=" * 60)
