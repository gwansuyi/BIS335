"""
Find the smallest positive integer m satisfying: 2 * [ln(m)]^2 < m

Analysis:
  Define f(m) = 2*[ln(m)]^2 - m.  The inequality holds when f(m) < 0.

  f has three real roots at m ≈ 0.583, 4.428, and 13.706.
  Therefore f(m) < 0 in two intervals: (0, 0.583) and (0.583, 4.428)
  Wait -- more precisely:
    - f(m) > 0 for m in (0, 0.583)    [ln(m) large negative, squared dominates]
    - f(m) < 0 for m in (0.583, 4.428) [includes integers 1, 2, 3, 4]
    - f(m) > 0 for m in (4.428, 13.706) [includes integers 5..13]
    - f(m) < 0 for m in (13.706, +inf)  [includes integers 14, 15, 16, ...]

  So 2*[ln(m)]^2 < m holds for integers:  m in {1, 2, 3, 4} and m >= 14.
  The smallest positive integer satisfying the inequality is m = 1.
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
    print("Find smallest integer m where: 2*[ln(m)]^2 < m")
    print("=" * 60)

    # --- Find all three continuous roots (where 2*[ln(m)]^2 = m) ---
    print("\n--- Roots of 2*[ln(m)]^2 = m (boundary points) ---\n")

    root1 = newton(bisection(0.01, 0.99))
    root2 = newton(bisection(4.0, 5.0))
    root3 = newton(bisection(13.0, 14.0))

    for i, root in enumerate([root1, root2, root3], 1):
        lhs = 2 * math.log(root) ** 2
        print(f"  Root {i}:  m = {root:.10f}")

    # --- Check integers and classify ---
    print("\n--- Integer evaluation: 2*[ln(m)]^2 vs m ---\n")
    print(f"  {'m':>4s}   {'2*[ln(m)]^2':>14s}   {'m':>6s}   {'2*ln^2 < m?':>12s}")
    print(f"  {'----':>4s}   {'--------------':>14s}   {'------':>6s}   {'----------':>12s}")

    satisfies = []
    violates = []
    for m in range(1, 21):
        lhs = 2 * math.log(m) ** 2
        holds = lhs < m
        label = "YES" if holds else "no"
        print(f"  {m:4d}   {lhs:14.6f}   {m:6d}   {label:>12s}")
        if holds:
            satisfies.append(m)
        else:
            violates.append(m)

    # --- Answer ---
    print("\n" + "=" * 60)
    print("ANSWER")
    print("=" * 60)

    print(f"\nThe equation 2*[ln(m)]^2 = m has roots at:")
    print(f"  m ≈ {root1:.4f},  m ≈ {root2:.4f},  m ≈ {root3:.4f}")
    print(f"\nThe inequality 2*[ln(m)]^2 < m holds when:")
    print(f"  - m is in ({root1:.4f}, {root2:.4f})  => integers {{1, 2, 3, 4}}")
    print(f"  - m is in ({root3:.4f}, +∞)       => integers {{14, 15, 16, ...}}")
    print(f"\nThe inequality FAILS for m in {{5, 6, 7, 8, 9, 10, 11, 12, 13}}.")
    print(f"\n>>> The smallest positive integer m satisfying 2*[ln(m)]^2 < m is:  m = {satisfies[0]}")
    print("=" * 60)
